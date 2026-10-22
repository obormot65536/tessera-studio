"""CalculiX (ccx) adapter: writes an ``.inp`` deck, runs ``ccx``, reads ``.frd`` and ``.dat``.

ADR-004: the adapter receives a discrete, solver-neutral ``SolveJob`` and returns a
``SolveResult``. Units are mm, N, MPa (ADR-002), so no conversion happens here.

Numbering: Tessera uses 0-based node and element indices; CalculiX uses 1-based ids. Node ``i``
is written as ccx node ``i + 1`` and element ``j`` as ccx element ``j + 1``. The model's tet10
connectivity is already in C3D10 order, so it is written unchanged (plus one).

Precision: the ASCII ``.frd`` file stores results as ``E12.5``, i.e. six significant digits.
Results read back from it carry a relative rounding error of up to 5e-6.

Integration points (ADR-008): ``*EL PRINT`` writes every element's stresses at its integration
points, and their coordinates, to ``job.dat`` (``E13.6``, seven significant digits). These are
ccx's own values, not extrapolated or averaged. ``.dat`` prints sxx syy szz sxy sxz syz, so its
stress columns ``[0, 1, 2, 3, 5, 4]`` are the neutral order xx yy zz xy yz zx.
"""

import os
import re
import shutil
import subprocess
from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

from tessera.model.mesh import Mesh
from tessera.model.solver import (
    AnalysisKind,
    PointStresses,
    SolveJob,
    SolverError,
    SolveResult,
)

JOB_NAME = "job"
CCX_ENV_VAR = "TESSERA_CCX"

# CalculiX names: at most 80 characters, case-insensitive (ccx upper-cases them), and no commas,
# spaces, '=' or ':' since those are keyword-line separators or special.
_UNSAFE_RE = re.compile(r"[^A-Za-z0-9_]")
_MAX_NAME = 80
_IDS_PER_LINE = 16  # ccx accepts at most 16 entries per data line in *NSET/*ELSET

ALL_ELSET = "EALL"
"""The set of every element, defined on ``*ELEMENT`` and printed to ``job.dat``. User sets are
always prefixed (``E_``), so none can extend it: ccx would print its elements twice."""
POINT_RULE = "tet4"
"""C3D10's 4-point rule; point k is the one nearest corner node k (ADR-008, finding 4)."""
_RULE_POINTS = 4


# --------------------------------------------------------------------------------------------
# Deck writer
# --------------------------------------------------------------------------------------------


def _safe_names(names: Iterable[str], prefix: str) -> dict[str, str]:
    """Map neutral names (e.g. ``part:beam``) to unique, deterministic CalculiX names.

    Unsafe characters become ``_``, names get ``prefix`` (so they never start with a digit and
    element sets, node sets and materials stay visibly apart), are cut to 80 characters and,
    when two names collide after upper-casing, get a ``_2``, ``_3``... suffix in input order.
    """
    mapping: dict[str, str] = {}
    used: set[str] = set()
    for name in names:
        if name in mapping:
            raise SolverError(f"Duplicate name {name!r}")
        base = (prefix + _UNSAFE_RE.sub("_", name))[:_MAX_NAME]
        candidate, n = base, 1
        while candidate.upper() in used:
            n += 1
            suffix = f"_{n}"
            candidate = base[: _MAX_NAME - len(suffix)] + suffix
        used.add(candidate.upper())
        mapping[name] = candidate
    return mapping


@dataclass(frozen=True)
class DeckNames:
    """Neutral name → CalculiX name, per namespace."""

    element_sets: dict[str, str]
    node_sets: dict[str, str]
    materials: dict[str, str]

    def as_metadata(self) -> dict[str, str]:
        return {
            f"ccx_name.{kind}.{neutral}": ccx
            for kind, table in (
                ("elset", self.element_sets),
                ("nset", self.node_sets),
                ("material", self.materials),
            )
            for neutral, ccx in table.items()
        }


def deck_names(job: SolveJob) -> DeckNames:
    """The CalculiX names ``write_deck`` uses for ``job``'s sets and materials."""
    return DeckNames(
        element_sets=_safe_names(job.mesh.element_sets, "E_"),
        node_sets=_safe_names(job.mesh.node_sets, "N_"),
        materials=_safe_names((m.name for m in job.materials), "M_"),
    )


def _num(x: float) -> str:
    # CalculiX numeric fields hold at most 20 characters; 13 significant digits in
    # scientific notation always fit, including a sign and a three-digit exponent.
    return f"{float(x):.12e}"


def _id_lines(ids: Iterable[int]) -> list[str]:
    ids = list(ids)
    return [
        ", ".join(str(i) for i in ids[k : k + _IDS_PER_LINE])
        for k in range(0, len(ids), _IDS_PER_LINE)
    ]


def _check_indices(kind: str, name: str, idx: NDArray[np.int64], upper: int) -> None:
    if idx.size and (int(idx.min()) < 0 or int(idx.max()) >= upper):
        raise SolverError(f"{kind} {name!r} has indices outside [0, {upper})")


def write_deck(job: SolveJob) -> str:
    """Return the CalculiX input deck (``.inp`` text) for ``job``.

    Neutral set and material names are mapped to CalculiX-safe names (see ``deck_names``);
    the mapping is written as comments at the top of the deck. Raises ``SolverError`` when the
    job cannot be expressed as a valid deck (duplicate material names, unknown sets or
    materials, out-of-range indices, unsupported analysis).
    """
    if job.analysis is not AnalysisKind.STATIC:
        raise SolverError(f"CalculiX adapter does not support analysis {job.analysis!r}")
    mesh = job.mesh
    n_nodes, n_elems = mesh.n_nodes, mesh.n_elements
    _check_indices("element connectivity", "elements", mesh.elements, n_nodes)

    for name, elems in mesh.element_sets.items():
        _check_indices("Element set", name, np.asarray(elems), n_elems)
    for name, nodes in mesh.node_sets.items():
        _check_indices("Node set", name, np.asarray(nodes), n_nodes)
    names = deck_names(job)
    elsets, nsets, mats = names.element_sets, names.node_sets, names.materials

    out: list[str] = ["** Tessera Studio: generated CalculiX deck"]
    out.extend(
        f"** {kind} {neutral!r} -> {ccx}"
        for kind, table in (("ELSET", elsets), ("NSET", nsets), ("MATERIAL", mats))
        for neutral, ccx in table.items()
    )
    out.append("*NODE")
    out.extend(
        f"{i + 1}, {_num(x)}, {_num(y)}, {_num(z)}" for i, (x, y, z) in enumerate(mesh.nodes)
    )
    out.append(f"*ELEMENT, TYPE=C3D10, ELSET={ALL_ELSET}")
    out.extend(
        f"{j + 1}, " + ", ".join(str(int(n) + 1) for n in conn)
        for j, conn in enumerate(mesh.elements)
    )
    for name, elems in mesh.element_sets.items():
        out.append(f"*ELSET, ELSET={elsets[name]}")
        out.extend(_id_lines(int(e) + 1 for e in np.asarray(elems).ravel()))
    for name, nodes in mesh.node_sets.items():
        out.append(f"*NSET, NSET={nsets[name]}")
        out.extend(_id_lines(int(n) + 1 for n in np.asarray(nodes).ravel()))

    for mat in job.materials:
        out.append(f"*MATERIAL, NAME={mats[mat.name]}")
        out.append("*ELASTIC")
        out.append(f"{_num(mat.youngs_modulus)}, {_num(mat.poisson_ratio)}")
        if mat.density is not None:
            out.append("*DENSITY")
            out.append(_num(mat.density))

    covered = np.zeros(n_elems, dtype=bool)
    for sec in job.sections:
        if sec.element_set not in mesh.element_sets:
            raise SolverError(f"Section refers to unknown element set {sec.element_set!r}")
        if sec.material not in mats:
            raise SolverError(f"Section refers to unknown material {sec.material!r}")
        covered[np.asarray(mesh.element_sets[sec.element_set])] = True
        out.append(
            f"*SOLID SECTION, ELSET={elsets[sec.element_set]}, MATERIAL={mats[sec.material]}"
        )
    if not covered.all():
        raise SolverError(
            f"{int((~covered).sum())} element(s) are not assigned a section (material)"
        )

    out.extend(["*STEP", "*STATIC"])
    if job.supports:
        out.append("*BOUNDARY")
        for sup in job.supports:
            if sup.node_set not in mesh.node_sets:
                raise SolverError(f"Support refers to unknown node set {sup.node_set!r}")
            for dof in sup.dofs:
                if dof not in (1, 2, 3):
                    raise SolverError(f"Support on {sup.node_set!r}: invalid DOF {dof}")
                out.append(f"{nsets[sup.node_set]}, {dof}, {dof}")

    totals = np.zeros((n_nodes, 3), dtype=np.float64)
    for load in job.loads:
        nodes = np.asarray(load.nodes, dtype=np.int64).ravel()
        forces = np.asarray(load.forces, dtype=np.float64).reshape(-1, 3)
        if nodes.shape[0] != forces.shape[0]:
            raise SolverError("NodalLoads: nodes and forces have different lengths")
        _check_indices("NodalLoads", "nodes", nodes, n_nodes)
        np.add.at(totals, nodes, forces)
    cloads = [
        f"{n + 1}, {d + 1}, {_num(totals[n, d])}"
        for n in range(n_nodes)
        for d in range(3)
        if totals[n, d] != 0.0
    ]
    if cloads:
        out.append("*CLOAD")
        out.extend(cloads)

    out.extend(["*NODE FILE", "U", "*EL FILE", "S"])
    out.extend([f"*EL PRINT, ELSET={ALL_ELSET}", "S", f"*EL PRINT, ELSET={ALL_ELSET}", "COORD"])
    out.extend(["*END STEP", ""])
    return "\n".join(out)


# --------------------------------------------------------------------------------------------
# FRD parser
# --------------------------------------------------------------------------------------------


@dataclass
class FrdResults:
    """Nodal result blocks of the last step in an ASCII ``.frd`` file, keyed by ccx node id."""

    blocks: dict[str, dict[int, list[float]]] = field(default_factory=dict)
    """Block name (e.g. ``DISP``, ``STRESS``) → {node id (1-based) → values}."""
    version: str | None = None


def _split_values(text: str) -> list[float]:
    """Split fixed-width ``E12.5`` fields (they may run together, e.g. ``1.0E+00-2.0E+00``)."""
    text = text.rstrip()
    return [float(text[k : k + 12]) for k in range(0, len(text), 12)]


def parse_frd(text: str) -> FrdResults:
    """Parse the nodal result blocks of an ASCII CalculiX ``.frd`` file.

    Later blocks with the same name overwrite earlier ones, so the result holds the last step /
    increment of each block.
    """
    result = FrdResults()
    lines = text.splitlines()
    i = 0
    id_width = 10
    while i < len(lines):
        line = lines[i]
        if line.startswith("    1UVERSION"):
            result.version = line[13:].strip() or None
        elif line.startswith("  100C"):
            # Last field of the result header is the format: 0 = short (I5 ids), 1 = long (I10).
            tokens = line.split()
            id_width = 5 if tokens and tokens[-1] == "0" else 10
        elif line.startswith(" -4"):
            name = line[3:].split()[0]
            values: dict[int, list[float]] = {}
            last: list[float] | None = None
            i += 1
            while i < len(lines) and not lines[i].startswith(" -3"):
                rec = lines[i]
                if rec.startswith(" -1"):
                    node = int(rec[3 : 3 + id_width])
                    last = _split_values(rec[3 + id_width :])
                    values[node] = last
                elif rec.startswith(" -2") and last is not None:
                    last.extend(_split_values(rec[3 + id_width :]))
                i += 1
            result.blocks[name] = values
        i += 1
    return result


def _to_array(block: dict[int, list[float]], n_nodes: int, n_comp: int) -> NDArray[np.float64]:
    """Map a ccx block onto 0-based node order; nodes absent from the file get NaN."""
    arr = np.full((n_nodes, n_comp), np.nan, dtype=np.float64)
    for node_id, vals in block.items():
        idx = node_id - 1
        if 0 <= idx < n_nodes:
            arr[idx, :] = vals[:n_comp]
    return arr


def check_results_complete(
    mesh: Mesh,
    displacements: NDArray[np.float64],
    stresses: NDArray[np.float64],
    point_stresses: PointStresses | None = None,
    log: str = "",
) -> None:
    """Raise if a node used by an element has no result (orphan nodes may stay NaN), or, when
    ``point_stresses`` are given, if an element lacks an integration-point result."""
    used = np.unique(mesh.elements)
    bad = used[np.isnan(displacements[used]).any(axis=1) | np.isnan(stresses[used]).any(axis=1)]
    if bad.size:
        raise SolverError(
            f"CalculiX returned no result for {bad.size} mesh node(s), e.g. node {bad[0]}", log
        )
    if point_stresses is None:
        return
    if point_stresses.n_elements != mesh.n_elements:
        raise SolverError(
            f"CalculiX returned integration-point results for {point_stresses.n_elements} "
            f"element(s); the mesh has {mesh.n_elements} element(s)",
            log,
        )
    missing = np.isnan(point_stresses.points).any(axis=(1, 2))
    missing |= np.isnan(point_stresses.stresses).any(axis=(1, 2))
    bad = np.flatnonzero(missing)
    if bad.size:
        raise SolverError(
            f"CalculiX returned no integration-point result for {bad.size} element(s), "
            f"e.g. element {bad[0]}",
            log,
        )


# --------------------------------------------------------------------------------------------
# DAT parser (integration points, ADR-008)
# --------------------------------------------------------------------------------------------

DAT_STRESSES = "stresses"
DAT_COORDINATES = "global coordinates"
_DAT_STRESS_COLUMNS = "sxx,syy,szz,sxy,sxz,syz"
_DAT_TO_NEUTRAL = [0, 1, 2, 3, 5, 4]
_DAT_HEADER = re.compile(
    r"^ (?P<name>[a-z][a-z ]*?) \(elem, integ\.pnt\.,(?P<columns>[^)]*)\) for set \S+ and time.*$",
    re.MULTILINE,
)
_BLANK_LINE = re.compile(r"\n[ \t]*\n")


@dataclass(frozen=True)
class DatBlock:
    """One ``*EL PRINT`` block of ``job.dat``: a row per element and integration point."""

    elements: NDArray[np.int64]
    """(K,) ccx element ids (1-based), ascending."""
    values: NDArray[np.float64]
    """(K, Q, C) values at points 1..Q of each element; stresses in the neutral order."""


@dataclass
class DatResults:
    """Integration-point blocks of a ``job.dat`` file, by name (``stresses``,
    ``global coordinates``); the last block of each name, i.e. the last increment."""

    blocks: dict[str, DatBlock] = field(default_factory=dict)


def parse_dat(text: str) -> DatResults:
    """Parse the integration-point blocks of a CalculiX ``.dat`` file.

    Each block has a header line, a blank line and one row per element and point: element id,
    point number, then the values, all separated by spaces. A later block of the same name
    replaces an earlier one (ccx writes a block per increment and step). Stresses are reordered
    from ``.dat``'s sxx syy szz sxy sxz syz to the neutral xx yy zz xy yz zx. Raises
    ``SolverError`` for a malformed row, a duplicate row, or an element that does not have
    exactly points 1..Q, Q being the same for every element of the block.
    """
    result = DatResults()
    headers = list(_DAT_HEADER.finditer(text))
    for k, header in enumerate(headers):
        name, columns = header["name"], header["columns"]
        end = headers[k + 1].start() if k + 1 < len(headers) else len(text)
        rows = _BLANK_LINE.split(text[header.end() : end].strip(), maxsplit=1)[0]
        n_rows = rows.count("\n") + 1 if rows else 0
        width = 2 + len(columns.split(","))
        tokens = rows.split()
        if name == DAT_STRESSES and columns != _DAT_STRESS_COLUMNS:
            raise SolverError(f"CalculiX job.dat: unexpected stress columns {columns!r}")
        try:
            if len(tokens) != n_rows * width:
                raise ValueError(f"expected {width} numbers per row")
            table = np.array(tokens, dtype=np.float64).reshape(n_rows, width)
        except ValueError as exc:
            raise SolverError(f"CalculiX job.dat block {name!r} is malformed: {exc}") from None
        values = _point_block(name, table)
        if name == DAT_STRESSES:
            values = DatBlock(values.elements, values.values[..., _DAT_TO_NEUTRAL])
        result.blocks[name] = values
    return result


def _point_block(name: str, table: NDArray[np.float64]) -> DatBlock:
    """Group the rows (element, point, values...) by element, checking points 1..Q."""
    ids, points = table[:, 0].astype(np.int64), table[:, 1].astype(np.int64)
    if not (np.array_equal(ids, table[:, 0]) and np.array_equal(points, table[:, 1])):
        raise SolverError(f"CalculiX job.dat block {name!r} is malformed: non-integer ids")
    q = int(points.max()) if points.size else 0
    keys = np.sort(ids * (q + 1) + points)
    repeated = keys[1:][keys[1:] == keys[:-1]]
    if repeated.size:
        element, point = divmod(int(repeated[0]), q + 1)
        raise SolverError(
            f"CalculiX job.dat block {name!r} has duplicate rows, e.g. element {element} "
            f"point {point}"
        )
    elements, counts = np.unique(ids, return_counts=True)
    if points.size and (points.min() < 1 or (counts != q).any()):
        bad = elements[counts != q][0] if (counts != q).any() else ids[points.argmin()]
        raise SolverError(
            f"CalculiX job.dat block {name!r}: every element needs points 1 to {q} once; "
            f"element {bad} does not have them"
        )
    order = np.lexsort((points, ids))
    values = table[order, 2:].reshape(len(elements), q, table.shape[1] - 2)
    return DatBlock(elements=elements, values=values)


def point_stresses(dat: DatResults, n_elements: int, n_points: int = _RULE_POINTS) -> PointStresses:
    """Integration-point stresses in mesh element order; elements ccx did not print are NaN."""
    arrays: list[NDArray[np.float64]] = []
    for name, width in ((DAT_COORDINATES, 3), (DAT_STRESSES, 6)):
        block = dat.blocks.get(name)
        if block is None:
            raise SolverError(f"CalculiX job.dat lacks the {name!r} block")
        if block.values.shape[1:] != (n_points, width):
            raise SolverError(
                f"CalculiX job.dat block {name!r} has {block.values.shape[1]} integration "
                f"points per element with {block.values.shape[2]} values; expected {n_points} "
                f"with {width}"
            )
        index = block.elements - 1
        if index.size and (index.min() < 0 or index.max() >= n_elements):
            raise SolverError(
                f"CalculiX job.dat block {name!r} has element ids outside 1..{n_elements}"
            )
        array = np.full((n_elements, n_points, width), np.nan)
        array[index] = block.values
        arrays.append(array)
    return PointStresses(points=arrays[0], stresses=arrays[1])


def read_point_stresses(path: Path, n_elements: int, log: str = "") -> PointStresses:
    """Read ``job.dat``. A missing or empty file is a ``SolverError``: ccx exits 0 with an empty
    ``.dat`` when ``*EL PRINT`` names no set (ADR-008, finding 7)."""
    if not path.is_file() or path.stat().st_size == 0:
        raise SolverError(
            f"CalculiX wrote no integration-point results ({path.name} is missing or empty)", log
        )
    try:
        return point_stresses(parse_dat(path.read_text()), n_elements)
    except SolverError as exc:
        raise SolverError(str(exc), log) from None


# --------------------------------------------------------------------------------------------
# Adapter
# --------------------------------------------------------------------------------------------


def find_ccx() -> str:
    """Return the ccx executable: ``$TESSERA_CCX`` if set, else ``ccx`` on PATH."""
    configured = os.environ.get(CCX_ENV_VAR)
    if configured:
        found = shutil.which(configured)
        if found is None:
            raise SolverError(
                f"{CCX_ENV_VAR}={configured!r} does not point to an executable CalculiX binary"
            )
        return found
    found = shutil.which("ccx")
    if found is None:
        raise SolverError(
            "CalculiX executable 'ccx' not found on PATH; install CalculiX 2.21+ "
            f"or set {CCX_ENV_VAR} to the ccx binary"
        )
    return found


class CalculixAdapter:
    """``SolverAdapter`` for linear static analysis with CalculiX."""

    name = "calculix"

    def __init__(self, timeout: float | None = None) -> None:
        self.timeout = timeout

    def capabilities(self) -> frozenset[AnalysisKind]:
        return frozenset({AnalysisKind.STATIC})

    def run(self, job: SolveJob, workdir: Path) -> SolveResult:
        deck = write_deck(job)
        ccx = find_ccx()
        workdir = Path(workdir)
        workdir.mkdir(parents=True, exist_ok=True)
        for suffix in (".frd", ".dat", ".sta", ".cvg"):
            (workdir / f"{JOB_NAME}{suffix}").unlink(missing_ok=True)
        (workdir / f"{JOB_NAME}.inp").write_text(deck)

        try:
            proc = subprocess.run(
                [ccx, "-i", JOB_NAME],
                cwd=workdir,
                capture_output=True,
                text=True,
                timeout=self.timeout,
                check=False,
            )
        except subprocess.TimeoutExpired as exc:
            log = _as_text(exc.stdout) + _as_text(exc.stderr)
            raise SolverError(f"CalculiX timed out after {self.timeout} s", log) from exc
        except OSError as exc:
            raise SolverError(f"Could not start CalculiX ({ccx}): {exc}") from exc

        log = proc.stdout + (f"\n--- stderr ---\n{proc.stderr}" if proc.stderr else "")
        frd = workdir / f"{JOB_NAME}.frd"
        if "*ERROR" in proc.stdout or "*ERROR" in proc.stderr:
            raise SolverError(f"CalculiX reported an error: {_first_error(log)}", log)
        if proc.returncode != 0:
            raise SolverError(f"CalculiX exited with code {proc.returncode}", log)
        if not frd.is_file() or frd.stat().st_size == 0:
            raise SolverError("CalculiX produced no results file (job.frd)", log)

        parsed = parse_frd(frd.read_text())
        n = job.mesh.n_nodes
        missing = [b for b in ("DISP", "STRESS") if b not in parsed.blocks]
        if missing:
            raise SolverError(f"CalculiX results file lacks block(s) {missing}", log)
        metadata = {"workdir": str(workdir), **deck_names(job).as_metadata()}
        if parsed.version:
            metadata["ccx_version"] = parsed.version
        displacements = _to_array(parsed.blocks["DISP"], n, 3)
        stresses = _to_array(parsed.blocks["STRESS"], n, 6)
        points = read_point_stresses(workdir / f"{JOB_NAME}.dat", job.mesh.n_elements, log)
        check_results_complete(job.mesh, displacements, stresses, points, log)
        metadata["point_rule"] = POINT_RULE
        return SolveResult(
            displacements=displacements,
            stresses=stresses,
            solver=self.name,
            log=log,
            metadata=metadata,
            point_stresses=points,
        )


def _as_text(data: str | bytes | None) -> str:
    if data is None:
        return ""
    return data.decode(errors="replace") if isinstance(data, bytes) else data


def _first_error(log: str) -> str:
    lines: Sequence[str] = log.splitlines()
    for k, line in enumerate(lines):
        if "*ERROR" in line:
            return " ".join(s.strip() for s in lines[k : k + 3] if s.strip())
    return "see log"
