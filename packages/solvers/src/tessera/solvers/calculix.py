"""CalculiX (ccx) adapter: writes an ``.inp`` deck, runs ``ccx``, reads ``.frd`` and ``.dat``.

ADR-004: the adapter receives a discrete, solver-neutral ``SolveJob`` and returns a
``SolveResult``. Units are mm, N, MPa (ADR-002), so no conversion happens here.

Numbering: Tessera uses 0-based node and element indices; CalculiX uses 1-based ids. Node ``i``
is written as ccx node ``i + 1`` and element ``j`` as ccx element ``j + 1``. The model's tet10
connectivity is already in C3D10 order, so it is written unchanged (plus one).

Precision: the ASCII ``.frd`` file stores results as ``E12.5``, i.e. six significant digits.
Results read back from it carry a relative rounding error of up to 5e-6.

Orthotropic materials (ADR-009 D2): ``*ELASTIC, TYPE=ENGINEERING CONSTANTS`` with E1 E2 E3 nu12
nu13 nu23 G12 G13 on one line and G23 on the next (major Poisson ratios, as in the model), and
one ``*ORIENTATION, SYSTEM=RECTANGULAR`` (a = material axis 1, b = material axis 2, in part
coordinates) per section with a ``MaterialFrame``, referenced by ``ORIENTATION=`` on its
``*SOLID SECTION``. Verified on ccx 2.21 (``test_calculix_orthotropic.py``):

1. ``.frd`` stresses (and displacements) are in GLOBAL axes with an ``*ORIENTATION`` active: a
   cube with material axis 1 along global -z, pulled along z, returns SZZ = 9.99990 and the other
   components below 1e-12 at every node; local axes would put sigma_11 in the SXX slot.
2. The card order ``*MATERIAL``, ``*ORIENTATION``, ``*SOLID SECTION`` is accepted. ccx 2.21 also
   gives identical results with ``*ORIENTATION`` before ``*MATERIAL`` or after the
   ``*SOLID SECTION`` that names it; the order written here is the one the ADR specifies.

Integration points (ADR-008), only when ``SolveJob.point_stresses`` asks for them:
``*EL PRINT`` writes every element's stresses at its integration points to ``job.dat``
(``E13.6``, seven significant digits). These are ccx's own values, not extrapolated or averaged.
``.dat`` prints sxx syy szz sxy sxz syz, so its stress columns ``[0, 1, 2, 3, 5, 4]`` are the
neutral order xx yy zz xy yz zx. The points' coordinates are computed from the mesh with the
C3D10 rule (``rule_points``), in full double precision: ccx's ``COORD`` output has seven
significant digits of the absolute position, too few for a small element far from the origin.

Stopping: ccx runs in its own process group (POSIX) or is stopped with ``taskkill /T``
(Windows), so a timeout or an interrupt also stops the processes it, or a wrapper script set in
``$TESSERA_CCX``, started. See ``_run_ccx`` for what is and is not covered.
"""

import contextlib
import ctypes
import math
import mmap
import os
import re
import shutil
import signal
import subprocess
import sys
import threading
import warnings
from collections.abc import Callable, Iterable, Iterator, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import IO, Any

import numpy as np
from numpy.typing import NDArray

from tessera.model.mesh import TET10_EDGES, Mesh
from tessera.model.solver import (
    AnalysisKind,
    OrthotropicMaterial,
    PointStresses,
    SolveJob,
    SolverError,
    SolveResult,
    SolverTimeoutError,
)

JOB_NAME = "job"
CCX_ENV_VAR = "TESSERA_CCX"

# CalculiX names: at most 80 characters, case-insensitive (ccx upper-cases them), and no commas,
# spaces, '=' or ':' since those are keyword-line separators or special.
_UNSAFE_RE = re.compile(r"[^A-Za-z0-9_]")
_MAX_NAME = 80
_IDS_PER_LINE = 16  # ccx accepts at most 16 entries per data line in *NSET/*ELSET

ALL_ELSET = "EALL"
"""The set of every element, defined on ``*ELEMENT`` and printed to ``job.dat`` when point
stresses are asked for. User sets are always prefixed (``E_``), so none can extend it: ccx
would print its elements twice."""
POINT_RULE = "tet4"
"""C3D10's 4-point rule; point k is the one nearest corner node k (ADR-008, finding 4)."""
_RULE_A = (5.0 - np.sqrt(5.0)) / 20.0
_RULE_B = (5.0 + 3.0 * np.sqrt(5.0)) / 20.0
_RULE_NATURAL = np.array(
    [
        [_RULE_A, _RULE_A, _RULE_A],
        [_RULE_B, _RULE_A, _RULE_A],
        [_RULE_A, _RULE_B, _RULE_A],
        [_RULE_A, _RULE_A, _RULE_B],
    ]
)
"""(xi, eta, zeta) of the four points, along the edges from corner 1 to corners 2, 3, 4."""


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
    orientations: dict[str, str] = field(default_factory=dict[str, str])
    """Element set of a section with a ``MaterialFrame`` → its ``*ORIENTATION`` name."""

    def as_metadata(self) -> dict[str, str]:
        return {
            f"ccx_name.{kind}.{neutral}": ccx
            for kind, table in (
                ("elset", self.element_sets),
                ("nset", self.node_sets),
                ("material", self.materials),
                ("orientation", self.orientations),
            )
            for neutral, ccx in table.items()
        }


def deck_names(job: SolveJob) -> DeckNames:
    """The CalculiX names ``write_deck`` uses for ``job``'s sets, materials and orientations."""
    return DeckNames(
        element_sets=_safe_names(job.mesh.element_sets, "E_"),
        node_sets=_safe_names(job.mesh.node_sets, "N_"),
        materials=_safe_names((m.name for m in job.materials), "M_"),
        orientations=_safe_names(
            (s.element_set for s in job.sections if s.orientation is not None), "O_"
        ),
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
    materials, out-of-range indices, unsupported analysis). An ``OrthotropicMaterial`` is written
    as ``*ELASTIC, TYPE=ENGINEERING CONSTANTS`` (ADR-009 D2) and each
    section with an orientation is written as ``*ORIENTATION`` plus ``ORIENTATION=`` on its
    ``*SOLID SECTION``; a job without orientations gives the same deck as before ADR-009.
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
    orients = names.orientations

    out: list[str] = ["** Tessera Studio: generated CalculiX deck"]
    out.extend(
        f"** {kind} {neutral!r} -> {ccx}"
        for kind, table in (
            ("ELSET", elsets),
            ("NSET", nsets),
            ("MATERIAL", mats),
            ("ORIENTATION", orients),
        )
        for neutral, ccx in table.items()
    )
    out.append("*NODE")
    out.extend(
        f"{i + 1}, {_num(x)}, {_num(y)}, {_num(z)}" for i, (x, y, z) in enumerate(mesh.nodes)
    )
    # Without point stresses the deck is exactly what it was before ADR-008.
    out.append(
        f"*ELEMENT, TYPE=C3D10, ELSET={ALL_ELSET}" if job.point_stresses else "*ELEMENT, TYPE=C3D10"
    )
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
        if isinstance(mat, OrthotropicMaterial):
            # ccx reads 8 constants per line: E1 E2 E3 nu12 nu13 nu23 G12 G13, then G23.
            first = (mat.e1, mat.e2, mat.e3, mat.nu12, mat.nu13, mat.nu23, mat.g12, mat.g13)
            out.append("*ELASTIC, TYPE=ENGINEERING CONSTANTS")
            out.append(", ".join(_num(v) for v in first))
            out.append(_num(mat.g23))
        else:
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
        section = f"*SOLID SECTION, ELSET={elsets[sec.element_set]}, MATERIAL={mats[sec.material]}"
        if sec.orientation is not None:
            name = orients[sec.element_set]
            a, b = sec.orientation.axis_1, sec.orientation.axis_2
            out.append(f"*ORIENTATION, NAME={name}, SYSTEM=RECTANGULAR")
            out.append(", ".join(_num(v) for v in (*a, *b)))
            section += f", ORIENTATION={name}"
        out.append(section)
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
    if job.point_stresses:
        out.extend([f"*EL PRINT, ELSET={ALL_ELSET}", "S"])
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
    finite = np.isfinite(displacements[used]).all(axis=1) & np.isfinite(stresses[used]).all(axis=1)
    bad = used[~finite]
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
    finite = np.isfinite(point_stresses.points).all(axis=(1, 2))
    finite &= np.isfinite(point_stresses.stresses).all(axis=(1, 2))
    bad = np.flatnonzero(~finite)
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
_DAT_MARK = b" (elem, integ.pnt.,"
_DAT_HEADER = re.compile(
    rb" (?P<name>[a-z][a-z ]*?) \(elem, integ\.pnt\.,(?P<columns>[^)\r\n]*)\) for set \S+ and time"
)
_DAT_FIRST_ROW = re.compile(rb"\S")
_DAT_BLANK_LINE = re.compile(rb"\n[ \t\r]*\n")


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


def parse_dat(data: str | bytes | mmap.mmap) -> DatResults:
    """Parse the integration-point blocks of a CalculiX ``.dat`` file.

    Each block has a header line, a blank line and one row per element and point: element id,
    point number, then the values, all separated by spaces. Only the last block of each name
    is parsed (ccx writes a block per increment and step), one block at a time, with numpy's
    C parser and without a Python token list, so ``data`` can be a memory-mapped file. Stresses
    are reordered from ``.dat``'s sxx syy szz sxy sxz syz to the neutral xx yy zz xy yz zx.
    Raises ``SolverError`` for a malformed row, a point number below 1, a duplicate row, or an
    element that does not have exactly points 1..Q, Q being the same for every element.
    """
    buffer = data.encode() if isinstance(data, str) else data
    headers = _dat_headers(buffer)
    last = {name: k for k, (name, _, _, _) in enumerate(headers)}
    result = DatResults()
    for name, k in last.items():
        columns, end = headers[k][1], headers[k][3]
        if name == DAT_STRESSES and columns != _DAT_STRESS_COLUMNS:
            raise SolverError(f"CalculiX job.dat: unexpected stress columns {columns!r}")
        limit = headers[k + 1][2] if k + 1 < len(headers) else len(buffer)
        table = _dat_table(name, buffer, end, limit, 2 + len(columns.split(",")))
        block = _point_block(name, table)
        if name == DAT_STRESSES:
            block = DatBlock(block.elements, block.values[..., _DAT_TO_NEUTRAL])
        result.blocks[name] = block
    return result


def _dat_headers(buffer: bytes | mmap.mmap) -> list[tuple[str, str, int, int]]:
    """(name, columns, line start, line end) of every block header, in file order.

    A plain substring search, not a multi-line regex over the whole file: the file is mostly
    data rows, and a regex anchored at each line start costs about as much as parsing them.
    """
    headers: list[tuple[str, str, int, int]] = []
    pos = 0
    while (mark := buffer.find(_DAT_MARK, pos)) >= 0:
        start = buffer.rfind(b"\n", 0, mark) + 1
        end = buffer.find(b"\n", mark)
        end = len(buffer) if end < 0 else end
        header = _DAT_HEADER.match(buffer[start:end])
        if header:
            headers.append((header["name"].decode(), header["columns"].decode(), start, end))
        pos = end
    return headers


def _dat_table(
    name: str, buffer: bytes | mmap.mmap, start: int, limit: int, width: int
) -> NDArray[np.float64]:
    """(R, width) rows of the block whose header line ends at ``start``, up to a blank line."""
    first = _DAT_FIRST_ROW.search(buffer, start, limit)  # skips the blank line after the header
    if first is None:
        return np.zeros((0, width))
    blank = _DAT_BLANK_LINE.search(buffer, first.start(), limit)
    end = blank.start() if blank else limit
    while end > first.start() and buffer[end - 1 : end].isspace():
        end -= 1
    rows = buffer[first.start() : end]  # the one copy: this block's text
    n_rows = rows.count(b"\n") + 1
    try:
        with warnings.catch_warnings():
            # numpy warns (DeprecationWarning) and stops at text it cannot read as a number.
            warnings.simplefilter("error", DeprecationWarning)
            numbers = np.fromstring(rows, sep=" ")  # type: ignore[call-overload]
    except (ValueError, DeprecationWarning):
        raise SolverError(f"CalculiX job.dat block {name!r} is malformed: not a number") from None
    del rows
    if numbers.size != n_rows * width:
        raise SolverError(
            f"CalculiX job.dat block {name!r} is malformed: expected {width} numbers per row"
        )
    return numbers.reshape(n_rows, width)


def _point_block(name: str, table: NDArray[np.float64]) -> DatBlock:
    """Group the rows (element, point, values...) by element, checking points 1..Q."""
    ids, points = table[:, 0].astype(np.int64), table[:, 1].astype(np.int64)
    if not (np.array_equal(ids, table[:, 0]) and np.array_equal(points, table[:, 1])):
        raise SolverError(f"CalculiX job.dat block {name!r} is malformed: non-integer ids")
    if points.size and points.min() < 1:
        raise SolverError(
            f"CalculiX job.dat block {name!r} has point number {points.min()}; "
            "point numbers start at 1"
        )
    q = int(points.max()) if points.size else 0
    keys = np.sort(ids * (q + 1) + points)
    repeated = keys[1:][keys[1:] == keys[:-1]]
    if repeated.size:
        element, point = divmod(int(repeated[0]), q + 1)
        raise SolverError(
            f"CalculiX job.dat block {name!r} has duplicate rows, e.g. element {element} "
            f"point {point}"
        )
    # No duplicates and every point in 1..Q: an element with Q rows has each point once.
    elements, counts = np.unique(ids, return_counts=True)
    if (counts != q).any():
        raise SolverError(
            f"CalculiX job.dat block {name!r}: every element needs points 1 to {q} once; "
            f"element {elements[counts != q][0]} does not have them"
        )
    order = np.lexsort((points, ids))
    values = table[order, 2:].reshape(len(elements), q, table.shape[1] - 2)
    return DatBlock(elements=elements, values=values)


def rule_points(mesh: Mesh) -> NDArray[np.float64]:
    """(M, 4, 3) physical coordinates (mm) of every element's C3D10 integration points.

    The isoparametric tet10 map of the 4-point rule, in ccx's point order (point k nearest
    corner node k), so it follows curved midside nodes. It equals ccx's ``COORD`` output to
    its print precision (pinned by the ccx tests), without the rounding.
    """
    xi = _RULE_NATURAL
    lam = np.column_stack([1.0 - xi.sum(axis=1), xi])  # (4, 4) barycentric
    corner = [lam[:, k] * (2.0 * lam[:, k] - 1.0) for k in range(4)]
    midside = [4.0 * lam[:, i] * lam[:, j] for i, j in TET10_EDGES]
    shape = np.column_stack(corner + midside)  # (4 points, 10 nodes)
    # One (M, 10) @ (10, 4) product per axis: less memory and faster than an einsum over
    # the gathered (M, 10, 3) node coordinates.
    return np.stack([mesh.nodes[mesh.elements, k] @ shape.T for k in range(3)], axis=-1)


def point_stresses(dat: DatResults, mesh: Mesh) -> PointStresses:
    """Integration-point stresses in mesh element order, at the points ``rule_points`` gives.

    Elements ccx did not print are NaN.
    """
    block = dat.blocks.get(DAT_STRESSES)
    if block is None:
        raise SolverError(f"CalculiX job.dat lacks the {DAT_STRESSES!r} block")
    n_points = _RULE_NATURAL.shape[0]
    if block.values.shape[1] != n_points:
        raise SolverError(
            f"CalculiX job.dat has {block.values.shape[1]} integration points per element; "
            f"C3D10 ({POINT_RULE}) has {n_points}"
        )
    index = block.elements - 1
    if index.size and (index.min() < 0 or index.max() >= mesh.n_elements):
        raise SolverError(
            f"CalculiX job.dat has element ids outside 1..{mesh.n_elements} (the mesh)"
        )
    stresses = np.full((mesh.n_elements, n_points, 6), np.nan)
    stresses[index] = block.values
    return PointStresses(points=rule_points(mesh), stresses=stresses)


def read_point_stresses(path: Path, mesh: Mesh, log: str = "") -> PointStresses:
    """Read ``job.dat`` (memory-mapped, parsed block by block).

    A missing or empty file is a ``SolverError``: ccx exits 0 with an empty ``.dat`` when
    ``*EL PRINT`` names no set (ADR-008, finding 7). Parse errors carry ``log`` too.
    """
    if not path.is_file() or path.stat().st_size == 0:
        raise SolverError(
            f"CalculiX wrote no integration-point results ({path.name} is missing or empty)", log
        )
    try:
        with path.open("rb") as file, mmap.mmap(file.fileno(), 0, access=mmap.ACCESS_READ) as data:
            return point_stresses(parse_dat(data), mesh)
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
        """``timeout``: seconds each ``run`` may take before ccx is stopped; ``None``: no limit.

        ``tessera.core.solve`` always passes one (``tessera.core.SOLVE_TIMEOUT`` by default).
        """
        if timeout is not None and not _positive_seconds(timeout):
            raise ValueError(f"timeout must be a positive number of seconds or None, got {timeout}")
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

        returncode, stdout, stderr = _run_ccx([ccx, "-i", JOB_NAME], workdir, self.timeout)
        log = _log(stdout, stderr)
        frd = workdir / f"{JOB_NAME}.frd"
        if "*ERROR" in stdout or "*ERROR" in stderr:
            raise SolverError(f"CalculiX reported an error: {_first_error(log)}", log)
        if returncode != 0:
            raise SolverError(f"CalculiX exited with code {returncode}", log)
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
        points = None
        if job.point_stresses:
            points = read_point_stresses(workdir / f"{JOB_NAME}.dat", job.mesh, log)
            metadata["point_rule"] = POINT_RULE
        check_results_complete(job.mesh, displacements, stresses, points, log)
        return SolveResult(
            displacements=displacements,
            stresses=stresses,
            solver=self.name,
            log=log,
            metadata=metadata,
            point_stresses=points,
        )


_DRAIN_SECONDS = 10.0
"""How long to wait for the output of a stopped ccx before giving up on it."""

_STOP_SIGNALS = ("SIGTERM", "SIGHUP")
"""Signals that end Tessera by default; while ccx runs they stop ccx first (POSIX, main thread)."""


class _Signalled(BaseException):
    """Raised by the temporary signal handler so the waiting code stops ccx before exiting."""

    def __init__(self, signum: int) -> None:
        super().__init__(signum)
        self.signum = signum


def _positive_seconds(value: object) -> bool:
    """A finite number of seconds above zero; ``True`` is not 1 s."""
    return (
        isinstance(value, int | float)
        and not isinstance(value, bool)
        and math.isfinite(value)
        and value > 0
    )


def _run_ccx(cmd: list[str], workdir: Path, timeout: float | None) -> tuple[int, str, str]:
    """Run ccx; return its exit code, stdout and stderr.

    ccx and the processes it starts are killed before this returns or raises when:

    - the run takes longer than ``timeout``: raises ``SolverTimeoutError``;
    - an exception interrupts the wait (Ctrl+C, a GUI cancel): re-raises it;
    - on POSIX, in the main thread, Tessera receives SIGTERM or SIGHUP: the signal is then
      handled as it would have been (by default the process ends);
    - on Linux, Tessera dies in any way, even SIGKILL: the kernel kills the direct child
      (``PR_SET_PDEATHSIG``), but not the processes that child started.

    Not covered: SIGTERM or SIGHUP while a non-main thread waits on macOS (Linux still has
    ``PR_SET_PDEATHSIG``), Tessera killed with SIGKILL on macOS or Windows, and processes
    that leave ccx's process group themselves (``setsid``).
    """
    try:
        proc = subprocess.Popen(
            cmd,
            cwd=workdir,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            # POSIX: a new session and process group, so the whole tree can be killed at once.
            start_new_session=sys.platform != "win32",
            preexec_fn=_die_with_parent() if sys.platform == "linux" else None,
        )
    except OSError as exc:
        raise SolverError(f"Could not start CalculiX ({cmd[0]}): {exc}") from exc
    try:
        with _signals_stop_ccx():
            return _wait(proc, timeout)
    except _Signalled as sig:
        if proc.returncode is None:  # the signal cut short the kill after a timeout
            _kill_tree(proc)
        # The handlers are restored: deliver the signal again to the previous handler, which
        # by default ends the process. If that handler returns, report the stop.
        signal.raise_signal(sig.signum)
        name = signal.Signals(sig.signum).name
        raise SolverError(f"CalculiX was stopped because Tessera received {name}") from None


def _wait(proc: subprocess.Popen[str], timeout: float | None) -> tuple[int, str, str]:
    try:
        stdout, stderr = proc.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        assert timeout is not None  # communicate only times out with a timeout
        _kill_tree(proc)
        stdout, stderr = _drain(proc)
        raise SolverTimeoutError(
            f"CalculiX did not finish within {timeout:g} s and was stopped",
            timeout,
            _log(stdout, stderr),
        ) from None
    except BaseException:
        _kill_tree(proc)
        _drain(proc)
        raise
    return proc.returncode, stdout, stderr


@contextlib.contextmanager
def _signals_stop_ccx() -> Iterator[None]:
    """While ccx runs, make SIGTERM and SIGHUP raise ``_Signalled`` instead of ending Tessera
    at once, so ``_wait`` kills ccx first. Only on POSIX and in the main thread (Python runs
    signal handlers there only); signals set to be ignored (``nohup``) stay ignored."""
    if sys.platform == "win32" or threading.current_thread() is not threading.main_thread():
        yield
        return
    raised = False

    def handler(signum: int, _frame: object) -> None:
        nonlocal raised
        if not raised:  # a second signal must not interrupt the cleanup of the first
            raised = True
            raise _Signalled(signum)

    previous: dict[signal.Signals, Any] = {}
    for name in _STOP_SIGNALS:
        signum = signal.Signals[name]
        old = signal.getsignal(signum)
        if old is signal.SIG_DFL or callable(old):
            previous[signum] = signal.signal(signum, handler)
    try:
        yield
    finally:
        for signum, old in previous.items():
            signal.signal(signum, old)


def _die_with_parent() -> Callable[[], None] | None:
    """A ``preexec_fn`` that makes Linux kill the child with SIGKILL when Tessera dies."""
    if sys.platform != "linux":
        return None
    try:
        libc = ctypes.CDLL(None, use_errno=True)
        prctl = libc.prctl
    except (OSError, AttributeError):
        return None
    parent = os.getpid()

    def preexec() -> None:
        prctl(_PR_SET_PDEATHSIG, int(signal.SIGKILL), 0, 0, 0)
        # If Tessera died before prctl ran, the child already has a new parent and would never
        # get the signal (the documented race), so stop now.
        if os.getppid() != parent:
            os.kill(os.getpid(), signal.SIGKILL)

    return preexec


_PR_SET_PDEATHSIG = 1  # <linux/prctl.h>


def _taskkill() -> str:
    """``taskkill.exe`` by absolute path when it exists, so PATH cannot redirect it."""
    system = Path(os.environ.get("SYSTEMROOT", r"C:\Windows")) / "System32" / "taskkill.exe"
    return str(system) if system.is_file() else "taskkill"


def _kill_tree(proc: subprocess.Popen[str]) -> None:
    """Kill ``proc`` and the processes it started (its process group / process tree)."""
    if sys.platform == "win32":
        # /T: the whole tree by parent pid; /F: forcefully. Works without a Job Object.
        with contextlib.suppress(OSError, subprocess.TimeoutExpired):
            subprocess.run(
                [_taskkill(), "/F", "/T", "/PID", str(proc.pid)],
                capture_output=True,
                check=False,
                timeout=_DRAIN_SECONDS,
            )
    else:
        with contextlib.suppress(ProcessLookupError):  # the group is already gone
            os.killpg(proc.pid, signal.SIGKILL)
    with contextlib.suppress(OSError):  # in case taskkill failed; the process may be gone
        proc.kill()


def _drain(proc: subprocess.Popen[str]) -> tuple[str, str]:
    """Collect what a killed ccx wrote. If a process that escaped the kill keeps the pipes
    open, give up after ``_DRAIN_SECONDS`` and return what was read so far."""
    try:
        return proc.communicate(timeout=_DRAIN_SECONDS)
    except subprocess.TimeoutExpired as exc:
        if sys.platform == "win32":
            # communicate's reader threads hold the pipes' locks, so closing them would block
            # until the escaped process exits. Leave the daemon threads; they read nothing
            # until end of file, so there is no partial output to return.
            return "", ""
        stdout, stderr = _decode(proc.stdout, exc.output), _decode(proc.stderr, exc.stderr)
        for pipe in (proc.stdout, proc.stderr):
            if pipe is not None:
                pipe.close()
        proc.wait()
        return stdout, stderr


def _decode(pipe: IO[str] | None, data: str | bytes | None) -> str:
    """Partial output from ``TimeoutExpired`` is bytes even in text mode."""
    if not data:
        return ""
    if isinstance(data, str):
        return data
    encoding = getattr(pipe, "encoding", None) or "utf-8"
    return data.decode(encoding, errors="replace")


def _log(stdout: str | None, stderr: str | None) -> str:
    stdout, stderr = stdout or "", stderr or ""
    return stdout + (f"\n--- stderr ---\n{stderr}" if stderr else "")


def _first_error(log: str) -> str:
    lines: Sequence[str] = log.splitlines()
    for k, line in enumerate(lines):
        if "*ERROR" in line:
            return " ".join(s.strip() for s in lines[k : k + 3] if s.strip())
    return "see log"
