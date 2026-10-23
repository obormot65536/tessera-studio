"""CalculiX (ccx) adapter: writes an ``.inp`` deck, runs ``ccx`` as a subprocess, reads ``.frd``.

ADR-004: the adapter receives a discrete, solver-neutral ``SolveJob`` and returns a
``SolveResult``. Units are mm, N, MPa (ADR-002), so no conversion happens here.

Numbering: Tessera uses 0-based node and element indices; CalculiX uses 1-based ids. Node ``i``
is written as ccx node ``i + 1`` and element ``j`` as ccx element ``j + 1``. The model's tet10
connectivity is already in C3D10 order, so it is written unchanged (plus one).

Precision: the ASCII ``.frd`` file stores results as ``E12.5``, i.e. six significant digits.
Results read back from it carry a relative rounding error of up to 5e-6.

Stopping: ccx runs in its own process group (POSIX) or is stopped with ``taskkill /T``
(Windows), so a timeout or an interrupt also stops the processes it, or a wrapper script set in
``$TESSERA_CCX``, started. See ``_run_ccx`` for what is and is not covered.
"""

import contextlib
import ctypes
import math
import os
import re
import shutil
import signal
import subprocess
import sys
import threading
from collections.abc import Callable, Iterable, Iterator, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import IO, Any

import numpy as np
from numpy.typing import NDArray

from tessera.model.mesh import Mesh
from tessera.model.solver import (
    AnalysisKind,
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
    out.append("*ELEMENT, TYPE=C3D10")
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

    out.extend(["*NODE FILE", "U", "*EL FILE", "S", "*END STEP", ""])
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
    mesh: Mesh, displacements: NDArray[np.float64], stresses: NDArray[np.float64]
) -> None:
    """Raise if a node used by an element has no result (orphan nodes may stay NaN)."""
    used = np.unique(mesh.elements)
    bad = used[np.isnan(displacements[used]).any(axis=1) | np.isnan(stresses[used]).any(axis=1)]
    if bad.size:
        raise SolverError(
            f"CalculiX returned no result for {bad.size} mesh node(s), e.g. node {bad[0]}"
        )


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
        check_results_complete(job.mesh, displacements, stresses)
        return SolveResult(
            displacements=displacements,
            stresses=stresses,
            solver=self.name,
            log=log,
            metadata=metadata,
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
