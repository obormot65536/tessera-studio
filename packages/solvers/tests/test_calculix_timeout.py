"""CalculiX timeouts with a fake ``ccx``: no real solves, runs on Linux and Windows (ADR-004)."""

import contextlib
import os
import signal
import subprocess
import sys
import threading
import time
from collections.abc import Iterator
from pathlib import Path

import numpy as np
import pytest

from tessera.model.mesh import Mesh
from tessera.model.solver import (
    FixedSupport,
    IsotropicMaterial,
    NodalLoads,
    Section,
    SolveJob,
    SolverError,
    SolverTimeoutError,
)
from tessera.solvers import CalculixAdapter
from tessera.solvers.calculix import CCX_ENV_VAR

HANG_LIMIT = 60.0
"""Seconds a test waits for ``run`` to return before calling it hung."""


def tiny_job() -> SolveJob:
    """One unit tet10, fixed on its x = 0 face, pulled at node 1: a valid deck for ccx."""
    nodes = [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [0.5, 0, 0]]
    nodes += [[0.5, 0.5, 0], [0, 0.5, 0], [0, 0, 0.5], [0.5, 0, 0.5], [0, 0.5, 0.5]]
    mesh = Mesh(
        nodes=np.array(nodes, dtype=np.float64),
        elements=np.arange(10, dtype=np.int64).reshape(1, 10),
        element_sets={"part": np.array([0], dtype=np.int64)},
        node_sets={"fixed": np.array([0, 2, 3, 6, 7, 9], dtype=np.int64)},
    )
    return SolveJob(
        mesh=mesh,
        materials=(IsotropicMaterial("resin", 2500.0, 0.35),),
        sections=(Section("part", "resin"),),
        supports=(FixedSupport("fixed"),),
        loads=(NodalLoads(np.array([1]), np.array([[1.0, 0.0, 0.0]])),),
    )


@pytest.fixture
def pids(tmp_path: Path) -> Iterator[Path]:
    """A file the fake processes append their pids to; whatever still runs is killed after."""
    path = tmp_path / "pids"
    yield path
    if path.is_file():
        for pid in path.read_text().split():
            with contextlib.suppress(OSError):
                os.kill(int(pid), signal.SIGTERM if sys.platform == "win32" else signal.SIGKILL)


def _ticker(heartbeat: Path, pids: Path) -> str:
    """Python code that records its pid, then appends to ``heartbeat`` every 50 ms forever."""
    return (
        "import os, time\n"
        f"with open({str(pids)!r}, 'a') as f: f.write(f'{{os.getpid()}}\\n')\n"
        "while True:\n"
        f"    with open({str(heartbeat)!r}, 'a') as f: f.write('.')\n"
        "    time.sleep(0.05)\n"
    )


def _fake_ccx(where: Path, body: str, wrapper: bool = True) -> Path:
    """An executable that runs ``body`` with this Python, by default through a shell wrapper.

    The wrapper does not ``exec``, so the Python process is a child of the process the adapter
    starts, like a user's wrapper script around ccx. Without it (POSIX only) the Python
    process is the one the adapter starts.
    """
    where.mkdir()
    script = where / "fake_ccx.py"
    script.write_text(body, encoding="utf-8")
    if not wrapper:
        script.write_text(f"#!{sys.executable}\n{body}", encoding="utf-8")
        script.chmod(0o755)
        return script
    if sys.platform == "win32":
        exe = where / "fake_ccx.cmd"
        exe.write_text(f'@"{sys.executable}" "{script}" %*\r\n', encoding="utf-8")
    else:
        exe = where / "fake_ccx"
        exe.write_text(f'#!/bin/sh\n"{sys.executable}" "{script}" "$@"\n', encoding="utf-8")
        exe.chmod(0o755)
    return exe


def _hanging_ccx(where: Path, heartbeat: Path, pids: Path, escape: bool = False) -> Path:
    """A ccx that starts a grandchild ticking ``heartbeat``, waits until it ticks, prints a
    line, then hangs. Waiting first means a kill always finds the grandchild running.

    The grandchild inherits the output pipes, so the run cannot finish until it is dead too.
    With ``escape`` it leaves ccx's process group (``setsid``), so killing the group misses it.
    """
    body = (
        "import os, subprocess, sys, time\n"
        f"with open({str(pids)!r}, 'a') as f: f.write(f'{{os.getpid()}}\\n')\n"
        f"subprocess.Popen([sys.executable, '-c', {_ticker(heartbeat, pids)!r}],"
        f" start_new_session={escape!r})\n"
        f"while not os.path.exists({str(heartbeat)!r}): time.sleep(0.05)\n"
        "print('fake ccx started', flush=True)\n"
        "time.sleep(600)\n"
    )
    return _fake_ccx(where, body)


def _run_with_limit(adapter: CalculixAdapter, workdir: Path) -> BaseException | None:
    """Run the adapter in a thread; fail the test if it is still running after HANG_LIMIT."""
    outcome: list[BaseException | None] = []

    def target() -> None:
        try:
            adapter.run(tiny_job(), workdir)
            outcome.append(None)
        except BaseException as exc:
            outcome.append(exc)

    thread = threading.Thread(target=target, daemon=True)
    thread.start()
    thread.join(HANG_LIMIT)
    assert not thread.is_alive(), f"CalculixAdapter.run still running after {HANG_LIMIT} s"
    return outcome[0]


def _assert_stopped(heartbeat: Path) -> None:
    assert heartbeat.is_file(), "the fake ccx's grandchild never started"
    time.sleep(0.3)  # let a dying process finish its last write
    size = heartbeat.stat().st_size
    time.sleep(0.5)
    assert heartbeat.stat().st_size == size, "a process started by ccx is still running"


def test_timeout_stops_ccx_and_the_processes_it_started(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, pids: Path
) -> None:
    heartbeat = tmp_path / "heartbeat"
    monkeypatch.setenv(CCX_ENV_VAR, str(_hanging_ccx(tmp_path / "bin", heartbeat, pids)))
    started = time.monotonic()
    error = _run_with_limit(CalculixAdapter(timeout=8), tmp_path / "work")
    assert time.monotonic() - started < 40
    assert isinstance(error, SolverTimeoutError)
    assert isinstance(error, SolverError)
    assert str(error) == "CalculiX did not finish within 8 s and was stopped"
    assert error.timeout == 8
    assert "fake ccx started" in error.log
    _assert_stopped(heartbeat)


def test_an_interrupted_run_stops_ccx_too(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, pids: Path
) -> None:
    # Ctrl+C in the CLI, or a GUI cancel that raises in the waiting thread.
    heartbeat = tmp_path / "heartbeat"
    monkeypatch.setenv(CCX_ENV_VAR, str(_hanging_ccx(tmp_path / "bin", heartbeat, pids)))
    real = subprocess.Popen.communicate
    calls: list[float | None] = []

    def interrupted(self: subprocess.Popen[str], input: str | None = None, timeout=None):
        calls.append(timeout)
        if len(calls) == 1:
            while not heartbeat.is_file():
                time.sleep(0.05)
            raise KeyboardInterrupt
        return real(self, input, timeout)

    monkeypatch.setattr(subprocess.Popen, "communicate", interrupted)
    error = _run_with_limit(CalculixAdapter(timeout=600), tmp_path / "work")
    assert isinstance(error, KeyboardInterrupt)
    _assert_stopped(heartbeat)


def test_ccx_errors_still_reach_the_user(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    body = "print(' *ERROR in calinput: unknown keyword')\nraise SystemExit(201)\n"
    monkeypatch.setenv(CCX_ENV_VAR, str(_fake_ccx(tmp_path / "bin", body)))
    with pytest.raises(SolverError, match="unknown keyword") as caught:
        CalculixAdapter(timeout=60).run(tiny_job(), tmp_path / "work")
    assert not isinstance(caught.value, SolverTimeoutError)


@pytest.mark.parametrize("timeout", [0.0, -1.0, float("nan"), float("inf"), True, "60"])
def test_timeout_must_be_a_positive_number_of_seconds(timeout: object) -> None:
    with pytest.raises(ValueError, match="timeout"):
        CalculixAdapter(timeout=timeout)  # type: ignore[arg-type]


def test_no_timeout_is_still_allowed_for_direct_use() -> None:
    assert CalculixAdapter().timeout is None


posix_only = pytest.mark.skipif(sys.platform == "win32", reason="POSIX signals and setsid")


@posix_only
def test_a_process_that_escaped_the_kill_does_not_hang_the_run_or_lose_the_log(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, pids: Path
) -> None:
    from tessera.solvers import calculix

    monkeypatch.setattr(calculix, "_DRAIN_SECONDS", 2.0)
    heartbeat = tmp_path / "heartbeat"
    fake = _hanging_ccx(tmp_path / "bin", heartbeat, pids, escape=True)
    monkeypatch.setenv(CCX_ENV_VAR, str(fake))
    started = time.monotonic()
    error = _run_with_limit(CalculixAdapter(timeout=3), tmp_path / "work")
    assert time.monotonic() - started < 15
    assert isinstance(error, SolverTimeoutError)
    assert "fake ccx started" in error.log


def _runner(tmp_path: Path, fake: Path) -> subprocess.Popen[bytes]:
    """A separate Python process (its own session, so signals stay away from pytest) that
    waits on ``fake`` through ``_run_ccx`` in its main thread."""
    script = tmp_path / "runner.py"
    script.write_text(
        "import sys\n"
        "from pathlib import Path\n"
        "from tessera.solvers.calculix import _run_ccx\n"
        "_run_ccx([sys.argv[1], '-i', 'job'], Path(sys.argv[2]), 600)\n",
        encoding="utf-8",
    )
    work = tmp_path / "work"
    work.mkdir()
    return subprocess.Popen(
        [sys.executable, str(script), str(fake), str(work)], start_new_session=True
    )


def _wait_for(path: Path) -> None:
    deadline = time.monotonic() + HANG_LIMIT
    while not path.is_file():
        assert time.monotonic() < deadline, f"{path.name} never appeared"
        time.sleep(0.05)


if sys.platform == "win32":
    _ENDING_SIGNALS: list[signal.Signals] = []  # no SIGHUP on Windows; POSIX-only test
else:
    _ENDING_SIGNALS = [signal.SIGTERM, signal.SIGHUP]


@posix_only
@pytest.mark.parametrize("signum", _ENDING_SIGNALS)
def test_tessera_ended_by_a_signal_stops_ccx_first(
    tmp_path: Path, pids: Path, signum: signal.Signals
) -> None:
    # The terminal closing (SIGHUP) or a supervisor stopping Tessera (SIGTERM).
    heartbeat = tmp_path / "heartbeat"
    runner = _runner(tmp_path, _hanging_ccx(tmp_path / "bin", heartbeat, pids))
    try:
        _wait_for(heartbeat)
        time.sleep(0.2)  # the fake prints after the heartbeat appears; let _run_ccx wait
        runner.send_signal(signum)
        assert runner.wait(HANG_LIMIT) == -signum  # ends the way the signal would have
    finally:
        runner.kill()
    _assert_stopped(heartbeat)


@pytest.mark.skipif(sys.platform != "linux", reason="PR_SET_PDEATHSIG is Linux-only")
def test_ccx_dies_with_tessera_even_on_sigkill(tmp_path: Path, pids: Path) -> None:
    # No handler runs on SIGKILL; the kernel stops the direct child (here: the ticking ccx).
    heartbeat = tmp_path / "heartbeat"
    fake = _fake_ccx(tmp_path / "bin", _ticker(heartbeat, pids), wrapper=False)
    runner = _runner(tmp_path, fake)
    try:
        _wait_for(heartbeat)
        runner.kill()
        runner.wait(HANG_LIMIT)
    finally:
        runner.kill()
    _assert_stopped(heartbeat)
