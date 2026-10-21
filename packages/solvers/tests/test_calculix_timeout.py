"""CalculiX timeouts with a fake ``ccx``: no real solves, runs on Linux and Windows (ADR-004)."""

import subprocess
import sys
import threading
import time
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


def _fake_ccx(where: Path, body: str) -> Path:
    """An executable that runs ``body`` with this Python through a shell wrapper.

    The wrapper does not ``exec``, so the Python process is a child of the process the adapter
    starts, like a user's wrapper script around ccx.
    """
    where.mkdir()
    script = where / "fake_ccx.py"
    script.write_text(body, encoding="utf-8")
    if sys.platform == "win32":
        exe = where / "fake_ccx.cmd"
        exe.write_text(f'@"{sys.executable}" "{script}" %*\r\n', encoding="utf-8")
    else:
        exe = where / "fake_ccx"
        exe.write_text(f'#!/bin/sh\n"{sys.executable}" "{script}" "$@"\n', encoding="utf-8")
        exe.chmod(0o755)
    return exe


def _hanging_ccx(where: Path, heartbeat: Path) -> Path:
    """A ccx that prints a line, starts a grandchild that ticks ``heartbeat``, then hangs.

    The grandchild inherits the output pipes, so the run cannot finish until it is dead too.
    """
    tick = (
        "import time\n"
        "while True:\n"
        f"    with open({str(heartbeat)!r}, 'a') as f: f.write('.')\n"
        "    time.sleep(0.05)\n"
    )
    body = (
        "import subprocess, sys, time\n"
        "print('fake ccx started', flush=True)\n"
        f"subprocess.Popen([sys.executable, '-c', {tick!r}])\n"
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
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    heartbeat = tmp_path / "heartbeat"
    monkeypatch.setenv(CCX_ENV_VAR, str(_hanging_ccx(tmp_path / "bin", heartbeat)))
    started = time.monotonic()
    error = _run_with_limit(CalculixAdapter(timeout=3), tmp_path / "work")
    assert time.monotonic() - started < 30
    assert isinstance(error, SolverTimeoutError)
    assert isinstance(error, SolverError)
    assert str(error) == "CalculiX did not finish within 3 s and was stopped"
    assert error.timeout == 3
    assert "fake ccx started" in error.log
    _assert_stopped(heartbeat)


def test_an_interrupted_run_stops_ccx_too(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    # Ctrl+C in the CLI, or a GUI cancel that raises in the waiting thread.
    heartbeat = tmp_path / "heartbeat"
    monkeypatch.setenv(CCX_ENV_VAR, str(_hanging_ccx(tmp_path / "bin", heartbeat)))
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


@pytest.mark.parametrize("timeout", [0.0, -1.0, float("nan"), float("inf")])
def test_timeout_must_be_a_positive_number_of_seconds(timeout: float) -> None:
    with pytest.raises(ValueError, match="timeout"):
        CalculixAdapter(timeout=timeout)


def test_no_timeout_is_still_allowed_for_direct_use() -> None:
    assert CalculixAdapter().timeout is None
