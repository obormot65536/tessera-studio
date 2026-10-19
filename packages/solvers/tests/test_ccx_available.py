import shutil
import subprocess

import pytest


@pytest.mark.solver
def test_calculix_runs() -> None:
    ccx = shutil.which("ccx")
    assert ccx is not None, "CalculiX (ccx) must be on PATH for solver tests"
    out = subprocess.run([ccx, "-v"], capture_output=True, text=True, timeout=30)
    assert "Version 2." in out.stdout
