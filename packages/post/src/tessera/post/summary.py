"""Scalar summaries of a solved study."""

from pathlib import Path

import numpy as np
from numpy.typing import NDArray

from tessera.model import Mesh, SolveResult, StudyResult


def von_mises(stresses: NDArray[np.float64]) -> NDArray[np.float64]:
    """Von Mises stress from (N, 6) components xx, yy, zz, xy, yz, zx (MPa)."""
    sxx, syy, szz, sxy, syz, szx = stresses.T
    return np.sqrt(
        0.5 * ((sxx - syy) ** 2 + (syy - szz) ** 2 + (szz - sxx) ** 2)
        + 3.0 * (sxy**2 + syz**2 + szx**2)
    )


def summarize(mesh: Mesh, solution: SolveResult, workdir: Path | None = None) -> StudyResult:
    """Summarize a solution. Nodes no element uses may have no result (NaN) and are skipped."""
    return StudyResult(
        mesh=mesh,
        solution=solution,
        max_displacement=float(np.nanmax(np.linalg.norm(solution.displacements, axis=1))),
        max_von_mises=float(np.nanmax(von_mises(solution.stresses))),
        workdir=workdir,
    )
