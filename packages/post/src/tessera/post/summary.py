"""Scalar summaries of a solved study."""

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


def summarize(mesh: Mesh, solution: SolveResult) -> StudyResult:
    return StudyResult(
        mesh=mesh,
        solution=solution,
        max_displacement=float(np.linalg.norm(solution.displacements, axis=1).max()),
        max_von_mises=float(von_mises(solution.stresses).max()),
    )
