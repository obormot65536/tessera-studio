"""Results of a solved study, as returned by ``tessera.core.solve``."""

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from tessera.model.mesh import Mesh
from tessera.model.solver import SolveResult


@dataclass(frozen=True)
class StudyResult:
    mesh: Mesh
    solution: SolveResult
    max_displacement: float
    """Largest nodal displacement magnitude, mm."""
    max_von_mises: float
    """Largest nodal von Mises stress, MPa."""
    workdir: Path | None = None
    """Where the solver's input and output files were kept."""

    @property
    def displacement_magnitude(self) -> np.ndarray:
        return np.linalg.norm(self.solution.displacements, axis=1)
