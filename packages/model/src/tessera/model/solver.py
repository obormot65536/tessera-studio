"""Solver adapter interface (ADR-004).

Physics turns user intent (a force on a face, a fixed face) into a discrete ``SolveJob``;
an adapter in ``tessera.solvers`` turns that job into a solver run and a ``SolveResult``.
Units: mm, N, MPa (ADR-002).
"""

from dataclasses import dataclass, field
from enum import StrEnum
from pathlib import Path
from typing import Protocol

import numpy as np
from numpy.typing import NDArray

from tessera.model.mesh import Mesh


class AnalysisKind(StrEnum):
    STATIC = "static"


@dataclass(frozen=True)
class IsotropicMaterial:
    name: str
    youngs_modulus: float
    """MPa"""
    poisson_ratio: float
    density: float | None = None
    """t/mm^3"""


@dataclass(frozen=True)
class Section:
    element_set: str
    material: str


@dataclass(frozen=True)
class FixedSupport:
    """Zero displacement on ``dofs`` (1 = x, 2 = y, 3 = z) for every node of ``node_set``."""

    node_set: str
    dofs: tuple[int, ...] = (1, 2, 3)
    line: bool = False
    """``True`` when the node set lies on geometry edges (a line support, ADR-006): the stress
    next to it is singular, so its peak value depends on the mesh."""


@dataclass(frozen=True)
class NodalLoads:
    nodes: NDArray[np.int64]
    """(K,) node indices"""
    forces: NDArray[np.float64]
    """(K, 3) forces in N"""


@dataclass(frozen=True)
class SolveJob:
    mesh: Mesh
    materials: tuple[IsotropicMaterial, ...]
    sections: tuple[Section, ...]
    supports: tuple[FixedSupport, ...] = ()
    loads: tuple[NodalLoads, ...] = ()
    analysis: AnalysisKind = AnalysisKind.STATIC


@dataclass(frozen=True)
class PointStresses:
    """Stresses at the integration points of every element, as the solver computed them (ADR-008).

    Row ``j`` belongs to mesh element ``j`` (0-based) and holds that element's own Q points:
    the values are not extrapolated and not averaged with neighbouring elements, so each element
    keeps the stress of its own material. The point rule and its order are the solver's; it is
    named by ``SolveResult.metadata["point_rule"]`` (``"tet4"`` for C3D10: four points, point k
    nearest corner node k). An entry the solver did not report is NaN.

    Integration points lie inside the element, so a maximum over them is biased low at surfaces.
    Reported peaks stay nodal; do not replace them with integration-point maxima (ADR-008,
    decision 4).
    """

    points: NDArray[np.float64]
    """(M, Q, 3) physical coordinates of the points in mm."""
    stresses: NDArray[np.float64]
    """(M, Q, 6) stresses in MPa, neutral order xx, yy, zz, xy, yz, zx."""

    def __post_init__(self) -> None:
        p, s = self.points.shape, self.stresses.shape
        if len(p) != 3 or p[2] != 3:
            raise ValueError(f"points must have shape (M, Q, 3), got {p}")
        if len(s) != 3 or s[2] != 6:
            raise ValueError(f"stresses must have shape (M, Q, 6), got {s}")
        if p[:2] != s[:2]:
            raise ValueError(f"points shape {p} and stresses shape {s} disagree on (M, Q)")

    @property
    def n_elements(self) -> int:
        """M, the number of elements."""
        return int(self.points.shape[0])

    @property
    def n_points(self) -> int:
        """Q, the number of integration points per element."""
        return int(self.points.shape[1])


@dataclass(frozen=True)
class SolveResult:
    displacements: NDArray[np.float64]
    """(N, 3) nodal displacements in mm"""
    stresses: NDArray[np.float64]
    """(N, 6) nodal stresses in MPa: xx, yy, zz, xy, yz, zx"""
    solver: str
    log: str = ""
    metadata: dict[str, str] = field(default_factory=dict)
    """Solver details as text. ``point_rule`` names the integration rule of ``point_stresses``
    (e.g. ``"tet4"``) when they are present."""
    point_stresses: PointStresses | None = None
    """Integration-point stresses (ADR-008); ``None`` when the solver does not provide them.
    Code that needs them (failure criteria) must refuse, not fall back to nodal ``stresses``."""


class SolverError(RuntimeError):
    def __init__(self, message: str, log: str = "") -> None:
        super().__init__(message)
        self.log = log


class SolverAdapter(Protocol):
    name: str

    def capabilities(self) -> frozenset[AnalysisKind]: ...

    def run(self, job: SolveJob, workdir: Path) -> SolveResult: ...
