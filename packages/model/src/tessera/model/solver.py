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
class SolveResult:
    displacements: NDArray[np.float64]
    """(N, 3) nodal displacements in mm"""
    stresses: NDArray[np.float64]
    """(N, 6) nodal stresses in MPa: xx, yy, zz, xy, yz, zx"""
    solver: str
    log: str = ""
    metadata: dict[str, str] = field(default_factory=dict)


class SolverError(RuntimeError):
    def __init__(self, message: str, log: str = "") -> None:
        super().__init__(message)
        self.log = log


class SolverTimeoutError(SolverError):
    """The solver did not finish within ``timeout`` seconds and was stopped (ADR-004)."""

    def __init__(self, message: str, timeout: float, log: str = "") -> None:
        super().__init__(message, log)
        self.timeout = timeout


class SolverAdapter(Protocol):
    name: str

    def capabilities(self) -> frozenset[AnalysisKind]: ...

    def run(self, job: SolveJob, workdir: Path) -> SolveResult: ...
