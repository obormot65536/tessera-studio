"""Solver adapter interface (ADR-004).

Physics turns user intent (a force on a face, a fixed face) into a discrete ``SolveJob``;
an adapter in ``tessera.solvers`` turns that job into a solver run and a ``SolveResult``.
Units: mm, N, MPa (ADR-002).
"""

import math
from dataclasses import dataclass, field
from enum import StrEnum
from pathlib import Path
from typing import Protocol, Self

import numpy as np
from numpy.typing import NDArray

from tessera.model.elasticity import check_orthotropic
from tessera.model.mesh import Mesh

Vector3 = tuple[float, float, float]

_FRAME_TOLERANCE = 1e-9


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
class OrthotropicMaterial:
    """Linear orthotropic material from engineering constants (ADR-009 D1).

    Material axes (see ``MaterialFrame``): 1 = raster/extrusion, 2 = across the road in the
    layer, 3 = build. Moduli in MPa. The Poisson ratios are the major ones,
    nu_ij = -eps_j / eps_i under uniaxial sigma_i (as in Jones, Abaqus and CalculiX). The minor
    ratios are derived, nu_ji = nu_ij E_j / E_i, and never stored. Validated with
    ``tessera.model.elasticity.check_orthotropic``: raises ``ValueError`` when a modulus is not
    positive and finite or the compliance is not positive definite.
    """

    name: str
    e1: float
    e2: float
    e3: float
    nu12: float
    nu13: float
    nu23: float
    g12: float
    g13: float
    g23: float
    density: float | None = None
    """t/mm^3"""

    def __post_init__(self) -> None:
        check_orthotropic(
            self.e1,
            self.e2,
            self.e3,
            self.nu12,
            self.nu13,
            self.nu23,
            self.g12,
            self.g13,
            self.g23,
        )


Material = IsotropicMaterial | OrthotropicMaterial


def _as_vector(v: Vector3 | NDArray[np.float64]) -> Vector3:
    x, y, z = (float(c) for c in v)
    return (x, y, z)


@dataclass(frozen=True)
class MaterialFrame:
    """Material axes 1 and 2 as unit vectors in part coordinates (ADR-009 D2).

    Both must be unit length and orthogonal within 1e-9; axis 3 = axis_1 x axis_2 (right-handed).
    Use ``from_axes`` to build one from any two non-parallel vectors.
    """

    axis_1: Vector3
    axis_2: Vector3

    def __post_init__(self) -> None:
        a, b = _as_vector(self.axis_1), _as_vector(self.axis_2)
        object.__setattr__(self, "axis_1", a)
        object.__setattr__(self, "axis_2", b)
        for name, v in (("axis_1", a), ("axis_2", b)):
            if not all(math.isfinite(c) for c in v):
                raise ValueError(f"MaterialFrame {name} must be finite, got {v!r}")
            if abs(math.hypot(*v) - 1.0) > _FRAME_TOLERANCE:
                raise ValueError(f"MaterialFrame {name} must be a unit vector, got {v!r}")
        if abs(sum(p * q for p, q in zip(a, b, strict=True))) > _FRAME_TOLERANCE:
            raise ValueError(f"MaterialFrame axes must be orthogonal, got {a!r} and {b!r}")

    @property
    def axis_3(self) -> Vector3:
        return _as_vector(np.cross(self.axis_1, self.axis_2))

    @property
    def rotation(self) -> NDArray[np.float64]:
        """3 x 3 matrix whose rows are the material axes 1, 2, 3 in part coordinates."""
        return np.array([self.axis_1, self.axis_2, self.axis_3], dtype=np.float64)

    @classmethod
    def from_axes(cls, a: Vector3, b: Vector3) -> Self:
        """Frame with axis 1 along ``a`` and axis 2 in the plane of ``a`` and ``b`` (Gram-Schmidt).

        Raises ``ValueError`` when ``a`` is zero or ``b`` is (nearly) parallel to it.
        """
        u = np.asarray(a, dtype=np.float64)
        v = np.asarray(b, dtype=np.float64)
        nu, nv = float(np.linalg.norm(u)), float(np.linalg.norm(v))
        if not (np.isfinite(u).all() and np.isfinite(v).all()) or nu == 0 or nv == 0:
            raise ValueError(
                f"MaterialFrame.from_axes needs two non-zero vectors, got {a!r}, {b!r}"
            )
        u = u / nu
        w = v - (v @ u) * u
        nw = float(np.linalg.norm(w))
        if nw <= _FRAME_TOLERANCE * nv:
            raise ValueError(f"MaterialFrame.from_axes: {b!r} is parallel to {a!r}")
        return cls(_as_vector(u), _as_vector(w / nw))


@dataclass(frozen=True)
class Section:
    element_set: str
    material: str
    orientation: MaterialFrame | None = None
    """Material axes for this element set; ``None`` means the part axes (global)."""


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
    materials: tuple[Material, ...]
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


class SolverAdapter(Protocol):
    name: str

    def capabilities(self) -> frozenset[AnalysisKind]: ...

    def run(self, job: SolveJob, workdir: Path) -> SolveResult: ...
