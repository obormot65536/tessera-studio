"""Descriptions of imported geometry, used for face selection by people and agents."""

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

Vec3 = tuple[float, float, float]


@dataclass(frozen=True)
class FaceInfo:
    index: int
    """Geometry face index, as used in ``FaceSet.faces``."""
    kind: str
    """Surface type reported by the CAD kernel, e.g. "Plane", "Cylinder"."""
    area: float
    """mm^2"""
    centroid: Vec3
    """mm"""
    normal: Vec3 | None
    """Outward unit normal for planar faces, ``None`` for curved faces."""


@dataclass(frozen=True, eq=False)
class TriangleSurface:
    """A closed triangulated solid boundary with Tessera-defined faces (ADR-005)."""

    vertices: NDArray[np.float64]
    """(N, 3) coordinates in mm."""
    triangles: NDArray[np.int64]
    """(M, 3) vertex indices, counter-clockwise seen from outside (normals point out)."""
    face_ids: NDArray[np.int64]
    """(M,) face index of each triangle, 1..n_faces, as used in ``FaceSet.faces``."""

    def __post_init__(self) -> None:
        if self.vertices.ndim != 2 or self.vertices.shape[1] != 3:
            raise ValueError(f"vertices must have shape (N, 3), got {self.vertices.shape}")
        if self.triangles.ndim != 2 or self.triangles.shape[1] != 3:
            raise ValueError(f"triangles must have shape (M, 3), got {self.triangles.shape}")
        if self.face_ids.shape != (self.triangles.shape[0],):
            raise ValueError("face_ids needs one entry per triangle")
        if self.triangles.size and (
            self.triangles.min() < 0 or self.triangles.max() >= self.vertices.shape[0]
        ):
            raise ValueError("triangles refer to a vertex index outside the vertices")
        used = np.unique(self.face_ids)
        if used.size and not np.array_equal(used, np.arange(1, used.size + 1)):
            raise ValueError("face_ids must number faces 1..n without gaps")

    @property
    def n_faces(self) -> int:
        return int(self.face_ids.max()) if self.face_ids.size else 0
