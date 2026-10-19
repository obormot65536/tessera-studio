"""Descriptions of imported geometry, used for face selection by people and agents."""

from dataclasses import dataclass

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
