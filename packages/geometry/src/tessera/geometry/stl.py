"""STL import with Tessera-defined faces (ADR-005)."""

from pathlib import Path

import numpy as np
import trimesh
from numpy.typing import NDArray

from tessera.model import FaceInfo, InputError, TriangleSurface

_PLANAR_TOLERANCE = 1e-6


def _components(n: int, pairs: NDArray[np.int64]) -> NDArray[np.int64]:
    """Connected-component labels of ``n`` items joined by ``pairs``, numbered 0.. in order of
    each component's lowest item."""
    labels = np.arange(n, dtype=np.int64)
    if pairs.size:
        a, b = pairs[:, 0], pairs[:, 1]
        while True:
            low = np.minimum(labels[a], labels[b])
            before = labels.copy()
            np.minimum.at(labels, a, low)
            np.minimum.at(labels, b, low)
            labels = labels[labels]
            if np.array_equal(labels, before):
                break
    _, numbered = np.unique(labels, return_inverse=True)
    return numbered.astype(np.int64)


def read_stl(path: Path, *, feature_angle_deg: float = 40.0) -> TriangleSurface:
    """Read a closed single-body STL (lengths in mm) and split it into faces at sharp edges."""
    if not path.is_file():
        raise FileNotFoundError(path)
    try:
        loaded = trimesh.load(str(path), file_type="stl", force="mesh", process=True)
    except Exception as exc:  # trimesh raises a variety of parser errors
        raise InputError(f"cannot read {path} as STL: {exc}") from exc
    if not isinstance(loaded, trimesh.Trimesh) or len(loaded.faces) == 0:
        raise InputError(f"{path} contains no triangles")
    mesh = loaded

    adjacency = np.asarray(mesh.face_adjacency, dtype=np.int64)
    if _components(len(mesh.faces), adjacency).max() > 0:
        raise InputError(f"{path} contains more than one body; v0.1 supports exactly one solid")
    if not mesh.is_watertight:
        raise InputError(f"{path} is not a closed solid (the surface has holes or open edges)")
    if not mesh.is_winding_consistent:
        raise InputError(f"{path} has inconsistent triangle winding")

    triangles = np.asarray(mesh.faces, dtype=np.int64)
    if mesh.volume < 0:
        triangles = triangles[:, ::-1].copy()  # inverted file: make normals point out

    smooth = np.asarray(mesh.face_adjacency_angles) <= np.radians(feature_angle_deg)
    face_ids = _components(len(triangles), adjacency[smooth]) + 1
    return TriangleSurface(
        vertices=np.asarray(mesh.vertices, dtype=np.float64),
        triangles=triangles,
        face_ids=face_ids,
    )


def describe_surface(surface: TriangleSurface) -> list[FaceInfo]:
    """Area, centroid and (for planar faces) outward normal of each face."""
    x = surface.vertices[surface.triangles]
    cross = np.cross(x[:, 1] - x[:, 0], x[:, 2] - x[:, 0])
    areas = 0.5 * np.linalg.norm(cross, axis=1)
    normals = cross / np.maximum(np.linalg.norm(cross, axis=1), 1e-300)[:, None]
    centroids = x.mean(axis=1)

    faces = []
    for face in range(1, surface.n_faces + 1):
        mask = surface.face_ids == face
        area = float(areas[mask].sum())
        centroid = (areas[mask, None] * centroids[mask]).sum(axis=0) / area
        mean = (areas[mask, None] * normals[mask]).sum(axis=0)
        length = float(np.linalg.norm(mean))
        if length > 0.0:  # a closed tube of facets has a zero summed normal
            mean /= length
        planar = length > 0.0 and bool((1.0 - normals[mask] @ mean).max() < _PLANAR_TOLERANCE)
        faces.append(
            FaceInfo(
                index=face,
                kind="Plane" if planar else "Curved",
                area=area,
                centroid=(float(centroid[0]), float(centroid[1]), float(centroid[2])),
                normal=(float(mean[0]), float(mean[1]), float(mean[2])) if planar else None,
            )
        )
    return faces
