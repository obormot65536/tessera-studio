"""STL import with Tessera-defined faces (ADR-005)."""

from pathlib import Path

import numpy as np
import trimesh
from numpy.typing import NDArray
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

from tessera.model import FaceInfo, InputError, TriangleSurface

_PLANAR_TOLERANCE = 1e-6
_SLIVER_RATIO = 1e-9
"""A triangle is zero-area when its area is below this fraction of its longest edge squared."""
MIN_FEATURE_ANGLE_DEG = 1.0


def _components(n: int, pairs: NDArray[np.int64]) -> NDArray[np.int64]:
    """Connected-component labels of ``n`` items joined by ``pairs``.

    Components are numbered 0.. in order of their lowest item, so labels are stable for a
    given input order.
    """
    if n == 0:
        return np.zeros(0, dtype=np.int64)
    graph = coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(n, n))
    count, labels = connected_components(graph, directed=False)
    first = np.full(count, n, dtype=np.int64)
    np.minimum.at(first, labels, np.arange(n, dtype=np.int64))
    rank = np.empty(count, dtype=np.int64)
    rank[np.argsort(first, kind="stable")] = np.arange(count)
    return rank[labels]


def read_stl(path: Path, *, feature_angle_deg: float = 40.0) -> TriangleSurface:
    """Read a closed single-body STL (lengths in mm) and split it into faces at sharp edges.

    Triangles with a repeated corner carry no area or topology and are dropped. Triangles with
    three distinct collinear corners would leave a crack if dropped, so the file is rejected.
    """
    if not MIN_FEATURE_ANGLE_DEG <= feature_angle_deg <= 180.0:
        raise InputError(
            f"feature angle must be between {MIN_FEATURE_ANGLE_DEG:g} and 180 degrees, "
            f"got {feature_angle_deg:g}"
        )
    if not path.is_file():
        raise FileNotFoundError(path)
    try:
        loaded = trimesh.load(str(path), file_type="stl", force="mesh", process=True)
    except Exception as exc:  # trimesh raises a variety of parser errors
        raise InputError(f"cannot read {path} as STL: {exc}") from exc
    if not isinstance(loaded, trimesh.Trimesh) or len(loaded.faces) == 0:
        raise InputError(f"{path} contains no triangles")
    mesh = loaded

    faces = np.asarray(mesh.faces, dtype=np.int64)
    needles = (
        (faces[:, 0] == faces[:, 1]) | (faces[:, 1] == faces[:, 2]) | (faces[:, 0] == faces[:, 2])
    )
    if needles.any():
        mesh.update_faces(~needles)
        mesh.remove_unreferenced_vertices()
        if len(mesh.faces) == 0:
            raise InputError(f"{path} contains no triangles with area")

    x = np.asarray(mesh.vertices)[np.asarray(mesh.faces)]
    area2 = np.linalg.norm(np.cross(x[:, 1] - x[:, 0], x[:, 2] - x[:, 0]), axis=1)
    longest = np.max(np.linalg.norm(x - np.roll(x, 1, axis=1), axis=2), axis=1)
    slivers = int((area2 <= _SLIVER_RATIO * longest**2).sum())
    if slivers:
        raise InputError(
            f"{path} contains {slivers} zero-area triangle(s) with three distinct corners; "
            "repair the mesh (e.g. in your CAD or mesh tool) and export again"
        )

    if not mesh.is_watertight:
        raise InputError(f"{path} is not a closed solid (the surface has holes or open edges)")
    if not mesh.is_winding_consistent:
        raise InputError(f"{path} has inconsistent triangle winding")
    adjacency = np.asarray(mesh.face_adjacency, dtype=np.int64)
    if _components(len(mesh.faces), adjacency).max() > 0:
        raise InputError(
            f"{path} has disconnected surfaces (more than one body or shell); "
            "v0.1 supports exactly one solid without internal cavities"
        )

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
    n = surface.n_faces
    index = surface.face_ids - 1
    x = surface.vertices[surface.triangles]
    cross = np.cross(x[:, 1] - x[:, 0], x[:, 2] - x[:, 0])
    areas = 0.5 * np.linalg.norm(cross, axis=1)
    normals = cross / np.maximum(np.linalg.norm(cross, axis=1), 1e-300)[:, None]
    centroids = x.mean(axis=1)

    def per_face(values: NDArray[np.float64]) -> NDArray[np.float64]:
        return np.stack(
            [np.bincount(index, weights=values[:, k], minlength=n) for k in range(3)], axis=1
        )

    face_area = np.bincount(index, weights=areas, minlength=n)
    face_centroid = per_face(areas[:, None] * centroids) / face_area[:, None]
    face_normal = per_face(areas[:, None] * normals)
    length = np.linalg.norm(face_normal, axis=1)
    unit = face_normal / np.maximum(length, 1e-300)[:, None]  # a closed tube sums to zero
    deviation = np.zeros(n)
    np.maximum.at(deviation, index, 1.0 - np.einsum("ij,ij->i", normals, unit[index]))
    planar = (length > 0.0) & (deviation < _PLANAR_TOLERANCE)

    return [
        FaceInfo(
            index=face + 1,
            kind="Plane" if planar[face] else "Curved",
            area=float(face_area[face]),
            centroid=(
                float(face_centroid[face, 0]),
                float(face_centroid[face, 1]),
                float(face_centroid[face, 2]),
            ),
            normal=(
                (float(unit[face, 0]), float(unit[face, 1]), float(unit[face, 2]))
                if planar[face]
                else None
            ),
        )
        for face in range(n)
    ]
