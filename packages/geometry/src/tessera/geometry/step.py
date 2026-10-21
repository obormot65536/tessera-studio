"""STEP import and face and edge queries (via gmsh's OpenCASCADE kernel)."""

from collections import defaultdict
from collections.abc import Iterator
from pathlib import Path

import gmsh
import numpy as np

from tessera.geometry._gmsh import gmsh_model
from tessera.geometry.edges import chains, ordered_ends, vec3
from tessera.model import EdgeInfo, FaceInfo, InputError

_ZERO_LENGTH = 1e-9
"""Curves shorter than this fraction of the part's size are degenerate and not listed."""


def _import_single_solid(path: Path) -> int:
    if not path.is_file():
        raise FileNotFoundError(path)
    try:
        gmsh.model.occ.importShapes(str(path))
    except Exception as exc:  # gmsh raises plain Exception
        raise InputError(f"cannot read {path} as STEP: {exc}") from exc
    gmsh.model.occ.synchronize()
    volumes = gmsh.model.getEntities(3)
    if len(volumes) != 1:
        raise InputError(f"{path} contains {len(volumes)} solids; v0.1 supports exactly one")
    return volumes[0][1]


def _face_points(face: int, n: int = 7) -> Iterator[np.ndarray]:
    """Points on the (trimmed) face, from a parameter-space grid."""
    (umin, vmin), (umax, vmax) = np.asarray(gmsh.model.getParametrizationBounds(2, face))
    for u in np.linspace(umin, umax, n + 2)[1:-1]:
        for v in np.linspace(vmin, vmax, n + 2)[1:-1]:
            xyz = np.asarray(gmsh.model.getValue(2, face, [u, v]), dtype=float)
            if gmsh.model.isInside(2, face, xyz.tolist()):
                yield xyz


def _outward_normal(face: int, volume: int, centroid: np.ndarray, size: float) -> np.ndarray:
    """Unit normal pointing out of the solid.

    Probes just in front of and behind points on the face; the centroid itself can lie off
    the face (e.g. inside a hole). Falls back to gmsh's face orientation, which is outward
    for valid solids, if no probe is decisive.
    """
    uv = gmsh.model.getParametrization(2, face, centroid.tolist())
    normal = np.asarray(gmsh.model.getNormal(face, uv), dtype=float)
    normal /= np.linalg.norm(normal)
    step = 1e-4 * size
    for point in _face_points(face):
        ahead = gmsh.model.isInside(3, volume, (point + step * normal).tolist())
        behind = gmsh.model.isInside(3, volume, (point - step * normal).tolist())
        if ahead != behind:
            return -normal if ahead else normal
    return normal


def describe_step(path: Path) -> tuple[list[FaceInfo], list[EdgeInfo]]:
    """Faces (by face index) and edges (by face pair) of a single-solid STEP file, read once."""
    with gmsh_model():
        volume = _import_single_solid(path)
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(3, volume)
        size = float(np.linalg.norm([xmax - xmin, ymax - ymin, zmax - zmin]))
        faces = []
        for _, tag in sorted(gmsh.model.getEntities(2)):
            kind = gmsh.model.getType(2, tag)
            centroid = np.asarray(gmsh.model.occ.getCenterOfMass(2, tag))
            normal = None
            if kind == "Plane":
                n = _outward_normal(tag, volume, centroid, size)
                normal = (float(n[0]), float(n[1]), float(n[2]))
            faces.append(
                FaceInfo(
                    index=tag,
                    kind=kind,
                    area=float(gmsh.model.occ.getMass(2, tag)),
                    centroid=(float(centroid[0]), float(centroid[1]), float(centroid[2])),
                    normal=normal,
                )
            )
        return faces, _step_edges(size)


def _step_edges(size: float) -> list[EdgeInfo]:
    """Curves between two different faces, grouped by face pair (ADR-006).

    Seams (the same face on both sides) and zero-length curves are left out.
    """
    by_pair: dict[tuple[int, int], list[tuple[int, float]]] = defaultdict(list)
    for _, curve in sorted(gmsh.model.getEntities(1)):
        up, _ = gmsh.model.getAdjacencies(1, curve)
        sides = sorted({int(f) for f in up})
        if len(sides) != 2:
            continue
        length = float(gmsh.model.occ.getMass(1, curve))
        if length > _ZERO_LENGTH * size:
            by_pair[(sides[0], sides[1])].append((curve, length))
    return [_step_edge(pair, curves, size) for pair, curves in sorted(by_pair.items())]


def _point(tag: int) -> tuple[float, float, float]:
    return vec3(np.asarray(gmsh.model.getValue(0, tag, []), dtype=float))


def _step_edge(pair: tuple[int, int], curves: list[tuple[int, float]], size: float) -> EdgeInfo:
    segments = []
    for curve, _ in curves:
        _, points = gmsh.model.getAdjacencies(1, curve)
        ends = [int(p) for p in points] or [-curve]  # a closed curve may have no vertex
        segments.append((ends[0], ends[-1]))
    connected = chains(segments)
    start = end = None
    if len(connected) == 1 and connected[0].ends is not None:
        p, q = connected[0].ends
        start, end = ordered_ends(_point(p), _point(q), size)
    kinds = {gmsh.model.getType(1, curve) for curve, _ in curves}
    longest = max(curves, key=lambda c: c[1])[0]
    lo, hi = gmsh.model.getParametrizationBounds(1, longest)
    middle = np.asarray(gmsh.model.getValue(1, longest, [0.5 * (lo[0] + hi[0])]), dtype=float)
    return EdgeInfo(
        faces=pair,
        kind=kinds.pop() if len(kinds) == 1 else "Mixed",
        length=sum(length for _, length in curves),
        start=start,
        end=end,
        point=vec3(middle),
        curves=len(connected),
    )
