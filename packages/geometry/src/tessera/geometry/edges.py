"""Edges between faces (ADR-006): where two faces meet, named by the face pair."""

from collections import defaultdict
from collections.abc import Sequence
from typing import NamedTuple

import numpy as np
from numpy.typing import NDArray

from tessera.model import EdgeInfo, TriangleSurface

Vec3 = tuple[float, float, float]
_STRAIGHT_TOLERANCE = 1e-6
"""A chain is straight when no vertex is further than this fraction of its length off the line."""
_END_ORDER_TOLERANCE = 1e-9
"""Coordinates of two ends closer than this fraction of the part size count as equal when
ordering the ends, so rounding noise cannot swap them between STEP and STL."""


class Chain(NamedTuple):
    """One connected curve made of segments."""

    ends: tuple[int, int] | None
    """The two end vertices, or ``None`` for a closed (or branching) curve."""
    vertices: list[int]
    """Every vertex on the curve, ascending."""


def chains(segments: Sequence[tuple[int, int]]) -> list[Chain]:
    """Group segments (pairs of vertex ids) into connected curves.

    A valid solid's face boundary does not branch. A segment from a vertex to itself is a
    closed curve on its own, e.g. a full circle with one vertex.
    """
    parent: dict[int, int] = {}

    def root(v: int) -> int:
        parent.setdefault(v, v)
        while parent[v] != v:
            parent[v] = parent[parent[v]]
            v = parent[v]
        return v

    degree: dict[int, int] = defaultdict(int)
    for a, b in segments:
        parent[root(a)] = root(b)
        degree[a] += 1
        degree[b] += 1
    members: dict[int, list[int]] = {}
    for v in sorted(parent):
        members.setdefault(root(v), []).append(v)
    found = []
    for vertices in members.values():
        odd = [v for v in vertices if degree[v] % 2]
        found.append(Chain((odd[0], odd[1]) if len(odd) == 2 else None, vertices))
    return found


def ordered_ends(a: Vec3, b: Vec3, size: float) -> tuple[Vec3, Vec3]:
    """The two ends in a fixed order, so STEP and STL list them alike.

    Ordered by x, then y, then z; coordinates within a tiny fraction of the part ``size`` of
    each other count as equal.
    """
    tolerance = _END_ORDER_TOLERANCE * size
    for p, q in zip(a, b, strict=True):
        if abs(p - q) > tolerance:
            return (a, b) if p < q else (b, a)
    return a, b


def vec3(x: NDArray[np.float64]) -> Vec3:
    return (float(x[0]), float(x[1]), float(x[2]))


def describe_surface_edges(surface: TriangleSurface) -> list[EdgeInfo]:
    """Edges between the Tessera faces of a closed triangulated surface, ordered by face pair.

    These are the triangle edges whose two triangles belong to different faces: exactly the
    face boundaries of ADR-005, at the same feature angle.
    """
    tris = surface.triangles
    m = tris.shape[0]
    edges = np.sort(np.concatenate([tris[:, [0, 1]], tris[:, [1, 2]], tris[:, [2, 0]]]), axis=1)
    owner = np.tile(np.arange(m), 3)
    order = np.lexsort((edges[:, 1], edges[:, 0]))
    edges, owner = edges[order], owner[order]
    # A closed surface has every edge exactly twice; the copies are neighbours after sorting.
    shared = np.flatnonzero(np.all(edges[:-1] == edges[1:], axis=1))
    fa = surface.face_ids[owner[shared]]
    fb = surface.face_ids[owner[shared + 1]]
    boundary = fa != fb
    lo, hi = np.minimum(fa, fb)[boundary], np.maximum(fa, fb)[boundary]
    segments = edges[shared[boundary]]
    span = surface.vertices.max(axis=0) - surface.vertices.min(axis=0) if m else np.zeros(3)
    size = float(np.linalg.norm(span))

    by_pair: dict[tuple[int, int], list[int]] = defaultdict(list)
    for k, (a, b) in enumerate(zip(lo.tolist(), hi.tolist(), strict=True)):
        by_pair[(a, b)].append(k)
    return [
        _surface_edge(surface.vertices, pair, segments[ks], size)
        for pair, ks in sorted(by_pair.items())
    ]


def _straight(vertices: NDArray[np.float64], chain: Chain) -> bool:
    """Whether an open chain lies on the line between its ends."""
    if chain.ends is None:
        return False
    p, q = vertices[chain.ends[0]], vertices[chain.ends[1]]
    chord = float(np.linalg.norm(q - p))
    off = np.linalg.norm(
        np.cross(vertices[chain.vertices] - p, (q - p) / max(chord, 1e-300)), axis=1
    )
    return bool(off.max() <= _STRAIGHT_TOLERANCE * chord)


def _surface_edge(
    vertices: NDArray[np.float64],
    pair: tuple[int, int],
    segments: NDArray[np.int64],
    size: float,
) -> EdgeInfo:
    a, b = vertices[segments[:, 0]], vertices[segments[:, 1]]
    lengths = np.linalg.norm(b - a, axis=1)
    length = float(lengths.sum())
    curves = chains([(int(s), int(t)) for s, t in segments])
    start = end = None
    if len(curves) == 1 and curves[0].ends is not None:
        p, q = curves[0].ends
        start, end = ordered_ends(vec3(vertices[p]), vec3(vertices[q]), size)
    kind = "Line" if all(_straight(vertices, c) for c in curves) else "Curved"
    # A point on the edge: the length-weighted centroid moved onto the nearest segment, which
    # is the middle of a straight edge.
    centroid = (lengths[:, None] * 0.5 * (a + b)).sum(axis=0) / max(length, 1e-300)
    d = b - a
    t = np.clip(
        np.einsum("ij,ij->i", centroid - a, d) / np.maximum((d * d).sum(axis=1), 1e-300), 0, 1
    )
    closest = a + t[:, None] * d
    point = closest[np.argmin(np.linalg.norm(closest - centroid, axis=1))]
    return EdgeInfo(
        faces=pair,
        kind=kind,
        length=length,
        start=start,
        end=end,
        point=vec3(point),
        curves=len(curves),
    )
