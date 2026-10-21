"""Edges between faces (ADR-006): where two faces meet, named by the face pair."""

from collections import defaultdict
from collections.abc import Sequence

import numpy as np
from numpy.typing import NDArray

from tessera.model import EdgeInfo, TriangleSurface

Vec3 = tuple[float, float, float]
_STRAIGHT_TOLERANCE = 1e-6
"""A chain is straight when no vertex is further than this fraction of its length off the line."""


def chains(segments: Sequence[tuple[int, int]]) -> list[tuple[int, int] | None]:
    """Group segments (pairs of vertex ids) into connected curves.

    Returns one entry per curve: its two end vertices, or ``None`` when the curve is closed
    (or branches, which a valid solid's face boundary does not). A segment from a vertex to
    itself is a closed curve on its own, e.g. a full circle with one vertex.
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
    ends: dict[int, list[int]] = {}
    for v in sorted(parent):
        odd = ends.setdefault(root(v), [])  # every curve gets an entry, closed ones too
        if degree[v] % 2:
            odd.append(v)
    return [(e[0], e[1]) if len(e) == 2 else None for e in ends.values()]


def ordered_ends(a: Vec3, b: Vec3) -> tuple[Vec3, Vec3]:
    """The two ends in a fixed order (lexicographic), so STEP and STL list them alike."""
    return (a, b) if a <= b else (b, a)


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

    by_pair: dict[tuple[int, int], list[int]] = defaultdict(list)
    for k, (a, b) in enumerate(zip(lo.tolist(), hi.tolist(), strict=True)):
        by_pair[(a, b)].append(k)
    return [
        _surface_edge(surface.vertices, pair, segments[ks]) for pair, ks in sorted(by_pair.items())
    ]


def _surface_edge(
    vertices: NDArray[np.float64], pair: tuple[int, int], segments: NDArray[np.int64]
) -> EdgeInfo:
    a, b = vertices[segments[:, 0]], vertices[segments[:, 1]]
    lengths = np.linalg.norm(b - a, axis=1)
    length = float(lengths.sum())
    curves = chains([(int(s), int(t)) for s, t in segments])
    start = end = None
    kind = "Curved"
    if len(curves) == 1 and curves[0] is not None:
        p, q = vertices[curves[0][0]], vertices[curves[0][1]]
        start, end = ordered_ends(vec3(p), vec3(q))
        used = vertices[np.unique(segments)]
        axis = (q - p) / max(float(np.linalg.norm(q - p)), 1e-300)
        off = np.linalg.norm(np.cross(used - p, axis), axis=1)
        if off.max() <= _STRAIGHT_TOLERANCE * length:
            kind = "Line"
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
