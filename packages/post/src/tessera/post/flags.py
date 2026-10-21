"""Flag nodes whose stress keeps changing as the mesh is refined (ADR-007).

In linear elasticity the stress has no finite value at the edges of a support that holds two or
three directions and at sharp inside (re-entrant) edges. Nodes within d = 2 h of those places
are flagged, h being the median corner-to-corner element edge length. Line supports (ADR-006)
are singular too and will be flagged once they exist.
The flags depend only on the mesh and the supports, so every stress result can use them.
"""

import numpy as np
from numpy.typing import NDArray
from scipy.spatial import KDTree

from tessera.model import TET10_EDGES, FixedSupport, FlagKind, Mesh, NodeFlags, SolveJob, StressFlag

# Defaults approved by the maintainer in ADR-007. They decide which number a user trusts:
# change them only with the maintainer's approval (hard rule 6).
FLAG_DISTANCE_FACTOR = 2.0
"""d = FLAG_DISTANCE_FACTOR x h."""
INSIDE_EDGE_ANGLE = 30.0
"""Degrees the outer surface must turn inward at an edge for the edge to be flagged."""

_SUPPORT_PREFIX = "support:"

# Faces of a C3D10 element as TRI6 (corners, then midside nodes of edges ab, bc, ca), and the
# corner opposite each face.
_FACE_TRI6 = np.array(
    [[0, 1, 2, 4, 5, 6], [0, 1, 3, 4, 8, 7], [1, 2, 3, 5, 9, 8], [0, 2, 3, 6, 9, 7]]
)
_FACE_OPPOSITE = np.array([3, 2, 0, 1])
_TRI6_FLIP = np.array([0, 2, 1, 5, 4, 3])
# Edge k of a TRI6 runs from corner k to corner (k + 1) % 3; its midside node is k + 3 and the
# corner opposite it is (k + 2) % 3.
_TRI_EDGE_CORNERS = np.array([[0, 1], [1, 2], [2, 0]])


def _tri6_derivatives() -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """(3, 6) TRI6 shape-function derivatives wrt (xi, eta) = (l2, l3) at each edge midpoint."""
    rows_xi, rows_eta = [], []
    for l1, l2, l3 in [(0.5, 0.5, 0.0), (0.0, 0.5, 0.5), (0.5, 0.0, 0.5)]:
        d1 = np.array([4 * l1 - 1, 0, 0, 4 * l2, 0, 4 * l3])
        d2 = np.array([0, 4 * l2 - 1, 0, 4 * l1, 4 * l3, 0])
        d3 = np.array([0, 0, 4 * l3 - 1, 0, 4 * l2, 4 * l1])
        rows_xi.append(d2 - d1)
        rows_eta.append(d3 - d1)
    return np.array(rows_xi, dtype=np.float64), np.array(rows_eta, dtype=np.float64)


_DXI, _DETA = _tri6_derivatives()


def support_is_singular(support: FixedSupport) -> bool:
    """Whether the stress near this support is flagged: the one place that decides it.

    Today: a support holding two or three directions clamps its face, and the stress at the
    face's edge is singular, so it is flagged. A support holding one direction is an exact
    symmetry plane or roller, so it is not. Line supports (ADR-006) will be flagged regardless
    of direction once the solve job marks them; that check belongs here.
    """
    return len(set(support.dofs)) >= 2


def _pair_keys(a: NDArray[np.int64], b: NDArray[np.int64], n_nodes: int) -> NDArray[np.int64]:
    """One int64 per unordered node pair (min/max, not np.sort: sorting tiny rows is slow)."""
    return np.minimum(a, b).astype(np.int64) * n_nodes + np.maximum(a, b)


def _groups(keys: NDArray[np.int64]) -> tuple[NDArray[np.int64], NDArray[np.int64]]:
    """Sort order of ``keys`` and the start of each run of equal keys in that order."""
    order = np.argsort(keys)
    ordered = keys[order]
    starts = np.flatnonzero(np.r_[True, ordered[1:] != ordered[:-1]])
    return order, starts


def median_edge_length(mesh: Mesh) -> float:
    """h: the median corner-to-corner length of the mesh's element edges, mm."""
    edges = np.array(TET10_EDGES)
    el = mesh.elements
    keys = np.sort(_pair_keys(el[:, edges[:, 0]], el[:, edges[:, 1]], mesh.n_nodes).ravel())
    keys = keys[np.r_[True, keys[1:] != keys[:-1]]]  # np.unique, but several times faster
    a, b = np.divmod(keys, mesh.n_nodes)
    return float(np.median(np.linalg.norm(mesh.nodes[a] - mesh.nodes[b], axis=1)))


def _within(points: NDArray[np.float64], sources: NDArray[np.float64], d: float) -> NDArray:
    """(N,) True where a point is within ``d`` of any source point."""
    if sources.shape[0] == 0:
        return np.zeros(points.shape[0], dtype=bool)
    bound = d * (1.0 + 1e-12)
    distance, _ = KDTree(sources).query(points, distance_upper_bound=bound)
    return np.asarray(distance <= bound)


def _outer_triangles(mesh: Mesh) -> NDArray[np.int64]:
    """(K, 6) element faces on the outer surface, wound so their normals point out."""
    faces = mesh.elements[:, _FACE_TRI6].reshape(-1, 6)
    opposite = mesh.elements[:, _FACE_OPPOSITE].ravel()
    # A face is outer when no other element face has the same three corners. Faces with the
    # same corners share their lowest and highest corner and the sum of all three.
    c = faces[:, :3]
    lo_hi = _pair_keys(c.min(axis=1), c.max(axis=1), mesh.n_nodes)
    total = c.sum(axis=1)
    order = np.lexsort((total, lo_hi))
    lo_hi, total = lo_hi[order], total[order]
    new = np.r_[True, (lo_hi[1:] != lo_hi[:-1]) | (total[1:] != total[:-1])]
    run = np.cumsum(new) - 1
    outer = np.empty(len(faces), dtype=bool)
    outer[order] = np.bincount(run)[run] == 1
    faces, opposite = faces[outer], opposite[outer]
    x = mesh.nodes
    a, b, c = (x[faces[:, k]] for k in range(3))
    inward = np.einsum("ij,ij->i", np.cross(b - a, c - a), x[opposite] - a) > 0
    faces[inward] = faces[inward][:, _TRI6_FLIP]
    return faces


def sharp_inside_edge_nodes(mesh: Mesh, angle: float = INSIDE_EDGE_ANGLE) -> NDArray[np.int64]:
    """Nodes (corners and midside) of outer-surface edges that turn inward by ``angle`` or more.

    The turn is measured between the normals of the two curved TRI6 faces at the midpoint of
    their shared edge, so a smooth curved surface (a hole, a fillet) meshed with curved elements
    does not turn there, however coarse the mesh. An edge is inside (concave) when each face
    rises toward the other's outward normal. An edge shared by more than two outer faces (two
    bodies touching along a line) is singular too and always counts as sharp.
    """
    tris = _outer_triangles(mesh)
    if tris.shape[0] == 0:
        return np.zeros(0, dtype=np.int64)
    x = mesh.nodes[tris]  # (K, 6, 3)
    normals = np.cross(np.einsum("pk,nkd->npd", _DXI, x), np.einsum("pk,nkd->npd", _DETA, x))
    normals = (normals / np.linalg.norm(normals, axis=2, keepdims=True)).reshape(-1, 3)

    # Entry 3t + k is edge k of triangle t; ``pairs`` holds its corners, lowest first.
    a, b = tris[:, _TRI_EDGE_CORNERS[:, 0]].ravel(), tris[:, _TRI_EDGE_CORNERS[:, 1]].ravel()
    pairs = np.stack([np.minimum(a, b), np.maximum(a, b)], axis=1)
    order, starts = _groups(_pair_keys(a, b, mesh.n_nodes))
    counts = np.diff(np.r_[starts, len(order)])
    shared = order[np.repeat(starts[counts > 2], counts[counts > 2])]  # non-manifold edges
    starts = starts[counts == 2]
    e1, e2 = order[starts], order[starts + 1]

    tri1, k1 = np.divmod(e1, 3)
    tri2, k2 = np.divmod(e2, 3)
    mid = mesh.nodes[tris[tri1, k1 + 3]]
    along = mesh.nodes[pairs[e1, 1]] - mesh.nodes[pairs[e1, 0]]
    along /= np.linalg.norm(along, axis=1, keepdims=True)

    def into(tri: NDArray[np.int64], k: NDArray[np.int64]) -> NDArray[np.float64]:
        """Unit direction from the edge into the face, across the edge."""
        v = mesh.nodes[tris[tri, (k + 2) % 3]] - mid
        v -= np.einsum("ij,ij->i", v, along)[:, None] * along
        return v / np.linalg.norm(v, axis=1, keepdims=True)

    n1, n2 = normals[e1], normals[e2]
    concave = (
        np.einsum("ij,ij->i", n1, into(tri2, k2)) + np.einsum("ij,ij->i", n2, into(tri1, k1)) > 0
    )
    turn = np.degrees(np.arccos(np.clip(np.einsum("ij,ij->i", n1, n2), -1.0, 1.0)))
    sharp = concave & (turn >= angle - 1e-9)
    edge = np.r_[e1[sharp], shared]
    tri, k = np.divmod(edge, 3)
    edge_nodes = np.concatenate([pairs[edge].ravel(), tris[tri, k + 3]])
    return np.unique(edge_nodes).astype(np.int64)


def flag_nodes(job: SolveJob) -> NodeFlags:
    """Flag nodes within d = 2 h of a singular support or a sharp inside edge (ADR-007).

    A support is flagged when ``support_is_singular`` says so; its flag names the support's id
    (its node set name without the ``support:`` prefix). Distances are measured to the nodes of
    the support or edge; along an edge they are at most h / 2 apart, which shortens the band
    by a negligible amount.
    """
    mesh = job.mesh
    h = median_edge_length(mesh)
    d = FLAG_DISTANCE_FACTOR * h
    masks: dict[StressFlag, NDArray[np.bool_]] = {}
    for support in job.supports:
        if not support_is_singular(support):
            continue
        flag = StressFlag(FlagKind.SUPPORT, support=support.node_set.removeprefix(_SUPPORT_PREFIX))
        mask = _within(mesh.nodes, mesh.nodes[mesh.node_sets[support.node_set]], d)
        if flag in masks:
            mask |= masks[flag]
        if mask.any():
            masks[flag] = mask
    edge = _within(mesh.nodes, mesh.nodes[sharp_inside_edge_nodes(mesh)], d)
    if edge.any():
        masks[StressFlag(FlagKind.SHARP_INSIDE_EDGE)] = edge
    return NodeFlags(masks=masks, n_nodes=mesh.n_nodes, distance=d, edge_length=h)
