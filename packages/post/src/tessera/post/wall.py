"""Wall thickness at a node, in mm and in elements (ADR-007 Amendment 2).

A mesh too coarse to show the stress rise at a singularity can make a flagged peak read settled
in the two-mesh check. How many elements span the wall there tells whether the check could have
seen the rise.

The wall is the diameter of the largest ball inside the part that touches the surface at the
node (a "shrinking ball", as used for medial axes): the ball's centre moves in from the node
along a face's inward normal, and its radius shrinks to r = |p - q|^2 / (2 d . (q - p)), the
ball through the node p and the surface point q nearest the centre, until no surface point is
inside it. Where several faces meet at the node (an edge or a corner), there is one ball per
face and the smallest counts.

The surface is the mesh's outer element faces, each TRI6 split into four flat triangles; their
centroids, with one k-d tree, are the surface samples, and the nearest point is refined to the
sample triangle's plane. A sample stops the ball only where the ball would leave the part
through it (its outward normal faces away from the ball's centre, or along the ball's
direction); so at a convex edge the neighbouring face, which the ball touches edge-on, does not
stop it, and at an inside edge the other face meeting there, coplanar with that face of the
node, is skipped. A support never doubles the wall: a roller cannot be told from a symmetry
plane, and doubling would let a thin wall pass the guard.
"""

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy.spatial import KDTree

from tessera.model import Mesh, WallSpan
from tessera.post.flags import _outer_triangles  # pyright: ignore[reportPrivateUsage]
from tessera.post.mesh_check import node_edge_length

FACE_CLUSTER_ANGLE = 30.0
"""Degrees: face normals at a node closer than this are one face, with one ball."""

# The four flat triangles of a TRI6 (corners 0-2, then the midside nodes of edges 01, 12, 20),
# wound like the TRI6.
_SUB_TRIANGLES = np.array([[0, 3, 5], [3, 1, 4], [5, 4, 2], [3, 4, 5]])
_FACING = 0.5
"""cos 60 degrees: how squarely the ball must meet a sample's surface to be stopped by it."""
_COPLANAR = 0.1
"""A sample this many local edge lengths from the plane of another face at the node lies on it."""
_MAX_STEPS = 200


@dataclass(frozen=True, eq=False)
class OuterSurface:
    """The outer surface of a mesh as flat triangles, with a k-d tree over their centroids."""

    mesh: Mesh
    triangles: NDArray[np.int64]
    """(T, 3) node indices, wound so the normals point out of the part."""
    centroids: NDArray[np.float64]
    """(T, 3) the surface samples."""
    normals: NDArray[np.float64]
    """(T, 3) unit outward normals."""
    sizes: NDArray[np.float64]
    """(T,) longest edge of each triangle."""
    max_size: float
    tree: KDTree
    """Over ``centroids``."""
    incident_start: NDArray[np.int64]
    """CSR offsets into ``incident``, one row per mesh node."""
    incident: NDArray[np.int64]
    """Triangles at each node."""
    surface_nodes: NDArray[np.int64]
    """The nodes on the surface."""
    extent: float
    """Diagonal of the mesh's bounding box, mm: the first ball's radius."""


def outer_surface(mesh: Mesh) -> OuterSurface:
    """Build the outer surface of ``mesh`` once, for any number of ``wall_span`` calls."""
    return _surface(mesh, _outer_triangles(mesh))


def _surface(mesh: Mesh, tri6: NDArray[np.int64]) -> OuterSurface:
    """The surface made of the given outward-wound TRI6 faces."""
    triangles = tri6[:, _SUB_TRIANGLES].reshape(-1, 3)
    corners = mesh.nodes[triangles]  # (T, 3, 3)
    normals = np.cross(corners[:, 1] - corners[:, 0], corners[:, 2] - corners[:, 0])
    normals /= np.linalg.norm(normals, axis=1, keepdims=True)
    edges = np.linalg.norm(corners - np.roll(corners, 1, axis=1), axis=2)
    centroids = corners.mean(axis=1)
    flat = triangles.ravel()
    counts = np.bincount(flat, minlength=mesh.n_nodes)
    span = mesh.nodes.max(axis=0) - mesh.nodes.min(axis=0) if mesh.n_nodes else np.zeros(3)
    return OuterSurface(
        mesh=mesh,
        triangles=triangles,
        centroids=centroids,
        normals=normals,
        sizes=edges.max(axis=1),
        max_size=float(edges.max()) if len(edges) else 0.0,
        tree=KDTree(centroids),
        incident_start=np.r_[0, np.cumsum(counts)],
        incident=np.argsort(flat, kind="stable") // 3,
        surface_nodes=np.flatnonzero(counts),
        extent=float(np.linalg.norm(span)),
    )


def samples_on(surface: OuterSurface, nodes: NDArray[np.int64]) -> NDArray[np.bool_]:
    """(T,) True for the surface samples whose triangle lies on ``nodes`` (e.g. a support's node
    set): all three of its nodes are in it."""
    on = np.zeros(surface.mesh.n_nodes, dtype=bool)
    on[nodes] = True
    return on[surface.triangles].all(axis=1)


def _face_directions(normals: NDArray[np.float64]) -> list[NDArray[np.float64]]:
    """Mean unit normal of each group of ``normals`` within ``FACE_CLUSTER_ANGLE`` of the group."""
    limit = np.cos(np.radians(FACE_CLUSTER_ANGLE))
    sums: list[NDArray[np.float64]] = []
    for n in normals:
        for total in sums:
            if n @ total > limit * np.linalg.norm(total):
                total += n
                break
        else:
            sums.append(n.copy())
    return [total / np.linalg.norm(total) for total in sums]


def _nearest_stop(
    surface: OuterSurface,
    p: NDArray[np.float64],
    d: NDArray[np.float64],
    others: NDArray[np.float64],
    centre: NDArray[np.float64],
    radius: float,
    edge: float,
    ignore: NDArray[np.bool_] | None,
) -> NDArray[np.float64] | None:
    """The surface point nearest ``centre`` that stops the ball, if one is inside it.

    Each sample's point is the foot of ``centre`` on the sample triangle's plane when that foot
    is within the triangle's size of the centroid (exact on flat faces), else the centroid.
    """
    total = len(surface.centroids)
    reach = radius + surface.max_size
    k = min(16, total)
    while True:
        distance, index = surface.tree.query(centre, k=k)
        distance, index = np.atleast_1d(distance), np.atleast_1d(index)
        q, n = surface.centroids[index], surface.normals[index]
        ahead = (q - p) @ d > 1e-6 * edge  # in front of the node's own tangent plane
        away = q - centre
        leaving = (np.einsum("ij,ij->i", away, n) > _FACING * np.linalg.norm(away, axis=1)) | (
            n @ d > _FACING
        )
        stop = ahead & leaving
        if ignore is not None:
            stop &= ~ignore[index]
        if len(others):
            near_other = (n @ others.T).max(axis=1) > np.cos(np.radians(FACE_CLUSTER_ANGLE))
            on_other = np.abs(np.einsum("ij,ij->i", q - p, n)) <= _COPLANAR * edge
            stop &= ~(near_other & on_other)
        foot = centre - np.einsum("ij,ij->i", centre - q, n)[:, None] * n
        use_foot = (np.linalg.norm(foot - q, axis=1) <= surface.sizes[index]) & (
            (foot - p) @ d > 1e-6 * edge
        )
        point = np.where(use_foot[:, None], foot, q)
        gap = np.linalg.norm(point - centre, axis=1)
        stop &= gap < radius * (1.0 - 1e-12)
        if stop.any():
            return point[np.flatnonzero(stop)[np.argmin(gap[stop])]]
        if distance[-1] >= reach or k == total:
            return None
        k = min(4 * k, total)


def _ball_radius(
    surface: OuterSurface,
    p: NDArray[np.float64],
    d: NDArray[np.float64],
    others: NDArray[np.float64],
    edge: float,
    ignore: NDArray[np.bool_] | None,
) -> float | None:
    """Radius of the largest ball touching the surface at ``p`` with its centre along ``d``;
    ``None`` if nothing stops it (an open surface)."""
    radius, stopped = surface.extent, False
    for _ in range(_MAX_STEPS):
        q = _nearest_stop(surface, p, d, others, p + radius * d, radius, edge, ignore)
        if q is None:
            break
        stopped = True
        shrunk = float((q - p) @ (q - p) / (2.0 * ((q - p) @ d)))
        if shrunk >= radius:
            break
        radius = shrunk
    return radius if stopped else None


def wall_span(
    surface: OuterSurface, node: int, ignore: NDArray[np.bool_] | None = None
) -> WallSpan | None:
    """The wall at ``node``: the diameter of the largest ball inside the part touching the
    surface there, in mm and in elements.

    One ball per distinct face at the node (normals within ``FACE_CLUSTER_ANGLE`` count as one),
    centred along that face's inward normal; the smallest counts. A node inside the part is
    measured at the nearest surface node. ``elements`` divides by the local median corner edge
    length at the node measured (``node_edge_length``). ``None`` if no ball is stopped by the
    surface (an open surface), or the mesh has no surface.

    ``ignore`` (from ``samples_on``) leaves samples out: for a spot flagged by a support, the
    faces it holds, where the part continues into the fixture.
    """
    mesh = surface.mesh
    if surface.incident_start[node] == surface.incident_start[node + 1]:
        if surface.surface_nodes.size == 0:
            return None
        gap = np.linalg.norm(mesh.nodes[surface.surface_nodes] - mesh.nodes[node], axis=1)
        node = int(surface.surface_nodes[int(np.argmin(gap))])
    p = mesh.nodes[node]
    edge = node_edge_length(mesh, node)
    faces = _face_directions(
        surface.normals[
            surface.incident[surface.incident_start[node] : surface.incident_start[node + 1]]
        ]
    )
    best: tuple[float, NDArray[np.float64]] | None = None
    for k, normal in enumerate(faces):
        others = np.array([f for j, f in enumerate(faces) if j != k]).reshape(-1, 3)
        radius = _ball_radius(surface, p, -normal, others, edge, ignore)
        if radius is not None and (best is None or radius < best[0]):
            best = (radius, -normal)
    if best is None:
        return None
    radius, direction = best
    return WallSpan(
        thickness=2.0 * radius,
        elements=2.0 * radius / edge,
        direction=(float(direction[0]), float(direction[1]), float(direction[2])),
    )
