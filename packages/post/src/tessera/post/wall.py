"""Wall thickness at a surface node, in mm and in elements (ADR-007 Amendment 2).

A mesh too coarse to show the stress rise at a singularity can make a flagged peak read settled
in the two-mesh check. How many elements span the wall there tells whether the check could have
seen the rise. The wall is measured by casting rays from the node into the part, along the
inward normal of each distinct face meeting there, against the mesh's outer surface (each TRI6
face split into four flat triangles), and taking the shortest distance to where a ray leaves the
part. A ray that leaves through a symmetry plane (a one-direction support normal to it) has its
distance doubled, because the real part continues on the other side.
"""

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from tessera.model import SolveJob, WallSpan
from tessera.post.flags import _outer_triangles  # pyright: ignore[reportPrivateUsage]
from tessera.post.mesh_check import local_edge_length

FACE_CLUSTER_ANGLE = 30.0
"""Degrees: face normals at a node closer than this are one face, measured along one ray."""
_RAY_START = 1e-3
"""The ray starts this many local edge lengths inside the surface."""
_MIRROR_ALIGNMENT = 0.99
"""|n . e| above which a hit face is normal to a support's direction."""

# The four flat triangles of a TRI6 (corners 0-2, then the midside nodes of edges 01, 12, 20),
# wound like the TRI6.
_SUB_TRIANGLES = np.array([[0, 3, 5], [3, 1, 4], [5, 4, 2], [3, 4, 5]])


@dataclass(frozen=True, eq=False)
class OuterSurface:
    """The outer surface of a solve job's mesh as flat triangles, ready for ray casts."""

    job: SolveJob
    triangles: NDArray[np.int64]
    """(T, 3) node indices, wound so the normals point out of the part."""
    normals: NDArray[np.float64]
    """(T, 3) unit outward normals."""
    mirror_axis: NDArray[np.int64]
    """(T,) 0, 1 or 2: the direction a one-direction support holds on the triangle's face (a
    symmetry plane), or -1."""
    incident_start: NDArray[np.int64]
    """CSR offsets into ``incident``, one row per node."""
    incident: NDArray[np.int64]
    """Triangles at each node."""


def outer_surface(job: SolveJob) -> OuterSurface:
    """Build the outer surface of ``job.mesh`` once, for any number of ``wall_span`` calls."""
    mesh = job.mesh
    tri6 = _outer_triangles(mesh)
    triangles = tri6[:, _SUB_TRIANGLES].reshape(-1, 3)
    x = mesh.nodes
    normals = np.cross(
        x[triangles[:, 1]] - x[triangles[:, 0]], x[triangles[:, 2]] - x[triangles[:, 0]]
    )
    normals /= np.linalg.norm(normals, axis=1, keepdims=True)

    mirror = np.full(len(tri6), -1, dtype=np.int64)
    for support in job.supports:
        if support.line or len(set(support.dofs)) != 1:
            continue
        on = np.zeros(mesh.n_nodes, dtype=bool)
        on[mesh.node_sets[support.node_set]] = True
        mirror[on[tri6[:, :3]].all(axis=1)] = support.dofs[0] - 1
    mirror_axis = np.repeat(mirror, len(_SUB_TRIANGLES))

    by_node = np.argsort(triangles.ravel(), kind="stable")
    counts = np.bincount(triangles.ravel(), minlength=mesh.n_nodes)
    return OuterSurface(
        job=job,
        triangles=triangles,
        normals=normals,
        mirror_axis=mirror_axis,
        incident_start=np.r_[0, np.cumsum(counts)],
        incident=by_node // 3,
    )


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


def _ray_hits(
    surface: OuterSurface, origin: NDArray[np.float64], direction: NDArray[np.float64]
) -> tuple[NDArray[np.float64], NDArray[np.int64]]:
    """Distances along the ray and indices of the triangles it crosses (Moller-Trumbore)."""
    x = surface.job.mesh.nodes
    tri = surface.triangles
    a = x[tri[:, 0]]
    e1, e2 = x[tri[:, 1]] - a, x[tri[:, 2]] - a
    p = np.cross(direction, e2)
    det = np.einsum("ij,ij->i", e1, p)
    scale = np.linalg.norm(e1, axis=1) * np.linalg.norm(e2, axis=1)
    ok = np.abs(det) > 1e-12 * scale
    inv = np.where(ok, 1.0 / np.where(ok, det, 1.0), 0.0)
    s = origin - a
    u = np.einsum("ij,ij->i", s, p) * inv
    q = np.cross(s, e1)
    v = (q @ direction) * inv
    t = np.einsum("ij,ij->i", e2, q) * inv
    tol = 1e-9
    inside = ok & (u >= -tol) & (v >= -tol) & (u + v <= 1.0 + tol)
    hit = np.flatnonzero(inside)
    return t[hit], hit


def wall_span(surface: OuterSurface, node: int) -> WallSpan | None:
    """The wall at ``node``: the shortest way through the part from it, mm and in elements.

    One ray per distinct face at the node (normals within ``FACE_CLUSTER_ANGLE`` count as one),
    into the part along that face's inward normal, from just inside the surface. A ray ends where
    it first leaves the part: the nearest crossing of a triangle facing along the ray. Leaving
    through a symmetry plane doubles the distance. ``elements`` divides by the local median
    corner edge length (``local_edge_length``). ``None`` for a node inside the part, or when no
    ray leaves the part (an open surface).
    """
    start, end = surface.incident_start[node], surface.incident_start[node + 1]
    if start == end:
        return None
    mesh = surface.job.mesh
    point = mesh.nodes[node]
    location = (float(point[0]), float(point[1]), float(point[2]))
    edge = local_edge_length(mesh, location)
    eps = _RAY_START * edge
    faces = _face_directions(surface.normals[surface.incident[start:end]])
    outward = np.sum(faces, axis=0)
    norm = np.linalg.norm(outward)

    best: tuple[float, NDArray[np.float64], bool] | None = None
    for normal in faces:
        direction = -normal
        # Start just inside, away from every face at the node, so no ray grazes them.
        offset = outward / norm if norm > 1e-6 else normal
        origin = point - eps * offset
        t, hit = _ray_hits(surface, origin, direction)
        leaving = (t > eps) & (surface.normals[hit] @ direction > 0.0)
        if not leaving.any():
            continue
        k = int(np.argmin(np.where(leaving, t, np.inf)))
        exit_point = origin + t[k] * direction
        distance = float(np.linalg.norm(exit_point - point))
        axis = int(surface.mirror_axis[hit[k]])
        mirrored = axis >= 0 and abs(float(surface.normals[hit[k], axis])) > _MIRROR_ALIGNMENT
        if mirrored:
            distance *= 2.0
        if best is None or distance < best[0]:
            best = (distance, direction, mirrored)
    if best is None:
        return None
    thickness, direction, mirrored = best
    return WallSpan(
        thickness=thickness,
        elements=thickness / edge,
        direction=(float(direction[0]), float(direction[1]), float(direction[2])),
        mirrored=mirrored,
    )
