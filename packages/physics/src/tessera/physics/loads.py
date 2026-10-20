"""Convert face and body loads into consistent nodal forces on quadratic elements."""

import numpy as np
from numpy.typing import NDArray

from tessera.model import TET10_EDGES, InputError

# Dunavant degree-4 rule on the triangle: barycentric points and weights (weights sum to 1).
_A, _B = 0.445948490915965, 0.091576213509771
_WA, _WB = 0.223381589678011, 0.109951743655322
_BARY = np.array(
    [
        [1 - 2 * _A, _A, _A],
        [_A, 1 - 2 * _A, _A],
        [_A, _A, 1 - 2 * _A],
        [1 - 2 * _B, _B, _B],
        [_B, 1 - 2 * _B, _B],
        [_B, _B, 1 - 2 * _B],
    ]
)
_W = np.array([_WA, _WA, _WA, _WB, _WB, _WB])


def _tri6_shape(l1: float, l2: float, l3: float) -> tuple[NDArray, NDArray, NDArray]:
    """Shape functions and derivatives wrt (xi, eta) = (l2, l3) for TRI6 node order."""
    n = np.array(
        [
            l1 * (2 * l1 - 1),
            l2 * (2 * l2 - 1),
            l3 * (2 * l3 - 1),
            4 * l1 * l2,
            4 * l2 * l3,
            4 * l3 * l1,
        ]
    )
    dxi = np.array([-(4 * l1 - 1), 4 * l2 - 1, 0.0, 4 * (l1 - l2), 4 * l3, -4 * l3])
    deta = np.array([-(4 * l1 - 1), 0.0, 4 * l3 - 1, -4 * l2, 4 * l2, 4 * (l1 - l3)])
    return n, dxi, deta


_SHAPES = [_tri6_shape(*p) for p in _BARY]


def tri6_load_weights(nodes: NDArray[np.float64], tris: NDArray[np.int64]) -> NDArray[np.float64]:
    """Return (K, 6) integrals of each shape function over each triangle (mm^2).

    Multiplying by a uniform traction gives the consistent nodal forces. The weights of one
    triangle sum to its area.
    """
    x = nodes[tris]  # (K, 6, 3)
    weights = np.zeros(tris.shape, dtype=np.float64)
    for (n, dxi, deta), w in zip(_SHAPES, _W, strict=True):
        jac = np.linalg.norm(np.cross(dxi @ x, deta @ x), axis=1)  # (K,)
        weights += np.outer(0.5 * w * jac, n)
    return weights


def uniform_face_load(
    nodes: NDArray[np.float64], tris: NDArray[np.int64], total_force: NDArray[np.float64]
) -> tuple[NDArray[np.int64], NDArray[np.float64]]:
    """Spread ``total_force`` (N) uniformly per unit area over ``tris``; return nodes and forces."""
    weights = tri6_load_weights(nodes, tris)
    area = weights.sum()
    node_ids, inverse = np.unique(tris.ravel(), return_inverse=True)
    per_node = np.bincount(inverse, weights=weights.ravel(), minlength=node_ids.size)
    forces = np.outer(per_node / area, total_force)
    return node_ids.astype(np.int64), forces


def _sum_per_node(
    connectivity: NDArray[np.int64], weights: NDArray[np.float64]
) -> tuple[NDArray[np.int64], NDArray[np.float64]]:
    """Add (K, n, ...) per-element node weights into (unique nodes, summed weights)."""
    node_ids, inverse = np.unique(connectivity.ravel(), return_inverse=True)
    flat = weights.reshape(connectivity.size, -1)
    summed = np.stack(
        [
            np.bincount(inverse, weights=flat[:, k], minlength=node_ids.size)
            for k in range(flat.shape[1])
        ],
        axis=1,
    )
    return node_ids.astype(np.int64), summed


def tri6_normal_weights(nodes: NDArray[np.float64], tris: NDArray[np.int64]) -> NDArray[np.float64]:
    """Return (K, 6, 3) integrals of each shape function times the area-weighted normal.

    The normal follows the triangle's winding (right hand). Exact for curved faces: the
    integrand has degree 4, which the Dunavant rule integrates exactly.
    """
    x = nodes[tris]  # (K, 6, 3)
    weights = np.zeros((*tris.shape, 3), dtype=np.float64)
    for (n, dxi, deta), w in zip(_SHAPES, _W, strict=True):
        normal = np.cross(dxi @ x, deta @ x)  # (K, 3), length = Jacobian
        weights += 0.5 * w * n[None, :, None] * normal[:, None, :]
    return weights


def pressure_face_load(
    nodes: NDArray[np.float64], tris: NDArray[np.int64], pressure: float
) -> tuple[NDArray[np.int64], NDArray[np.float64]]:
    """Consistent nodal forces (N) of a pressure (MPa) on outward-oriented ``tris``.

    Positive pressure pushes into the body, against the outward normal, at every point of a
    curved face.
    """
    return _sum_per_node(tris, -pressure * tri6_normal_weights(nodes, tris))


# Corners of the four faces of a C3D10 element, and the corner opposite each face.
_TET_FACE_CORNERS = np.array([[0, 1, 2], [0, 1, 3], [1, 2, 3], [0, 2, 3]])
_TET_OPPOSITE = np.array([3, 2, 0, 1])
_TRI6_FLIP = np.array([0, 2, 1, 5, 4, 3])


def orient_outward(
    nodes: NDArray[np.float64], elements: NDArray[np.int64], tris: NDArray[np.int64]
) -> NDArray[np.int64]:
    """Reorder each boundary triangle so its normal points out of the element it bounds.

    ``elements`` should be the elements of one part. A triangle must be a face of exactly one
    of them: a face between two elements has no outside.
    """
    corners = np.unique(tris[:, :3])
    near = elements[np.isin(elements[:, :4], corners).sum(axis=1) >= 3]
    face_corners = near[:, _TET_FACE_CORNERS].reshape(-1, 3)
    opposite = near[:, _TET_OPPOSITE].ravel()
    triples = np.sort(np.vstack([face_corners, tris[:, :3]]), axis=1)
    _, keys = np.unique(triples, axis=0, return_inverse=True)  # one id per corner triple
    keys = keys.ravel()
    face_keys, tri_keys = keys[: face_corners.shape[0]], keys[face_corners.shape[0] :]

    order = np.argsort(face_keys, kind="stable")
    sorted_keys = face_keys[order]
    first = np.searchsorted(sorted_keys, tri_keys, side="left")
    count = np.searchsorted(sorted_keys, tri_keys, side="right") - first
    if (count == 0).any():
        bad = tris[np.flatnonzero(count == 0)[0], :3].tolist()
        raise ValueError(f"triangle with corners {bad} is not a face of any element")
    if (count > 1).any():
        bad = tris[np.flatnonzero(count > 1)[0], :3].tolist()
        raise InputError(
            f"triangle with corners {bad} lies between two elements, so it has no outside; "
            "pressure needs a boundary face"
        )
    owner = order[first]

    x = nodes[tris[:, :3]]
    normal = np.cross(x[:, 1] - x[:, 0], x[:, 2] - x[:, 0])
    inward = np.einsum("ij,ij->i", normal, nodes[opposite[owner]] - x[:, 0]) > 0
    oriented = tris.copy()
    oriented[inward] = tris[inward][:, _TRI6_FLIP]
    return oriented


def _tet10_shape(
    l1: float, l2: float, l3: float
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """C3D10 shape functions and their derivatives wrt (xi, eta, zeta) = (l1, l2, l3)."""
    l0 = 1.0 - l1 - l2 - l3
    lam = np.array([l0, l1, l2, l3])
    dlam = np.array([[-1.0, -1.0, -1.0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])  # (4, 3)
    n = [lam[i] * (2 * lam[i] - 1) for i in range(4)]
    dn = [(4 * lam[i] - 1) * dlam[i] for i in range(4)]
    for i, j in TET10_EDGES:
        n.append(4 * lam[i] * lam[j])
        dn.append(4 * (lam[j] * dlam[i] + lam[i] * dlam[j]))
    return np.array(n), np.array(dn)


# Three-point Gauss-Jacobi rules on [-1, 1] for the weights (1 - x)^2 and (1 - x), from
# scipy.special.roots_jacobi(3, alpha, 0); each is exact to degree 5 (tested).
_JACOBI_2 = (
    np.array([-0.8540119518537006, -0.30599246792329643, 0.41000441977699675]),
    np.array([1.2570908885190917, 1.169970154078929, 0.23960562406864572]),
)
_JACOBI_1 = (
    np.array([-0.8228240809745921, -0.1810662711185305, 0.5753189235216941]),
    np.array([0.8037276549558384, 0.9169644254383448, 0.2793079196058167]),
)


def _tet_rule() -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """27-point conical product rule on the unit tetrahedron, exact to degree 5.

    The collapsed map (u, (1 - u) v, (1 - u)(1 - v) w) has Jacobian (1 - u)^2 (1 - v); Gauss-
    Jacobi points in u and v absorb those factors, Gauss-Legendre points serve w.
    """
    (xu, wu), (xv, wv) = _JACOBI_2, _JACOBI_1
    xw, ww = np.polynomial.legendre.leggauss(3)
    u, v, w = (0.5 * (x + 1.0) for x in np.meshgrid(xu, xv, xw, indexing="ij"))
    weights = np.einsum("i,j,k->ijk", wu / 8, wv / 4, ww / 2).ravel()
    points = np.stack([u, (1 - u) * v, (1 - u) * (1 - v) * w], axis=-1).reshape(-1, 3)
    return points, weights


# Degree 5 is enough: shape function (2) times the curved element's Jacobian (3).
_TET_POINTS, _TET_WEIGHTS = _tet_rule()
_TET_N, _TET_DN = (np.array(a) for a in zip(*(_tet10_shape(*p) for p in _TET_POINTS), strict=True))
# (27, 10) shape functions and the stacked (27 * 3, 10) derivative rows of all points.
_TET_DN_ROWS = _TET_DN.transpose(0, 2, 1).reshape(-1, 10)


def _det3(j: NDArray[np.float64]) -> NDArray[np.float64]:
    return (
        j[..., 0, 0] * (j[..., 1, 1] * j[..., 2, 2] - j[..., 1, 2] * j[..., 2, 1])
        - j[..., 0, 1] * (j[..., 1, 0] * j[..., 2, 2] - j[..., 1, 2] * j[..., 2, 0])
        + j[..., 0, 2] * (j[..., 1, 0] * j[..., 2, 1] - j[..., 1, 1] * j[..., 2, 0])
    )


def tet10_volume_weights(
    nodes: NDArray[np.float64], elements: NDArray[np.int64]
) -> NDArray[np.float64]:
    """Return (M, 10) integrals of each shape function over each element (mm^3)."""
    weights = np.empty(elements.shape, dtype=np.float64)
    for chunk in range(0, elements.shape[0], 100_000):  # bounds memory on large meshes
        x = nodes[elements[chunk : chunk + 100_000]]  # (m, 10, 3)
        jac = (_TET_DN_ROWS @ x).reshape(x.shape[0], -1, 3, 3)  # (m, 27, 3, 3)
        weights[chunk : chunk + 100_000] = (_det3(jac) * _TET_WEIGHTS) @ _TET_N
    return weights


def body_force_load(
    nodes: NDArray[np.float64],
    elements: NDArray[np.int64],
    density: float,
    acceleration: tuple[float, float, float],
) -> tuple[NDArray[np.int64], NDArray[np.float64]]:
    """Consistent nodal forces (N) of a uniform acceleration (mm/s^2) on mass (t/mm^3)."""
    node_ids, volume = _sum_per_node(elements, tet10_volume_weights(nodes, elements))
    return node_ids, density * volume * np.asarray(acceleration, dtype=float)[None, :]
