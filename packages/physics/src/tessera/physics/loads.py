"""Convert face loads into consistent nodal forces on quadratic triangles."""

import numpy as np
from numpy.typing import NDArray

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
