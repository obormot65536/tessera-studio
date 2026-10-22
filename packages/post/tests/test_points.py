"""Stresses from integration points (ADR-008): element-nodal values and point von Mises.

The oracle is an affine stress field: sampled at any four or more points of an element that do
not lie in one plane, the fit must give back the field exactly at the element's ten nodes.
"""

import numpy as np
import pytest
from numpy.typing import NDArray

from tessera.model import TET10_EDGES, Mesh, PointStresses
from tessera.post import element_nodal_stresses, point_von_mises, von_mises

# C3D10 4-point rule: point k has barycentric coordinate B at corner k and A at the other three.
A = (5.0 - np.sqrt(5.0)) / 20.0
B = (5.0 + 3.0 * np.sqrt(5.0)) / 20.0
RULE = np.full((4, 4), A) + np.eye(4) * (B - A)

OFFSET = np.array([1000.0, -500.0, 250.0])
"""The mesh sits far from the origin, as a part on a print bed may."""


def two_curved_tets() -> Mesh:
    """Two distorted tet10 elements sharing the face (0, 1, 2), with curved midside nodes."""
    corners = np.array(
        [[0, 0, 0], [2.0, 0.2, 0], [0.3, 1.8, 0.1], [0.4, 0.5, 1.6], [0.9, 0.6, -1.4]]
    )
    rng = np.random.default_rng(8)
    nodes = [*corners]
    midside: dict[tuple[int, int], int] = {}
    elements = []
    for tet in ([0, 1, 2, 3], [0, 2, 1, 4]):
        conn = list(tet)
        for i, j in TET10_EDGES:
            key = (min(tet[i], tet[j]), max(tet[i], tet[j]))
            if key not in midside:
                midside[key] = len(nodes)
                nodes.append((corners[key[0]] + corners[key[1]]) / 2 + rng.uniform(-0.08, 0.08, 3))
            conn.append(midside[key])
        elements.append(conn)
    return Mesh(nodes=np.array(nodes) + OFFSET, elements=np.array(elements, dtype=np.int64))


def tet10_map(nodes: NDArray[np.float64], barycentric: NDArray[np.float64]) -> NDArray[np.float64]:
    """Isoparametric map of barycentric coordinates (P, 4) through one element's (10, 3) nodes."""
    lam = barycentric
    shape = [lam[:, k] * (2 * lam[:, k] - 1) for k in range(4)]
    shape += [4 * lam[:, i] * lam[:, j] for i, j in TET10_EDGES]
    return np.column_stack(shape) @ nodes


def mapped_points(mesh: Mesh, barycentric: NDArray[np.float64]) -> NDArray[np.float64]:
    return np.stack([tet10_map(mesh.nodes[conn], barycentric) for conn in mesh.elements])


rng = np.random.default_rng(2026)
C = rng.uniform(-50, 50, 6)
G = rng.uniform(-5, 5, (6, 3))


def affine(x: NDArray[np.float64]) -> NDArray[np.float64]:
    """An affine stress field (MPa) of position (mm), all six components varying."""
    return C + (x - OFFSET) @ G.T


def test_affine_field_is_exact_at_the_nodes() -> None:
    mesh = two_curved_tets()
    points = mapped_points(mesh, RULE)
    nodal = element_nodal_stresses(mesh, PointStresses(points, affine(points)))
    assert nodal.shape == (2, 10, 6)
    np.testing.assert_allclose(nodal, affine(mesh.nodes[mesh.elements]), rtol=0, atol=1e-9)


def test_more_than_four_points_are_fitted_by_least_squares() -> None:
    mesh = two_curved_tets()
    centroid = np.full((1, 4), 0.25)
    points = mapped_points(mesh, np.vstack([RULE, centroid]))  # (2, 5, 3)
    nodal = element_nodal_stresses(mesh, PointStresses(points, affine(points)))
    np.testing.assert_allclose(nodal, affine(mesh.nodes[mesh.elements]), rtol=0, atol=1e-9)


def test_elements_are_not_averaged() -> None:
    # Each element keeps its own stress at a node it shares with the other (an interface
    # between two materials keeps both sides).
    mesh = two_curved_tets()
    points = mapped_points(mesh, RULE)
    own = np.array([[1.0, 2, 3, 4, 5, 6], [-6.0, -5, -4, -3, -2, -1]])
    stresses = np.repeat(own[:, None, :], 4, axis=1)
    nodal = element_nodal_stresses(mesh, PointStresses(points, stresses))
    np.testing.assert_allclose(nodal, np.repeat(own[:, None, :], 10, axis=1), atol=1e-12)
    assert mesh.elements[0, 0] == mesh.elements[1, 0]  # node 0 is shared ...
    assert not np.allclose(nodal[0, 0], nodal[1, 0])  # ... with two different stresses


def test_an_element_without_results_stays_nan() -> None:
    mesh = two_curved_tets()
    points = mapped_points(mesh, RULE)
    stresses = affine(points)
    stresses[1, 2, 4] = np.nan
    nodal = element_nodal_stresses(mesh, PointStresses(points, stresses))
    assert np.isnan(nodal[1]).all()
    np.testing.assert_allclose(nodal[0], affine(mesh.nodes[mesh.elements[0]]), atol=1e-9)


def test_points_must_match_the_mesh() -> None:
    mesh = two_curved_tets()
    points = mapped_points(mesh, RULE)
    with pytest.raises(ValueError, match="elements"):
        element_nodal_stresses(mesh, PointStresses(points[:1], affine(points[:1])))
    with pytest.raises(ValueError, match="at least 4"):
        element_nodal_stresses(mesh, PointStresses(points[:, :3], affine(points[:, :3])))


def test_points_in_one_plane_are_refused() -> None:
    mesh = two_curved_tets()
    flat = mapped_points(mesh, RULE)
    flat[0, :, 2] = flat[0, 0, 2]  # element 0's points in one plane
    with pytest.raises(ValueError, match="plane"):
        element_nodal_stresses(mesh, PointStresses(flat, affine(flat)))


def test_point_von_mises() -> None:
    mesh = two_curved_tets()
    points = mapped_points(mesh, RULE)
    ps = PointStresses(points, affine(points))
    vm = point_von_mises(ps)
    assert vm.shape == (2, 4)
    np.testing.assert_allclose(vm, von_mises(ps.stresses.reshape(-1, 6)).reshape(2, 4))
