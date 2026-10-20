"""Pressure on (curved) quadratic faces and gravity on quadratic tetrahedra."""

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from tessera.model import TET10_EDGES
from tessera.physics import (
    body_force_load,
    orient_outward,
    pressure_face_load,
    tet10_volume_weights,
    tri6_normal_weights,
)

# The four faces of a C3D10 tetrahedron as TRI6 (corners, then midsides), outward normals.
_TET10_FACES = np.array(
    [
        [0, 2, 1, 6, 5, 4],
        [0, 1, 3, 4, 8, 7],
        [1, 2, 3, 5, 9, 8],
        [2, 0, 3, 6, 7, 9],
    ]
)
_FLIP = np.array([0, 2, 1, 5, 4, 3])


def _tet10(corners: np.ndarray, bulge: np.ndarray | None = None) -> np.ndarray:
    mids = np.array([0.5 * (corners[i] + corners[j]) for i, j in TET10_EDGES])
    if bulge is not None:
        mids = mids + bulge
    return np.vstack([corners, mids])


_UNIT = np.array([[0.0, 0, 0], [2, 0, 0], [0, 3, 0], [0, 0, 4]])  # volume 4

small = st.floats(min_value=-0.15, max_value=0.15, allow_nan=False)
bulges = st.lists(st.tuples(small, small, small), min_size=6, max_size=6)


def _divergence_volume(nodes: np.ndarray, faces: np.ndarray) -> float:
    """V = 1/3 of the flux of x through the closed boundary (exact for these faces)."""
    from tessera.physics.loads import _SHAPES, _W

    x = nodes[faces]
    total = 0.0
    for (n, dxi, deta), w in zip(_SHAPES, _W, strict=True):
        point = n @ x
        normal = np.cross(dxi @ x, deta @ x)
        total += 0.5 * w * np.einsum("ij,ij->", point, normal)
    return total / 3.0


def test_flat_face_pressure_pushes_against_the_outward_normal() -> None:
    nodes = _tet10(_UNIT)
    bottom = _TET10_FACES[:1]  # z = 0, outward normal -z, area 3
    ids, forces = pressure_face_load(nodes, bottom, 2.0)
    np.testing.assert_allclose(forces.sum(axis=0), [0.0, 0.0, 6.0], atol=1e-12)
    assert set(ids.tolist()) == set(bottom[0].tolist())


def test_normal_weights_of_a_flat_face_sum_to_area_times_normal() -> None:
    nodes = _tet10(_UNIT)
    w = tri6_normal_weights(nodes, _TET10_FACES[:1])[0]
    np.testing.assert_allclose(w.sum(axis=0), [0.0, 0.0, -3.0], atol=1e-12)
    np.testing.assert_allclose(w[:3], 0.0, atol=1e-12)


@settings(max_examples=50)
@given(bulges)
def test_pressure_on_a_closed_curved_surface_has_no_resultant(bulge: list) -> None:
    nodes = _tet10(_UNIT, np.array(bulge))
    _, forces = pressure_face_load(nodes, _TET10_FACES, 5.0)
    np.testing.assert_allclose(forces.sum(axis=0), 0.0, atol=1e-10)


def test_orient_outward_flips_inward_faces_only() -> None:
    nodes = _tet10(_UNIT)
    elements = np.arange(10)[None, :]
    mixed = _TET10_FACES.copy()
    mixed[1] = mixed[1][_FLIP]
    mixed[3] = mixed[3][_FLIP]
    np.testing.assert_array_equal(orient_outward(nodes, elements, mixed), _TET10_FACES)


def test_orient_outward_rejects_faces_that_are_not_on_an_element() -> None:
    nodes = _tet10(_UNIT)
    with pytest.raises(ValueError, match="not a face"):
        orient_outward(nodes, np.arange(10)[None, :], np.array([[0, 1, 9, 4, 8, 7]]))


def test_straight_tet10_volume_weights_are_minus_a_twentieth_and_a_fifth() -> None:
    w = tet10_volume_weights(_tet10(_UNIT), np.arange(10)[None, :])[0]
    np.testing.assert_allclose(w[:4], -4.0 / 20, rtol=1e-12)
    np.testing.assert_allclose(w[4:], 4.0 / 5, rtol=1e-12)


@settings(max_examples=50)
@given(bulges)
def test_curved_tet10_volume_weights_sum_to_its_volume(bulge: list) -> None:
    nodes = _tet10(_UNIT, np.array(bulge))
    w = tet10_volume_weights(nodes, np.arange(10)[None, :])[0]
    assert w.sum() == pytest.approx(_divergence_volume(nodes, _TET10_FACES), rel=1e-10)


def test_body_force_totals_density_times_volume_times_acceleration() -> None:
    nodes = _tet10(_UNIT)
    # Density of steel in t/mm^3, so the total is in N.
    ids, forces = body_force_load(nodes, np.arange(10)[None, :], 7.85e-9, (0.0, 0.0, -9810.0))
    np.testing.assert_allclose(forces.sum(axis=0), [0.0, 0.0, -7.85e-9 * 4.0 * 9810.0])
    assert ids.tolist() == list(range(10))


@pytest.mark.parametrize(
    "a,b,c", [(a, b, c) for a in range(6) for b in range(6 - a) for c in range(6 - a - b)]
)
def test_tet_rule_integrates_every_monomial_up_to_degree_5(a: int, b: int, c: int) -> None:
    from math import factorial

    from tessera.physics.loads import _tet_rule

    points, weights = _tet_rule()
    exact = factorial(a) * factorial(b) * factorial(c) / factorial(a + b + c + 3)
    got = weights @ (points[:, 0] ** a * points[:, 1] ** b * points[:, 2] ** c)
    assert got == pytest.approx(exact, rel=1e-12)


@settings(max_examples=50)
@given(bulges)
def test_pressure_on_a_closed_curved_surface_has_no_net_moment(bulge: list) -> None:
    # Unlike the resultant, the moment depends on how the load is shared between nodes.
    nodes = _tet10(_UNIT, np.array(bulge))
    ids, forces = pressure_face_load(nodes, _TET10_FACES, 5.0)
    np.testing.assert_allclose(np.cross(nodes[ids], forces).sum(axis=0), 0.0, atol=1e-9)


@settings(max_examples=50)
@given(bulges)
def test_reversed_face_gets_exactly_the_opposite_load(bulge: list) -> None:
    nodes = _tet10(_UNIT, np.array(bulge))
    face = _TET10_FACES[2:3]
    ids, forces = pressure_face_load(nodes, face, 3.0)
    ids_r, forces_r = pressure_face_load(nodes, face[:, _FLIP], 3.0)
    np.testing.assert_array_equal(ids, ids_r)
    np.testing.assert_allclose(forces_r, -forces, atol=1e-12)


@settings(max_examples=50)
@given(bulges)
def test_curved_tet10_volume_weights_give_its_centroid(bulge: list) -> None:
    # First moment of volume, from the flux of x^2 / 2 through the boundary (exact here).
    from tessera.physics.loads import _SHAPES, _W

    nodes = _tet10(_UNIT, np.array(bulge))
    w = tet10_volume_weights(nodes, np.arange(10)[None, :])[0]
    x = nodes[_TET10_FACES]
    flux = np.zeros(3)
    for (n, dxi, deta), wq in zip(_SHAPES, _W, strict=True):
        point = n @ x
        normal = np.cross(dxi @ x, deta @ x)
        flux += 0.5 * wq * (0.5 * point**2 * normal).sum(axis=0)
    # x^2 n has degree 6 on a curved face: compare with a tolerance, not exactly.
    np.testing.assert_allclose(w @ nodes, flux, rtol=1e-3, atol=1e-3)


def test_a_face_between_two_elements_has_no_outside() -> None:
    from tessera.model import InputError

    nodes = _tet10(_UNIT)
    mirror = nodes.copy()
    mirror[3] = [0.0, 0.0, -4.0]  # second element below the shared bottom face
    for k, (i, j) in enumerate(TET10_EDGES):
        mirror[4 + k] = 0.5 * (mirror[i] + mirror[j])
    both = np.vstack([nodes, mirror[[3, 7, 8, 9]]])
    second = np.array([0, 1, 2, 10, 4, 5, 6, 11, 12, 13])
    elements = np.vstack([np.arange(10), second])
    with pytest.raises(InputError, match="between two elements"):
        orient_outward(both, elements, _TET10_FACES[:1])
