"""Moment loads spread over faces like a distributing coupling."""

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from tessera.model import TRI6_EDGES, InputError
from tessera.physics import moment_face_load


def _disc(radius: float = 5.0, n: int = 24) -> tuple[np.ndarray, np.ndarray]:
    """A flat fan of tri6 around the origin in z = 0 (straight edges)."""
    rim = [
        (radius * np.cos(t), radius * np.sin(t), 0.0)
        for t in np.linspace(0, 2 * np.pi, n, endpoint=False)
    ]
    nodes = [(0.0, 0.0, 0.0), *rim]
    tris = []
    for k in range(n):
        corners = [0, 1 + k, 1 + (k + 1) % n]
        mids = []
        for i, j in TRI6_EDGES:
            nodes.append(tuple(0.5 * (np.array(nodes[corners[i]]) + np.array(nodes[corners[j]]))))
            mids.append(len(nodes) - 1)
        tris.append(corners + mids)
    return np.array(nodes), np.array(tris)


coord = st.floats(min_value=-10, max_value=10, allow_nan=False)
vec = st.tuples(coord, coord, coord)


@settings(max_examples=50)
@given(vec, vec)
def test_forces_have_no_resultant_and_exactly_the_moment(moment: tuple, shift: tuple) -> None:
    nodes, tris = _disc()
    nodes = nodes + np.array(shift)  # the centroid moves with the face
    ids, forces = moment_face_load(nodes, tris, np.array(moment))
    np.testing.assert_allclose(forces.sum(axis=0), 0.0, atol=1e-9)
    centroid = np.array(shift)
    lever = nodes[ids] - centroid
    np.testing.assert_allclose(np.cross(lever, forces).sum(axis=0), moment, atol=1e-9)


def test_torque_on_a_disc_is_tangential_and_grows_with_radius() -> None:
    nodes, tris = _disc()
    ids, forces = moment_face_load(nodes, tris, np.array([0.0, 0.0, 100.0]))
    r = nodes[ids]
    radial = np.einsum("ij,ij->i", forces, r)
    np.testing.assert_allclose(radial, 0.0, atol=1e-12)
    np.testing.assert_allclose(forces[:, 2], 0.0, atol=1e-12)
    assert np.all(np.cross(r, forces)[:, 2] >= -1e-12)  # every node turns the same way


def test_a_line_cannot_carry_a_moment() -> None:
    nodes = np.array([[0.0, 0, 0], [1, 0, 0], [2, 0, 0], [0.5, 0, 0], [1.5, 0, 0], [1, 0, 0]])
    with pytest.raises(InputError, match="span an area"):
        moment_face_load(nodes, np.array([[0, 1, 2, 3, 4, 5]]), np.array([0.0, 0, 1]))


def _strip(length: float, width: float, origin: tuple[float, float, float]) -> tuple:
    """Two flat tri6 forming a length x width rectangle in z = origin z."""
    corners = np.array([[0, 0, 0], [length, 0, 0], [length, width, 0], [0, width, 0]], float)
    corners += origin
    nodes = list(corners)
    tris = []
    for a, b, c in ((0, 1, 2), (0, 2, 3)):
        mids = []
        for i, j in ((a, b), (b, c), (c, a)):
            nodes.append(0.5 * (corners[i] + corners[j]))
            mids.append(len(nodes) - 1)
        tris.append([a, b, c, *mids])
    return np.array(nodes), np.array(tris)


@pytest.mark.parametrize("origin", [(0.0, 0.0, 0.0), (1e4, -1e4, 5e3), (1e5, 1e5, 1e5)])
def test_thin_strip_far_from_the_origin_keeps_zero_net_force(origin: tuple) -> None:
    nodes, tris = _strip(200.0, 0.01, origin)
    _, forces = moment_face_load(nodes, tris, np.array([1e3, 1e3, 1e3]))
    assert np.linalg.norm(forces.sum(axis=0)) <= 1e-12 * np.abs(forces).max()


def test_disjoint_faces_share_one_moment_about_their_common_centroid() -> None:
    a_nodes, a_tris = _strip(10.0, 10.0, (0.0, 0.0, 0.0))
    b_nodes, b_tris = _strip(10.0, 10.0, (30.0, 0.0, 0.0))
    nodes = np.vstack([a_nodes, b_nodes])
    tris = np.vstack([a_tris, b_tris + len(a_nodes)])
    ids, forces = moment_face_load(nodes, tris, np.array([0.0, 0.0, 500.0]))
    np.testing.assert_allclose(forces.sum(axis=0), 0.0, atol=1e-10)
    lever = nodes[ids] - [20.0, 5.0, 0.0]
    np.testing.assert_allclose(np.cross(lever, forces).sum(axis=0), [0, 0, 500.0], atol=1e-9)


def test_a_needle_too_thin_to_carry_a_moment_is_refused() -> None:
    nodes, tris = _strip(1.0, 1e-7, (0.0, 0.0, 0.0))
    with pytest.raises(InputError, match="span an area"):
        moment_face_load(nodes, tris, np.array([1.0, 0.0, 0.0]))
