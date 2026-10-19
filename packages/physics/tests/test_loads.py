import numpy as np
import pytest
from hypothesis import given
from hypothesis import strategies as st

from tessera.model import TRI6_EDGES
from tessera.physics import tri6_load_weights, uniform_face_load


def _flat_tri6(p0, p1, p2) -> np.ndarray:
    corners = np.array([p0, p1, p2], dtype=float)
    mids = [0.5 * (corners[i] + corners[j]) for i, j in TRI6_EDGES]
    return np.vstack([corners, mids])


coord = st.floats(min_value=-100, max_value=100, allow_nan=False)
point = st.tuples(coord, coord, coord)


@given(point, point, point)
def test_flat_tri6_weights_are_zero_at_corners_and_a_third_at_midsides(p0, p1, p2) -> None:
    nodes = _flat_tri6(p0, p1, p2)
    area = 0.5 * np.linalg.norm(np.cross(nodes[1] - nodes[0], nodes[2] - nodes[0]))
    w = tri6_load_weights(nodes, np.arange(6)[None, :])[0]
    np.testing.assert_allclose(w[:3], 0.0, atol=1e-9 * max(area, 1.0))
    np.testing.assert_allclose(w[3:], area / 3, rtol=1e-9, atol=1e-9)


def test_uniform_load_sums_to_total_force() -> None:
    nodes = np.vstack([_flat_tri6((0, 0, 0), (2, 0, 0), (0, 1, 0)), [[2.0, 1.0, 0.0]]])
    # Second triangle shares edge (1)-(2) and its midside node (index 4).
    mid_a = 0.5 * (nodes[1] + nodes[6])
    mid_b = 0.5 * (nodes[6] + nodes[2])
    nodes = np.vstack([nodes, mid_a, mid_b])
    tris = np.array([[0, 1, 2, 3, 4, 5], [1, 6, 2, 7, 8, 4]])
    ids, forces = uniform_face_load(nodes, tris, np.array([0.0, 0.0, -100.0]))
    np.testing.assert_allclose(forces.sum(axis=0), [0.0, 0.0, -100.0])
    assert set(ids.tolist()) == set(range(9))
    assert forces[:, 2].min() == pytest.approx(-100.0 / 3)  # shared midside node gets two thirds
