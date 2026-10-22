"""Integration-point stresses in the solver result (ADR-008)."""

import dataclasses

import numpy as np
import pytest

import tessera.model as model
from tessera.model import PointStresses, SolveResult


def point_stresses(m: int = 3, q: int = 4) -> PointStresses:
    rng = np.random.default_rng(8)
    return PointStresses(points=rng.random((m, q, 3)), stresses=rng.random((m, q, 6)))


def test_shapes_and_counts() -> None:
    ps = point_stresses(3, 4)
    assert ps.points.shape == (3, 4, 3)
    assert ps.stresses.shape == (3, 4, 6)
    assert ps.n_elements == 3
    assert ps.n_points == 4


def test_no_elements_is_allowed() -> None:
    ps = PointStresses(points=np.zeros((0, 4, 3)), stresses=np.zeros((0, 4, 6)))
    assert ps.n_elements == 0


@pytest.mark.parametrize(
    ("points", "stresses"),
    [
        ((3, 4, 2), (3, 4, 6)),  # points not 3-D vectors
        ((3, 4, 3), (3, 4, 5)),  # stresses not 6 components
        ((12, 3), (12, 6)),  # flattened, no point axis
        ((3, 4, 3), (2, 4, 6)),  # different element counts
        ((3, 4, 3), (3, 1, 6)),  # different point counts
        ((3, 4, 3, 1), (3, 4, 6)),  # too many axes
    ],
)
def test_wrong_shapes_are_rejected(points: tuple[int, ...], stresses: tuple[int, ...]) -> None:
    with pytest.raises(ValueError, match="shape"):
        PointStresses(points=np.zeros(points), stresses=np.zeros(stresses))


def test_point_stresses_are_frozen() -> None:
    ps = point_stresses()
    with pytest.raises(dataclasses.FrozenInstanceError):
        ps.points = np.zeros((3, 4, 3))  # type: ignore[misc]


def test_solve_result_has_no_point_stresses_by_default() -> None:
    result = SolveResult(np.zeros((2, 3)), np.zeros((2, 6)), "synthetic")
    assert result.point_stresses is None


def test_point_stresses_is_an_additive_last_field() -> None:
    # ADR-004: SolveResult grows additively, so existing positional calls keep working.
    names = [f.name for f in dataclasses.fields(SolveResult)]
    assert names == ["displacements", "stresses", "solver", "log", "metadata", "point_stresses"]
    ps = point_stresses()
    result = SolveResult(np.zeros((2, 3)), np.zeros((2, 6)), "synthetic", point_stresses=ps)
    assert result.point_stresses is ps


def test_exported_from_the_package() -> None:
    assert "PointStresses" in model.__all__
    assert model.PointStresses is PointStresses
