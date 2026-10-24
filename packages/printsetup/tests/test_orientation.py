"""Material frame from a print setup (ADR-009 D2)."""

import math

import numpy as np
import pytest

from tessera.model import MaterialFrame
from tessera.printsetup import material_frame


def _axes(frame: MaterialFrame) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    return np.array(frame.axis_1), np.array(frame.axis_2), np.array(frame.axis_3)


def _assert_frame(
    frame: MaterialFrame,
    e1: tuple[float, float, float],
    e2: tuple[float, float, float],
    e3: tuple[float, float, float],
    atol: float = 1e-12,
) -> None:
    a1, a2, a3 = _axes(frame)
    np.testing.assert_allclose(a1, e1, atol=atol)
    np.testing.assert_allclose(a2, e2, atol=atol)
    np.testing.assert_allclose(a3, e3, atol=atol)


def test_build_along_z_gives_identity() -> None:
    _assert_frame(material_frame((0.0, 0.0, 1.0)), (1, 0, 0), (0, 1, 0), (0, 0, 1))


def test_build_along_x_matches_adr_example() -> None:
    _assert_frame(material_frame((1.0, 0.0, 0.0)), (0, 0, -1), (0, 1, 0), (1, 0, 0))


def test_build_along_minus_z_is_half_turn_about_x() -> None:
    _assert_frame(material_frame((0.0, 0.0, -1.0)), (1, 0, 0), (0, -1, 0), (0, 0, -1))


def test_raster_90_about_z_puts_axis_1_along_y() -> None:
    _assert_frame(
        material_frame((0.0, 0.0, 1.0), raster_angle=90.0), (0, 1, 0), (-1, 0, 0), (0, 0, 1)
    )


@pytest.mark.parametrize("raster_angle", [0.0, 30.0, -135.0])
def test_arbitrary_tilt_is_orthonormal_right_handed_with_e3_along_build(
    raster_angle: float,
) -> None:
    b = np.array([0.3, -0.5, 0.8])
    b /= np.linalg.norm(b)
    frame = material_frame((float(b[0]), float(b[1]), float(b[2])), raster_angle=raster_angle)
    r = frame.rotation
    np.testing.assert_allclose(r @ r.T, np.eye(3), atol=1e-12)
    assert np.linalg.det(r) == pytest.approx(1.0, abs=1e-12)
    np.testing.assert_allclose(frame.axis_3, b, atol=1e-12)
    a1, _, _ = _axes(frame)
    a1_zero, _, _ = _axes(material_frame((float(b[0]), float(b[1]), float(b[2]))))
    assert float(a1 @ a1_zero) == pytest.approx(math.cos(math.radians(raster_angle)), abs=1e-12)


def test_minimal_tilt_keeps_axis_1_off_the_tilt_axis() -> None:
    # Tilting +Z toward +y rotates about x, so axis 1 stays on x (smallest rotation).
    c, s = math.cos(0.4), math.sin(0.4)
    _assert_frame(material_frame((0.0, -s, c)), (1, 0, 0), (0, c, s), (0, -s, c))


def test_non_unit_build_direction_is_normalized() -> None:
    a = material_frame((0.0, 0.0, 7.5))
    _assert_frame(a, (1, 0, 0), (0, 1, 0), (0, 0, 1))
    b = material_frame((2.0, 0.0, 0.0), raster_angle=20.0)
    c = material_frame((1.0, 0.0, 0.0), raster_angle=20.0)
    np.testing.assert_allclose(b.rotation, c.rotation, atol=1e-15)


@pytest.mark.parametrize(
    "build",
    [
        (0.0, 0.0, 0.0),
        (math.nan, 0.0, 1.0),
        (0.0, math.inf, 1.0),
    ],
)
def test_zero_or_non_finite_build_direction_raises(build: tuple[float, float, float]) -> None:
    with pytest.raises(ValueError, match="build direction"):
        material_frame(build)


@pytest.mark.parametrize("raster_angle", [math.nan, math.inf])
def test_non_finite_raster_angle_raises(raster_angle: float) -> None:
    with pytest.raises(ValueError, match="raster angle"):
        material_frame((0.0, 0.0, 1.0), raster_angle=raster_angle)


@pytest.mark.parametrize("eps", [1e-4, 1e-8, 1e-12, 1e-15])
def test_continuous_approaching_minus_z_along_the_x_half_turn(eps: float) -> None:
    # Approaching -Z by tilting about printer X (the path the -Z convention picks) converges
    # to the -Z frame: no division blow-up just off the antipode.
    exact = material_frame((0.0, 0.0, -1.0)).rotation
    for sign in (1.0, -1.0):
        near = material_frame((0.0, sign * math.sin(eps), -math.cos(eps))).rotation
        np.testing.assert_allclose(near, exact, atol=4 * eps + 1e-15)


@pytest.mark.parametrize("eps", [1e-6, 1e-10, 1e-14, 1e-200])
def test_just_off_minus_z_stays_orthonormal(eps: float) -> None:
    for build in ((eps, 0.0, -1.0), (eps, eps, -1.0), (-eps, 3 * eps, -1.0)):
        r = material_frame(build).rotation
        assert np.isfinite(r).all()
        np.testing.assert_allclose(r @ r.T, np.eye(3), atol=1e-12)
        assert np.linalg.det(r) == pytest.approx(1.0, abs=1e-12)
