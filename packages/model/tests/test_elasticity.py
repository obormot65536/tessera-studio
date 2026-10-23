"""Orthotropic elasticity helpers (ADR-009 D1, D2)."""

import math

import numpy as np
import pytest

from tessera.model.elasticity import check_orthotropic, orthotropic_compliance, rotate_stress
from tessera.model.solver import MaterialFrame

# The catalog's TO constants (MPa): 1 = raster, 2 = across the road, 3 = build.
TO = {
    "e1": 3000.0,
    "e2": 2400.0,
    "e3": 1800.0,
    "nu12": 0.35,
    "nu13": 0.30,
    "nu23": 0.30,
    "g12": 1000.0,
    "g13": 700.0,
    "g23": 650.0,
}


def test_compliance_matches_hand_built_matrix() -> None:
    e1, e2, e3 = 3000.0, 2400.0, 1800.0
    expected = np.zeros((6, 6))
    expected[0, 0], expected[1, 1], expected[2, 2] = 1 / e1, 1 / e2, 1 / e3
    expected[0, 1] = expected[1, 0] = -0.35 / e1
    expected[0, 2] = expected[2, 0] = -0.30 / e1
    expected[1, 2] = expected[2, 1] = -0.30 / e2
    # Voigt xx yy zz xy yz zx with engineering shear strains.
    expected[3, 3], expected[4, 4], expected[5, 5] = 1 / 1000.0, 1 / 650.0, 1 / 700.0

    s = orthotropic_compliance(**TO)

    assert s.shape == (6, 6)
    np.testing.assert_allclose(s, expected, rtol=1e-15, atol=0)
    np.testing.assert_array_equal(s, s.T)


def test_compliance_reciprocity_and_minor_ratios() -> None:
    s = orthotropic_compliance(**TO)
    # nu_ji = -S_ij E_j: the derived minor ratios of the catalog (nu21, nu31, nu32).
    assert -s[0, 1] * 2400.0 == pytest.approx(0.28, rel=1e-12)
    assert -s[0, 2] * 1800.0 == pytest.approx(0.18, rel=1e-12)
    assert -s[1, 2] * 1800.0 == pytest.approx(0.225, rel=1e-12)


def test_catalog_set_is_accepted_with_delta() -> None:
    check_orthotropic(**TO)
    nu21, nu31, nu32 = 0.28, 0.18, 0.225
    delta = 1 - 0.35 * nu21 - 0.30 * nu32 - 0.30 * nu31 - 2 * nu21 * nu32 * 0.30
    assert delta == pytest.approx(0.7427, abs=1e-12)
    # det(S) = Delta / (E1 E2 E3 G12 G13 G23): the compliance agrees with the closed form.
    s = orthotropic_compliance(**TO)
    expected_det = delta / (3000.0 * 2400.0 * 1800.0 * 1000.0 * 700.0 * 650.0)
    assert np.linalg.det(s) == pytest.approx(expected_det, rel=1e-10)


def test_isotropic_limit_inverts_to_lame_constants() -> None:
    e, nu = 2000.0, 0.3
    g = e / (2 * (1 + nu))
    lam = e * nu / ((1 + nu) * (1 - 2 * nu))
    s = orthotropic_compliance(e, e, e, nu, nu, nu, g, g, g)

    c = np.linalg.inv(s)

    expected = np.zeros((6, 6))
    expected[:3, :3] = lam
    expected[np.arange(3), np.arange(3)] = lam + 2 * g
    expected[np.arange(3, 6), np.arange(3, 6)] = g
    np.testing.assert_allclose(c, expected, rtol=1e-12, atol=1e-9)


def _with(**changes: float) -> dict[str, float]:
    return {**TO, **changes}


@pytest.mark.parametrize(
    "constants",
    [
        pytest.param(_with(e2=3000.0, nu12=1.2), id="nu12-above-sqrt-bound"),
        pytest.param(
            _with(e1=2000.0, e2=2000.0, e3=2000.0, nu12=0.6, nu13=0.6, nu23=0.6),
            id="delta-negative",
        ),
    ],
)
def test_not_positive_definite_is_rejected_with_bounds_in_message(
    constants: dict[str, float],
) -> None:
    with pytest.raises(ValueError, match="positive definite") as info:
        check_orthotropic(**constants)
    message = str(info.value)
    assert "|nu_ij| < sqrt(E_i/E_j)" in message
    assert "Delta" in message
    assert "> 0" in message


@pytest.mark.parametrize(
    ("name", "value"),
    [
        ("e1", 0.0),
        ("e2", -2400.0),
        ("e3", math.inf),
        ("g12", 0.0),
        ("g13", -1.0),
        ("g23", math.nan),
        ("e1", math.nan),
    ],
)
def test_non_positive_or_non_finite_modulus_is_rejected(name: str, value: float) -> None:
    with pytest.raises(ValueError, match=rf"{name}\b.*positive.*finite"):
        check_orthotropic(**_with(**{name: value}))


@pytest.mark.parametrize("name", ["nu12", "nu13", "nu23"])
def test_non_finite_poisson_ratio_is_rejected(name: str) -> None:
    with pytest.raises(ValueError, match=rf"{name}\b.*finite"):
        check_orthotropic(**_with(**{name: math.nan}))


def _frame_about_z(degrees: float) -> MaterialFrame:
    a = math.radians(degrees)
    return MaterialFrame.from_axes((math.cos(a), math.sin(a), 0.0), (-math.sin(a), math.cos(a), 0))


def test_rotate_stress_identity_frame_is_unchanged() -> None:
    stress = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
    frame = MaterialFrame((1.0, 0.0, 0.0), (0.0, 1.0, 0.0))
    np.testing.assert_array_equal(rotate_stress(stress, frame), stress)


def test_rotate_stress_uniaxial_at_45_degrees() -> None:
    stress = np.array([10.0, 0, 0, 0, 0, 0])
    local = rotate_stress(stress, _frame_about_z(45.0))
    np.testing.assert_allclose(local, [5.0, 5.0, 0, -5.0, 0, 0], atol=1e-12)


def test_rotate_stress_build_along_x_permutes_components() -> None:
    # ADR-009's example: build (1, 0, 0) gives e1 = -z, e2 = y, e3 = x.
    frame = MaterialFrame((0.0, 0.0, -1.0), (0.0, 1.0, 0.0))
    stress = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])  # xx yy zz xy yz zx
    local = rotate_stress(stress, frame)
    # s11 = zz, s22 = yy, s33 = xx, s12 = -zy, s23 = yx, s31 = -xz.
    np.testing.assert_allclose(local, [3.0, 2.0, 1.0, -5.0, 4.0, -6.0], atol=1e-15)


def test_rotate_stress_round_trips_and_keeps_invariants() -> None:
    rng = np.random.default_rng(9)
    stress = rng.normal(size=(5, 6)) * 10
    frame = MaterialFrame.from_axes((1.0, 2.0, -0.5), (0.3, -1.0, 2.0))
    back = MaterialFrame(tuple(frame.rotation[:, 0]), tuple(frame.rotation[:, 1]))  # type: ignore[arg-type]

    local = rotate_stress(stress, frame)

    assert local.shape == (5, 6)
    np.testing.assert_allclose(rotate_stress(local, back), stress, atol=1e-12)
    np.testing.assert_allclose(local[:, :3].sum(axis=1), stress[:, :3].sum(axis=1), atol=1e-12)


def test_rotate_stress_leaves_isotropic_state_unchanged() -> None:
    stress = np.array([-7.0, -7.0, -7.0, 0.0, 0.0, 0.0])
    frame = MaterialFrame.from_axes((1.0, 2.0, -0.5), (0.3, -1.0, 2.0))
    np.testing.assert_allclose(rotate_stress(stress, frame), stress, atol=1e-12)


def test_rotate_stress_rejects_wrong_shape() -> None:
    frame = MaterialFrame((1.0, 0.0, 0.0), (0.0, 1.0, 0.0))
    with pytest.raises(ValueError, match="6"):
        rotate_stress(np.zeros(5), frame)
