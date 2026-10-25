"""Failure criteria and strength ratios (ADR-010 D2), on stresses alone."""

import math

import numpy as np
import pytest
from numpy.typing import NDArray

from tessera.model import (
    FailureCriterion,
    FailureMode,
    IsotropicMaterial,
    IsotropicStrength,
    MaterialFrame,
    OrthotropicMaterial,
    OrthotropicStrength,
)
from tessera.post import (
    TSAI_WU_INTERACTION,
    criterion,
    failure_mode,
    isotropic_ratio,
    strength_ratio,
    tsai_wu_ratio,
    von_mises,
)

# Made-up TO strengths (as the catalog's TO stiffnesses, not data for a real filament).
TO_STRENGTH = OrthotropicStrength(
    t1=40.0, t2=30.0, t3=15.0, s12=20.0, s13=12.0, s23=12.0, c1=50.0, c2=45.0, c3=40.0
)
TO_STIFFNESS = dict(
    e1=3000.0, e2=2400.0, e3=1800.0, nu12=0.35, nu13=0.30, nu23=0.30, g12=1000.0, g13=700.0,
    g23=650.0,
)  # fmt: skip


def voigt(**components: float) -> NDArray[np.float64]:
    """A stress (6,) in 11 22 33 12 23 31 (or xx yy zz xy yz zx) from named components."""
    order = ("s11", "s22", "s33", "s12", "s23", "s31")
    return np.array([components.get(name, 0.0) for name in order])


def test_interaction_coefficient_is_the_approved_one() -> None:
    assert TSAI_WU_INTERACTION == -0.5


@pytest.mark.parametrize(
    ("component", "strength"),
    [
        ("s11", 40.0),
        ("s22", 30.0),
        ("s33", 15.0),
    ],
)
@pytest.mark.parametrize("sigma", [5.0, 12.5])
def test_uniaxial_tension_along_a_material_axis(
    component: str, strength: float, sigma: float
) -> None:
    ratio = tsai_wu_ratio(voigt(**{component: sigma}), TO_STRENGTH)
    assert ratio == pytest.approx(strength / sigma, rel=1e-12)


@pytest.mark.parametrize(("component", "strength"), [("s11", 50.0), ("s22", 45.0), ("s33", 40.0)])
def test_uniaxial_compression_along_a_material_axis(component: str, strength: float) -> None:
    ratio = tsai_wu_ratio(voigt(**{component: -10.0}), TO_STRENGTH)
    assert ratio == pytest.approx(strength / 10.0, rel=1e-12)


@pytest.mark.parametrize(("component", "strength"), [("s12", 20.0), ("s23", 12.0), ("s31", 12.0)])
@pytest.mark.parametrize("tau", [4.0, -4.0])
def test_pure_shear(component: str, strength: float, tau: float) -> None:
    ratio = tsai_wu_ratio(voigt(**{component: tau}), TO_STRENGTH)
    assert ratio == pytest.approx(strength / abs(tau), rel=1e-12)


def test_ratio_is_the_load_factor_that_reaches_the_envelope() -> None:
    rng = np.random.default_rng(10)
    stress = rng.normal(scale=8.0, size=(200, 6))
    ratio = tsai_wu_ratio(stress, TO_STRENGTH)
    # Tsai-Wu index at R sigma is exactly 1.
    t = np.array([40.0, 30.0, 15.0])
    c = np.array([50.0, 45.0, 40.0])
    f = 1.0 / (t * c)
    q = np.diag(f)
    for i, j in ((0, 1), (0, 2), (1, 2)):
        q[i, j] = q[j, i] = -0.5 * math.sqrt(f[i] * f[j])
    s = ratio[:, None] * stress
    index = (
        s[:, :3] @ (1 / t - 1 / c)
        + np.einsum("ni,ij,nj->n", s[:, :3], q, s[:, :3])
        + (s[:, 3] / 20.0) ** 2
        + (s[:, 4] / 12.0) ** 2
        + (s[:, 5] / 12.0) ** 2
    )
    assert np.allclose(index, 1.0, rtol=0, atol=1e-12)


def test_ratio_scales_inversely_with_the_load() -> None:
    """R(k sigma) = R(sigma) / k for k > 0, linear terms included: the ratio is a load factor,
    unlike the Tsai-Wu index, which is not proportional to the load."""
    sigma = voigt(s11=10.0, s33=-5.0, s12=3.0)
    one, two, half = tsai_wu_ratio(np.stack([sigma, 2 * sigma, 0.5 * sigma]), TO_STRENGTH)
    assert (two, half) == (pytest.approx(one / 2, rel=1e-12), pytest.approx(2 * one, rel=1e-12))
    assert tsai_wu_ratio(-sigma, TO_STRENGTH) != pytest.approx(one, rel=1e-6)


def test_equal_strengths_reduce_tsai_wu_to_von_mises() -> None:
    x = 100.0
    s = x / math.sqrt(3.0)
    iso = OrthotropicStrength(t1=x, t2=x, t3=x, s12=s, s13=s, s23=s)
    stress = np.random.default_rng(3).normal(scale=30.0, size=(500, 6))
    assert np.allclose(tsai_wu_ratio(stress, iso), x / von_mises(stress), rtol=1e-12, atol=0)


def test_hydrostatic_stress_with_equal_strengths_never_fails() -> None:
    x = 100.0
    s = x / math.sqrt(3.0)
    iso = OrthotropicStrength(t1=x, t2=x, t3=x, s12=s, s13=s, s23=s)
    assert np.isinf(tsai_wu_ratio(voigt(s11=50.0, s22=50.0, s33=50.0), iso))


def test_zero_stress_has_an_infinite_ratio() -> None:
    assert np.isinf(tsai_wu_ratio(np.zeros(6), TO_STRENGTH))
    assert np.isinf(isotropic_ratio(np.zeros(6), IsotropicStrength(250.0)))
    assert np.isinf(isotropic_ratio(np.zeros(6), IsotropicStrength(50.0, brittle=True)))


def test_missing_compressive_strengths_are_the_tensile_ones() -> None:
    tension_only = OrthotropicStrength(t1=40.0, t2=30.0, t3=15.0, s12=20.0, s13=12.0, s23=12.0)
    assert tsai_wu_ratio(voigt(s33=-10.0), tension_only) == pytest.approx(1.5, rel=1e-12)


def test_ductile_isotropic_is_yield_over_von_mises() -> None:
    stress = np.random.default_rng(4).normal(scale=50.0, size=(50, 6))
    ratio = isotropic_ratio(stress, IsotropicStrength(250.0))
    assert np.allclose(ratio, 250.0 / von_mises(stress), rtol=1e-14, atol=0)


def test_brittle_isotropic_uses_the_principal_stresses() -> None:
    resin = IsotropicStrength(50.0, compressive=80.0, brittle=True)
    # Pure shear 20 MPa: principal stresses +20 and -20, tension governs.
    assert isotropic_ratio(voigt(s12=20.0), resin) == pytest.approx(2.5, rel=1e-12)
    # Uniaxial compression 40 MPa: 80 / 40.
    assert isotropic_ratio(voigt(s22=-40.0), resin) == pytest.approx(2.0, rel=1e-12)
    # Missing compressive strength: the tensile one.
    assert isotropic_ratio(
        voigt(s22=-40.0), IsotropicStrength(50.0, brittle=True)
    ) == pytest.approx(1.25, rel=1e-12)


@pytest.mark.parametrize(
    ("stress", "mode"),
    [
        (voigt(s11=30.0), FailureMode.ALONG_RASTER_TENSION),
        (voigt(s11=-30.0), FailureMode.ALONG_RASTER_COMPRESSION),
        (voigt(s22=20.0), FailureMode.ACROSS_ROADS_TENSION),
        (voigt(s22=-20.0), FailureMode.ACROSS_ROADS_COMPRESSION),
        (voigt(s11=20.0, s33=10.0), FailureMode.ACROSS_LAYERS_TENSION),
        (voigt(s33=-20.0), FailureMode.ACROSS_LAYERS_COMPRESSION),
        (voigt(s12=10.0), FailureMode.IN_LAYER_SHEAR),
        (voigt(s23=-10.0), FailureMode.INTERLAYER_SHEAR),
        (voigt(s31=10.0, s11=20.0), FailureMode.INTERLAYER_SHEAR),
    ],
)
def test_tsai_wu_mode_names_the_largest_maximum_stress_term(
    stress: NDArray[np.float64], mode: FailureMode
) -> None:
    material = OrthotropicMaterial("to", **TO_STIFFNESS, strength=TO_STRENGTH)
    assert failure_mode(stress, material) is mode


def test_isotropic_modes() -> None:
    steel = IsotropicMaterial("steel", 210000.0, 0.3, strength=IsotropicStrength(250.0))
    resin = IsotropicMaterial(
        "resin", 2500.0, 0.35, strength=IsotropicStrength(50.0, compressive=80.0, brittle=True)
    )
    assert failure_mode(voigt(s11=-100.0), steel) is FailureMode.VON_MISES
    assert failure_mode(voigt(s11=10.0), resin) is FailureMode.MAX_PRINCIPAL_TENSION
    assert failure_mode(voigt(s11=-100.0), resin) is FailureMode.MAX_PRINCIPAL_COMPRESSION


def test_criterion_by_material() -> None:
    assert (
        criterion(OrthotropicMaterial("to", **TO_STIFFNESS, strength=TO_STRENGTH))
        is FailureCriterion.TSAI_WU
    )
    steel = IsotropicMaterial("s", 210000.0, 0.3, strength=IsotropicStrength(250.0))
    resin = IsotropicMaterial("r", 2500.0, 0.35, strength=IsotropicStrength(50.0, brittle=True))
    assert criterion(steel) is FailureCriterion.VON_MISES
    assert criterion(resin) is FailureCriterion.MAX_PRINCIPAL
    with pytest.raises(ValueError, match="material 'x' has no strength"):
        criterion(IsotropicMaterial("x", 1.0, 0.3))


def test_strength_ratio_turns_part_stresses_into_the_material_axes() -> None:
    """Printed on its side (build along x): a pull along x pulls the layers apart."""
    material = OrthotropicMaterial("to", **TO_STIFFNESS, strength=TO_STRENGTH)
    side = MaterialFrame((0.0, 0.0, -1.0), (0.0, 1.0, 0.0))  # axis 3 = +x
    pull_x = voigt(s11=10.0)[None, :]
    assert strength_ratio(pull_x, material, side)[0] == pytest.approx(1.5, rel=1e-12)
    assert strength_ratio(pull_x, material)[0] == pytest.approx(4.0, rel=1e-12)
    assert failure_mode(pull_x[0], material, side) is FailureMode.ACROSS_LAYERS_TENSION


@pytest.mark.parametrize(
    ("angle", "expected"),
    [(0.0, 4.0), (30.0, 3.380006), (45.0, 3.096069), (60.0, 2.982983), (90.0, 3.0)],
)
def test_off_axis_values_in_the_catalog(angle: float, expected: float) -> None:
    """TO-06's closed-form values: raster turned by ``angle`` about Z, 10 MPa along x."""
    t = math.radians(angle)
    frame = MaterialFrame((math.cos(t), math.sin(t), 0.0), (-math.sin(t), math.cos(t), 0.0))
    material = OrthotropicMaterial("to", **TO_STIFFNESS, strength=TO_STRENGTH)
    ratio = strength_ratio(voigt(s11=10.0)[None, :], material, frame)[0]
    assert ratio == pytest.approx(expected, abs=5e-7)
