"""Failure criteria and strength ratios (ADR-010).

The strength ratio R at a point is how many times every load could be multiplied before the
criterion predicts first failure there (linear analysis, so exact for the criterion); ``inf``
where the stress is zero. Orthotropic materials use Tsai-Wu on stresses in the material axes,
isotropic ones von Mises (ductile) or the maximum principal stress (brittle). Every function
here works on stresses alone; ``summarize`` decides where they are taken (ADR-010 D4).
"""

import numpy as np
from numpy.typing import NDArray

from tessera.model import (
    FailureCriterion,
    FailureMode,
    IsotropicStrength,
    Material,
    MaterialFrame,
    OrthotropicMaterial,
    OrthotropicStrength,
)
from tessera.model.elasticity import rotate_stress
from tessera.post.stress import von_mises

TSAI_WU_INTERACTION = -0.5
"""F*_ij in F_ij = F*_ij sqrt(F_ii F_jj) for the normal-stress pairs (ADR-010 D2, approved by
the maintainer): with equal strengths in every direction and S = T / sqrt(3), Tsai-Wu is then
exactly von Mises."""

_NORMAL_MODES = (
    (FailureMode.ALONG_RASTER_TENSION, FailureMode.ALONG_RASTER_COMPRESSION),
    (FailureMode.ACROSS_ROADS_TENSION, FailureMode.ACROSS_ROADS_COMPRESSION),
    (FailureMode.ACROSS_LAYERS_TENSION, FailureMode.ACROSS_LAYERS_COMPRESSION),
)
_SHEAR_MODES = (
    FailureMode.IN_LAYER_SHEAR,
    FailureMode.INTERLAYER_SHEAR,
    FailureMode.INTERLAYER_SHEAR,
)


def criterion(material: Material) -> FailureCriterion:
    """The criterion ADR-010 D2 assigns to ``material``; it must have a strength."""
    strength = material.strength
    if strength is None:
        raise ValueError(f"material {material.name!r} has no strength")
    if isinstance(strength, OrthotropicStrength):
        return FailureCriterion.TSAI_WU
    return FailureCriterion.MAX_PRINCIPAL if strength.brittle else FailureCriterion.VON_MISES


def _ratio(a: NDArray[np.float64], b: NDArray[np.float64]) -> NDArray[np.float64]:
    """Positive root R of a R^2 + b R = 1, written so a = 0 needs no special case: 1/b for
    b > 0, inf for b <= 0 (and a = b = 0, no stress)."""
    a = np.maximum(a, 0.0)  # a >= 0 in exact arithmetic (the quadratic form is semidefinite)
    with np.errstate(divide="ignore"):
        return 2.0 / (b + np.sqrt(b * b + 4.0 * a))


def tsai_wu_ratio(
    stress: NDArray[np.float64], strength: OrthotropicStrength
) -> NDArray[np.float64]:
    """R for stresses (..., 6) in the material axes, order 11 22 33 12 23 31 (MPa)."""
    t = np.array([strength.t1, strength.t2, strength.t3])
    c = np.array(strength.compressive)
    linear = 1.0 / t - 1.0 / c
    f = 1.0 / (t * c)
    normal = np.diag(f)
    for i, j in ((0, 1), (0, 2), (1, 2)):
        normal[i, j] = normal[j, i] = TSAI_WU_INTERACTION * np.sqrt(f[i] * f[j])
    s = np.asarray(stress, dtype=np.float64)
    sn = s[..., :3]
    shear = s[..., 3:] / np.array([strength.s12, strength.s23, strength.s13])
    a = np.einsum("...i,ij,...j->...", sn, normal, sn) + (shear**2).sum(axis=-1)
    return _ratio(a, sn @ linear)


def tsai_wu_mode(stress: NDArray[np.float64], strength: OrthotropicStrength) -> FailureMode:
    """The term with the largest maximum-stress ratio |sigma_k| / strength_k for one stress (6,)
    in the material axes: it names how the spot fails and never changes R (ADR-010 D2)."""
    s = np.asarray(stress, dtype=np.float64)
    t = (strength.t1, strength.t2, strength.t3)
    c = strength.compressive
    terms = [
        (s[k] / t[k] if s[k] >= 0 else -s[k] / c[k], _NORMAL_MODES[k][int(s[k] < 0)])
        for k in range(3)
    ]
    shear_strength = (strength.s12, strength.s23, strength.s13)
    terms += [(abs(s[3 + k]) / shear_strength[k], _SHEAR_MODES[k]) for k in range(3)]
    return max(terms, key=lambda term: term[0])[1]


def _principal(stress: NDArray[np.float64]) -> NDArray[np.float64]:
    """Principal stresses (..., 3), ascending, of Voigt stresses xx yy zz xy yz zx."""
    s = np.asarray(stress, dtype=np.float64)
    xx, yy, zz, xy, yz, zx = (s[..., k] for k in range(6))
    tensor = np.stack(
        [np.stack([xx, xy, zx], -1), np.stack([xy, yy, yz], -1), np.stack([zx, yz, zz], -1)], -2
    )
    return np.linalg.eigvalsh(tensor)


def isotropic_ratio(
    stress: NDArray[np.float64], strength: IsotropicStrength
) -> NDArray[np.float64]:
    """R for stresses (..., 6) in any axes: yield / von Mises for a ductile material, the lower
    of T / sigma_1 and C / -sigma_3 for a brittle one."""
    s = np.asarray(stress, dtype=np.float64)
    if not strength.brittle:
        with np.errstate(divide="ignore"):
            return strength.tensile / von_mises(s)
    principal = _principal(s)
    with np.errstate(divide="ignore"):
        tension = np.where(principal[..., 2] > 0, strength.tensile / principal[..., 2], np.inf)
        compression = np.where(
            principal[..., 0] < 0, strength.compressive_or_tensile / -principal[..., 0], np.inf
        )
    return np.minimum(tension, compression)


def isotropic_mode(stress: NDArray[np.float64], strength: IsotropicStrength) -> FailureMode:
    if not strength.brittle:
        return FailureMode.VON_MISES
    low, _, high = _principal(np.asarray(stress, dtype=np.float64))
    tension = strength.tensile / high if high > 0 else np.inf
    compression = strength.compressive_or_tensile / -low if low < 0 else np.inf
    if compression < tension:
        return FailureMode.MAX_PRINCIPAL_COMPRESSION
    return FailureMode.MAX_PRINCIPAL_TENSION


def strength_ratio(
    stress: NDArray[np.float64], material: Material, frame: MaterialFrame | None = None
) -> NDArray[np.float64]:
    """R for stresses (N, 6) in part axes (xx yy zz xy yz zx) of ``material``, whose axes are
    ``frame`` (``None``: the part axes)."""
    strength = material.strength
    if strength is None:
        raise ValueError(f"material {material.name!r} has no strength")
    if isinstance(material, OrthotropicMaterial):
        assert isinstance(strength, OrthotropicStrength)
        return tsai_wu_ratio(_material_axes(stress, frame), strength)
    assert isinstance(strength, IsotropicStrength)
    return isotropic_ratio(stress, strength)


def failure_mode(
    stress: NDArray[np.float64], material: Material, frame: MaterialFrame | None = None
) -> FailureMode:
    """How one stress (6,) in part axes fails ``material`` first (ADR-010 D2)."""
    strength = material.strength
    if strength is None:
        raise ValueError(f"material {material.name!r} has no strength")
    if isinstance(material, OrthotropicMaterial):
        assert isinstance(strength, OrthotropicStrength)
        return tsai_wu_mode(_material_axes(stress, frame), strength)
    assert isinstance(strength, IsotropicStrength)
    return isotropic_mode(stress, strength)


def _material_axes(stress: NDArray[np.float64], frame: MaterialFrame | None) -> NDArray[np.float64]:
    s = np.asarray(stress, dtype=np.float64)
    if frame is None:
        return s
    return rotate_stress(s.reshape(-1, 6), frame).reshape(s.shape)
