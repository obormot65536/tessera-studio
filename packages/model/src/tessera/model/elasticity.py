"""Orthotropic linear elasticity (ADR-009 D1).

Voigt order everywhere: xx, yy, zz, xy, yz, zx (the order of ``SolveResult.stresses``), with
engineering shear strains (gamma = 2 eps). Material axes: 1 = raster, 2 = across the road in
the layer, 3 = build. Poisson ratios are the major ones, nu_ij = -eps_j / eps_i under uniaxial
sigma_i; the minor ones follow from nu_ji = nu_ij E_j / E_i. Units: MPa (ADR-002).
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np
from numpy.typing import ArrayLike, NDArray

if TYPE_CHECKING:
    from tessera.model.solver import MaterialFrame

_BOUNDS = (
    "the usual bounds are |nu_ij| < sqrt(E_i/E_j) and "
    "Delta = 1 - nu12 nu21 - nu23 nu32 - nu13 nu31 - 2 nu21 nu32 nu13 > 0"
)


def orthotropic_compliance(
    e1: float,
    e2: float,
    e3: float,
    nu12: float,
    nu13: float,
    nu23: float,
    g12: float,
    g13: float,
    g23: float,
) -> NDArray[np.float64]:
    """Return the 6 x 6 compliance S (1/MPa) with eps = S sigma, in Voigt xx yy zz xy yz zx.

    No validation; see ``check_orthotropic``.
    """
    s = np.zeros((6, 6), dtype=np.float64)
    s[0, 0], s[1, 1], s[2, 2] = 1.0 / e1, 1.0 / e2, 1.0 / e3
    s[0, 1] = s[1, 0] = -nu12 / e1
    s[0, 2] = s[2, 0] = -nu13 / e1
    s[1, 2] = s[2, 1] = -nu23 / e2
    s[3, 3], s[4, 4], s[5, 5] = 1.0 / g12, 1.0 / g23, 1.0 / g13
    return s


def check_orthotropic(
    e1: float,
    e2: float,
    e3: float,
    nu12: float,
    nu13: float,
    nu23: float,
    g12: float,
    g13: float,
    g23: float,
) -> None:
    """Raise ``ValueError`` unless the constants describe a stable orthotropic material.

    Every modulus must be finite and greater than 0, every Poisson ratio finite, and the
    compliance positive definite (smallest eigenvalue greater than 0).
    """
    moduli = {"e1": e1, "e2": e2, "e3": e3, "g12": g12, "g13": g13, "g23": g23}
    for name, value in moduli.items():
        if not (math.isfinite(value) and value > 0):
            raise ValueError(f"{name} must be positive and finite (MPa), got {value!r}")
    for name, value in {"nu12": nu12, "nu13": nu13, "nu23": nu23}.items():
        if not math.isfinite(value):
            raise ValueError(f"{name} must be finite, got {value!r}")
    s = orthotropic_compliance(e1, e2, e3, nu12, nu13, nu23, g12, g13, g23)
    if np.linalg.eigvalsh(s).min() <= 0:
        raise ValueError(
            "orthotropic constants do not give a positive definite compliance "
            f"(nu12={nu12!r}, nu13={nu13!r}, nu23={nu23!r}); {_BOUNDS}, "
            "with major ratios nu_ij and nu_ji = nu_ij E_j / E_i"
        )


def rotate_stress(stress: ArrayLike, frame: MaterialFrame) -> NDArray[np.float64]:
    """Express Voigt stresses given in part axes in the material axes of ``frame``.

    ``stress`` is (6,) or (N, 6) in xx yy zz xy yz zx; the result has the same shape, in
    11 22 33 12 23 31: sigma' = R sigma R^T with R = ``frame.rotation``.
    """
    sig = np.asarray(stress, dtype=np.float64)
    if sig.shape[-1:] != (6,) or sig.ndim not in (1, 2):
        raise ValueError(f"stress must have shape (6,) or (N, 6), got {sig.shape}")
    flat = sig.reshape(-1, 6)
    xx, yy, zz, xy, yz, zx = flat.T
    tensor = np.stack(
        [
            np.stack([xx, xy, zx], axis=-1),
            np.stack([xy, yy, yz], axis=-1),
            np.stack([zx, yz, zz], axis=-1),
        ],
        axis=-2,
    )
    r = frame.rotation
    local = r @ tensor @ r.T
    out = np.stack(
        [
            local[:, 0, 0],
            local[:, 1, 1],
            local[:, 2, 2],
            local[:, 0, 1],
            local[:, 1, 2],
            local[:, 2, 0],
        ],
        axis=-1,
    )
    return out.reshape(sig.shape)
