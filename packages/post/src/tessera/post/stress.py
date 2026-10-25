"""Scalar measures of a stress tensor."""

import numpy as np
from numpy.typing import NDArray


def von_mises(stresses: NDArray[np.float64]) -> NDArray[np.float64]:
    """Von Mises stress from (..., 6) components xx, yy, zz, xy, yz, zx (MPa), shape (...).

    (N, 6) nodal stresses give (N,); (M, Q, 6) integration-point stresses give (M, Q).
    """
    if stresses.shape[-1:] != (6,):
        raise ValueError(f"stresses must have shape (..., 6), got {stresses.shape}")
    sxx, syy, szz, sxy, syz, szx = (stresses[..., k] for k in range(6))
    return np.sqrt(
        0.5 * ((sxx - syy) ** 2 + (syy - szz) ** 2 + (szz - sxx) ** 2)
        + 3.0 * (sxy**2 + syz**2 + szx**2)
    )
