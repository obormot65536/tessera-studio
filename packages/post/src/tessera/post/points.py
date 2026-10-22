"""Stresses from integration points (ADR-008).

Works on ``PointStresses`` only and knows nothing about the solver that produced them. Every
reported peak stays nodal (ADR-008, decision 4): integration points lie inside the elements, so
a maximum over them reads low at surfaces. Do not build peaks on these functions without a new
ADR and the maintainer's approval (hard rule 6).
"""

import warnings

import numpy as np
from numpy.typing import NDArray

from tessera.model import Mesh, PointStresses
from tessera.post.summary import von_mises

_FLAT = 1e-8
"""Smallest to largest singular value of the (centred, scaled) fit matrix below which an
element's points count as lying in one plane."""


class DegenerateElementsWarning(UserWarning):
    """Some elements' integration points lie in one plane (or on a line), so no affine field
    passes through them; those elements got NaN. ``elements`` holds their indices."""

    def __init__(self, elements: NDArray[np.int64]) -> None:
        self.elements = elements
        super().__init__(
            f"{elements.size} element(s) have integration points in one plane and got NaN, "
            f"e.g. element {elements[0]}"
        )


def element_nodal_stresses(mesh: Mesh, point_stresses: PointStresses) -> NDArray[np.float64]:
    """(M, 10, 6) stresses at each element's ten nodes, from that element's points only (MPa).

    For each element, the affine field (constant plus linear in x, y, z) through its Q >= 4
    integration points, in physical coordinates (mm), is evaluated at its ten nodes: exactly
    through four points, by least squares through more. Elements are not averaged, so a node
    shared by two elements has two values, one per element, each in that element's material.
    The result is exact for a linear stress field. Components in the neutral order xx, yy, zz,
    xy, yz, zx. An element with a NaN or infinite point or stress gets NaN.

    An element whose points lie in one plane (a flat or collapsed element) gets NaN too, and a
    ``DegenerateElementsWarning`` lists such elements; the other elements are still fitted.
    Raises ``ValueError`` if ``point_stresses`` does not have one row per mesh element or has
    fewer than four points per element.
    """
    m, q = point_stresses.n_elements, point_stresses.n_points
    if m != mesh.n_elements:
        raise ValueError(f"point stresses are for {m} elements; the mesh has {mesh.n_elements}")
    if q < 4:
        raise ValueError(f"an affine fit needs at least 4 points per element, got {q}")
    out = np.full((m, mesh.elements.shape[1], 6), np.nan)
    ok = np.isfinite(point_stresses.points).all(axis=(1, 2))
    ok &= np.isfinite(point_stresses.stresses).all(axis=(1, 2))
    if not ok.any():
        return out

    # Centre and scale each element's points, so the fit stays well conditioned for a small
    # element far from the origin.
    points = point_stresses.points[ok]
    centre = points.mean(axis=1, keepdims=True)
    scale = np.abs(points - centre).max(axis=(1, 2), keepdims=True)
    scale[scale == 0] = 1.0
    fit = _affine_basis((points - centre) / scale)  # (K, Q, 4)
    singular = np.linalg.svd(fit, compute_uv=False)
    flat = singular[:, -1] <= _FLAT * singular[:, 0]
    fitted = np.flatnonzero(ok)
    if flat.any():
        warnings.warn(DegenerateElementsWarning(fitted[flat]), stacklevel=2)
        good = ~flat
        fitted, fit, centre, scale = fitted[good], fit[good], centre[good], scale[good]
    coefficients = np.linalg.pinv(fit) @ point_stresses.stresses[fitted]  # (K, 4, 6)
    nodes = mesh.nodes[mesh.elements[fitted]]
    out[fitted] = _affine_basis((nodes - centre) / scale) @ coefficients
    return out


def point_von_mises(point_stresses: PointStresses) -> NDArray[np.float64]:
    """(M, Q) von Mises stress at every integration point, MPa. For diagnostics, not peaks."""
    return von_mises(point_stresses.stresses)


def _affine_basis(x: NDArray[np.float64]) -> NDArray[np.float64]:
    """(..., P, 3) coordinates -> (..., P, 4) rows [1, x, y, z]."""
    return np.concatenate([np.ones((*x.shape[:-1], 1)), x], axis=-1)
