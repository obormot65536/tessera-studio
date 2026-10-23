"""Thin ring under diametral compression, shortening of the loaded diameter (TS-10; gate 3 %).

Problem (mm, N, MPa per ADR-002)
--------------------------------
Circular steel ring in the xy plane: mean radius R = 50, radial thickness h = 5 (R/h = 10),
axial width b = 5 (z), so a square section; inner radius R_i = 47.5, outer R_o = 52.5.
E = 210 000 MPa, nu = 0.30 (``steel-generic``). Two opposite forces P = 100 N pinch it along the
y axis, at the top (0, R_o) pointing -y and at the bottom (0, -R_o) pointing +y. Linear static,
quadratic tetrahedra (C3D10) in CalculiX.

Modelled as one eighth: the quarter x >= 0, y >= 0 of the catalog, also cut at mid-width, so
0 <= z <= b/2. Every cut is a symmetry plane of geometry and load, so each carries an exact
symmetry roller (``Support.dofs``):

- x = 0 (the crown section under the load): u_x = 0;
- y = 0 (the horizontal section; it replaces the bottom force): u_y = 0;
- z = 0 (mid-width): u_z = 0.

The three rollers also remove all rigid-body motion. No line support (ADR-006) is needed: the
physical "support" of a pinched ring is the opposite force, which the y = 0 plane stands for.

The load: a narrow land, not a point
------------------------------------
Tessera loads faces, so the top force sits on a small flat face, a land, that is part of the
geometry: the ring is cut by the plane y = sqrt(R_o^2 - c^2), which makes a flat face of
half-width c = h/4 = 1.25 mm (|x| <= c) on top of the outer surface. The cut is at most
c^2 / (2 R_o) = 0.015 mm deep (0.3 % of h) over an arc of 2c = 2.5 mm (0.8 % of the outer
circumference), so it does not change the ring's stiffness measurably. In the eighth model the
land is c x b/2 and carries ``ForceLoad`` (0, -P/4, 0), spread as a uniform traction of
P / (2 c b) = 8 MPa.

Why spreading the load over the land does not change the measured quantity:

- Beam theory. With the force spread over |x| <= c instead of at x = 0, the moment at every
  section beyond the land changes by the constant F c / 2 (F = P/2, the load on a quarter; its
  resultant moves c/2 towards the section). The compatibility condition (crown rotation zero,
  int M dphi fixed, see below) takes that constant up in the crown moment M0, so M beyond the
  land is unchanged; only the moment inside the land (|phi| < c/R) differs, by O(F c). The
  shortening changes by O((c/R)^2) = 6e-4 relative, below 0.1 %.
- Continuum. A load on a short land is a local, self-equilibrated perturbation of the same
  resultant (Saint-Venant); its stresses decay over about one section depth h. The measured
  point is the centroid of the crown section, h/2 below the land, so it does pick up the local
  compression of the material between the land and the centroid. That is of order
  P / (pi E b) * ln(h / c) ~ 4e-5 mm per crown, about 0.05 % of the shortening (a
  Flamant-type estimate, inferred and not derived). The force is a traction on a finite face,
  so no displacement is singular, and the measured value converges (table below).
- Check on the through-thickness picture: the case also prints the inner-surface node of the
  crown section. Plane-section bending alone predicts it shortens less than the centroid by
  the Poisson strain of the crown moment, 2 nu M0 h^2 / (8 E I) = 2.7e-4 mm = 0.16 %; measured
  0.155 % at every mesh, so the section behaves as beam theory says right under the land.

Measured quantity
-----------------
Delta_V: shortening of the loaded diameter between the centroids of the two crown sections,
Delta_V = -2 u_y at (0, R, 0); the y = 0 plane is held in y, so the bottom crown moves by the
same amount the other way. That node exists on every mesh: the mid-width face z = 0 is split
along the centroidal arc r = R (a CAD split line imprinted with gmsh's ``fragment``; the part is
still one solid), so (0, R, 0) is a geometry vertex, the end of that arc on the x = 0 face.
The other end of the arc, (R, 0, 0), gives the increase of the horizontal diameter,
Delta_H = 2 u_x, which is printed as a cross-check and not gated.

Reference: curved bar by Castigliano, with hoop-force, shear and curvature terms
-------------------------------------------------------------------------------
Quarter ring, angle phi from the crown, section centroid at (R sin phi, R cos phi), load
F = P/2 at the crown. By symmetry the crown section has no shear and no rotation, and
horizontal equilibrium with the shear-free y = 0 section gives it no hoop force, so the crown
carries only the unknown moment M0. With M positive when it reduces the curvature (Timoshenko's
convention) and N positive in tension:

    M(phi) = M0 - F R sin(phi),   N(phi) = -F sin(phi),   V(phi) = F cos(phi).

Strain energy of a curved bar with plane sections (Winkler-Bach), per unit length ds = R dphi:

    dU/ds = M^2 / (2 A E e R) + N^2 / (2 A E) - M N / (A E R) + alpha V^2 / (2 A G)

with A = b h, G = E / (2 (1 + nu)), e = R - r_n the shift of the neutral axis
(r_n = h / ln(R_o / R_i) for a rectangle), and alpha = 6/5, the energy shear factor of a
rectangle. For a thin bar A e R -> I = b h^3 / 12 and the coupling term vanishes. Source: S. P.
Timoshenko, *Strength of Materials*, Part I, the chapter on curved bars, deflection of curved
bars by the strain energy (article and equation numbers **unverified**, formula quoted from
memory). Castigliano:

    dU/dM0 = 0 (crown rotation)  ->  M0 = F (R - e) * 2 / pi       (thin: P R / pi)
    Delta_V = dU_ring/dP = 2 dU_quarter/dF

    Delta_V = P R^2 / (A E e) (pi/4 - 2/pi)                         bending (curved bar)
            + P R / (A E) (4/pi - pi/4 - 2 e / (pi R))              hoop force + coupling
            + alpha pi P R / (4 A G)                                 shear

With I* = A e R, a = e/R = I* / (A R^2) and beta = alpha E I* / (G A R^2) this is
Delta_V = P R^3 / (E I*) [K1 pi/4 - 2 K2^2 / pi], K1 = 1 - a + beta, K2 = 1 - a, which is the
form of W. C. Young & R. G. Budynas, *Roark's Formulas for Stress and Strain*, 7th ed., table
of circular-ring formulas, case 1 (table and case **unverified**; form recalled from memory, not
opened). The thin-ring limit (a = beta = 0, I* = I) is the catalog's
(pi/4 - 2/pi) P R^3 / (E I) = 0.14878 P R^3 / (E I). The test derives every step again by
quadrature and finite differences of U, independently of the algebra above.

Horizontal diameter (cross-check, not gated), by the unit-load method with a pair of unit
forces in x at the two end sections (m = -R (1 - cos phi), n = cos phi, v = sin phi):

    Delta_H = P R^2 / (A E e) (2/pi - 1/2) + P R / (A E) (1/2 - 4/pi + 2 e / (pi R))
              + alpha P R / (2 A G)

**Derived here**, not checked against a book; the test re-derives it by Castigliano with a real
horizontal force Q at the side end (crown hoop force Q from equilibrium, M0 eliminated by zero
crown rotation). The thin limit (2/pi - 1/2) P R^3 / (E I) = 0.13662 P R^3 / (E I) is the
classical value. The shear term widens the ring (+0.95 % of Delta_H); with its sign flipped
the measured Delta_H would sit 1.9 % off instead of 0.07 %, so the FE run confirms the sign.

Numbers: I = 52.083 mm^4, e = 0.041694 mm (A e R / I = 1.00067), M0 = 1.5902e3 N mm.

    Delta_V thin ring (catalog formula)   0.170032 mm
    curvature (I -> A e R)                -0.067 %
    hoop force + coupling                 +0.273 %
    shear (alpha = 6/5)                   +1.373 %
    Delta_V reference                     0.172717 mm  (+1.58 % over the thin ring)
    Delta_H reference                     0.156783 mm  (thin 0.156137 mm, +0.41 %)

Delta_V / R = 0.35 %, so linear theory holds. Cowper's kappa (as in TS-01, TS-05) instead of
1/alpha would lower the reference by 0.03 %. Plane sections and sigma_z = 0 are the remaining
beam-theory approximations, O((h/R)^2) and well below the gate.

Measured through Tessera's public API (gmsh mesher, CalculiX 2.21, October 2026; 4 cores):

    size [mm]  el. through h  nodes   Delta_V vs ref.  vs thin ring  Delta_H vs ref.  runtime
    2.5             2          1528      -0.199 %       +1.377 %       -0.092 %       ~1 s
    1.25            4          5967      -0.172 %       +1.404 %       -0.070 %       ~2 s  (gate)
    0.8333          6         18069      -0.169 %       +1.407 %       -0.068 %       ~6 s

The error falls monotonically and has settled by the gate mesh; the -0.17 % left is the
O((h/R)^2) beam-theory approximation (plane sections, sigma_r neglected), not discretisation.
Against the thin-ring formula alone the error would read +1.4 %, inside a 3 % gate: like TS-05,
the gate by itself does not prove the shear and hoop-force terms are present. A 1 % gate would.

Gate: 4 elements through the wall (size h/4); one solve takes about 2 s.
"""

from __future__ import annotations

import itertools
import math
import shutil
from collections.abc import Callable
from pathlib import Path
from typing import Any, NamedTuple

import numpy as np
import pytest

pytestmark = [
    pytest.mark.benchmark,
    pytest.mark.solver,
    pytest.mark.skipif(shutil.which("ccx") is None, reason="CalculiX (ccx) not on PATH"),
]

REFERENCE = "thin_ring"

# Geometry (mm), load (N), material (MPa)
MEAN_RADIUS = 50.0  # R
THICKNESS = 5.0  # h, radial
WIDTH = 5.0  # b, along z (the model is the half 0 <= z <= b/2)
LAND_HALF_WIDTH = THICKNESS / 4  # c: the flat land carrying the load, |x| <= c
FORCE = 100.0  # P, each of the two pinching forces
YOUNGS_MODULUS = 210_000.0  # steel-generic
POISSON_RATIO = 0.30
SHEAR_FACTOR = 6.0 / 5.0  # alpha, energy shear factor of a rectangle

GATE_MESH_SIZE = THICKNESS / 4  # mm: 4 quadratic elements through the ring wall
CONVERGENCE_MESH_SIZES = (THICKNESS / 2, THICKNESS / 4, THICKNESS / 6)  # mm, coarse to fine


class RingResult(NamedTuple):
    """Diameter changes (mm) from one solve."""

    delta_v: float  # shortening of the loaded diameter, centroid of the crown section
    delta_v_inner: float  # the same from the inner-surface node of the crown section
    delta_h: float  # increase of the horizontal diameter, centroid of the side section
    n_nodes: int


def inner_radius() -> float:
    return MEAN_RADIUS - THICKNESS / 2


def outer_radius() -> float:
    return MEAN_RADIUS + THICKNESS / 2


def section_area() -> float:
    """A = b h (mm^2)."""
    return WIDTH * THICKNESS


def second_moment() -> float:
    """I = b h^3 / 12 (mm^4)."""
    return WIDTH * THICKNESS**3 / 12.0


def shear_modulus() -> float:
    return YOUNGS_MODULUS / (2.0 * (1.0 + POISSON_RATIO))


def neutral_axis_shift() -> float:
    """e = R - r_n, r_n = h / ln(R_o / R_i): the Winkler-Bach neutral axis of a rectangle."""
    return MEAN_RADIUS - THICKNESS / math.log(outer_radius() / inner_radius())


def thin_ring_shortening() -> float:
    """Catalog thin-ring value (pi/4 - 2/pi) P R^3 / (E I) (mm)."""
    return (math.pi / 4 - 2 / math.pi) * FORCE * MEAN_RADIUS**3 / (YOUNGS_MODULUS * second_moment())


def shortening_terms(shear_factor: float = SHEAR_FACTOR) -> tuple[float, float, float]:
    """(bending, hoop force + coupling, shear) parts of Delta_V (mm), curved bar."""
    p, r, e = FORCE, MEAN_RADIUS, neutral_axis_shift()
    ae = section_area() * YOUNGS_MODULUS
    bending = p * r**2 / (ae * e) * (math.pi / 4 - 2 / math.pi)
    hoop = p * r / ae * (4 / math.pi - math.pi / 4 - 2 * e / (math.pi * r))
    shear = shear_factor * math.pi * p * r / (4 * section_area() * shear_modulus())
    return bending, hoop, shear


def diameter_shortening() -> float:
    """Reference Delta_V (mm): curved bar with hoop-force, coupling and shear energy."""
    return sum(shortening_terms())


def diameter_increase() -> float:
    """Cross-check Delta_H (mm), same energy, unit-load method (derived in this module)."""
    p, r, e = FORCE, MEAN_RADIUS, neutral_axis_shift()
    ae = section_area() * YOUNGS_MODULUS
    return (
        p * r**2 / (ae * e) * (2 / math.pi - 0.5)
        + p * r / ae * (0.5 - 4 / math.pi + 2 * e / (math.pi * r))
        + SHEAR_FACTOR * p * r / (2 * section_area() * shear_modulus())
    )


def write_ring_step(path: Path) -> Path:
    """Eighth ring x, y >= 0, R_i <= r <= R_o, 0 <= z <= b/2, with the land, as STEP.

    The land is the flat face y = sqrt(R_o^2 - c^2), |x| <= c. The face z = 0 is split along the
    centroidal arc r = R (a CAD split line), so (0, R, 0) and (R, 0, 0) are geometry vertices.
    """
    import gmsh

    owns_session = not gmsh.isInitialized()
    if owns_session:
        gmsh.initialize(interruptible=False)
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add(f"benchmark-{path.stem}")
    try:
        occ = gmsh.model.occ
        half_width = WIDTH / 2
        r_o, r_i = outer_radius(), inner_radius()
        outer = occ.addCylinder(0.0, 0.0, 0.0, 0.0, 0.0, half_width, r_o, angle=math.pi / 2)
        hole = occ.addCylinder(0.0, 0.0, -1.0, 0.0, 0.0, half_width + 2.0, r_i)
        (ring,), _ = occ.cut([(3, outer)], [(3, hole)])
        land_y = math.sqrt(r_o**2 - LAND_HALF_WIDTH**2)
        cap = occ.addBox(-1.0, land_y, -1.0, r_o + 2.0, 2.0 * THICKNESS, half_width + 2.0)
        (ring,), _ = occ.cut([ring], [(3, cap)])
        # Inner half of the z = 0 face: the quarter annulus R_i <= r <= R.
        disk = occ.addDisk(0.0, 0.0, 0.0, MEAN_RADIUS, MEAN_RADIUS)
        inner_disk = occ.addDisk(0.0, 0.0, 0.0, r_i, r_i)
        quadrant = occ.addRectangle(0.0, 0.0, 0.0, r_o, r_o)
        (annulus,), _ = occ.cut([(2, disk)], [(2, inner_disk)])
        (split,), _ = occ.intersect([annulus], [(2, quadrant)])
        occ.fragment([ring], [split])  # imprint: still one solid, z = 0 face split at r = R
        occ.synchronize()
        assert len(gmsh.model.getEntities(3)) == 1
        gmsh.write(str(path))
    finally:
        gmsh.model.remove()
        if owns_session:
            gmsh.finalize()
    return path


def _plane_faces(faces: list[Any], normal: tuple[float, float, float]) -> list[Any]:
    return [f for f in faces if f.normal is not None and np.allclose(f.normal, normal, atol=1e-6)]


def _solve_ring(step_path: Path, mesh_size: float, workdir: Path) -> RingResult:
    """Run the eighth ring through the public API."""
    from tessera.core import describe_faces, edit, new_project, solve

    p = new_project("thin-ring")
    edit.add_part(p, step_path, part_id="ring", material_id="steel-generic")
    faces = describe_faces(p, "ring")
    assert len(faces) == 8, [(f.index, f.kind, f.normal) for f in faces]

    (land,) = _plane_faces(faces, (0.0, 1.0, 0.0))
    assert land.area == pytest.approx(LAND_HALF_WIDTH * WIDTH / 2)
    (crown_cut,) = _plane_faces(faces, (-1.0, 0.0, 0.0))  # x = 0
    (side_cut,) = _plane_faces(faces, (0.0, -1.0, 0.0))  # y = 0
    mid_width = _plane_faces(faces, (0.0, 0.0, -1.0))  # z = 0, split at r = R
    assert len(mid_width) == 2
    assert sum(f.area for f in mid_width) == pytest.approx(
        math.pi / 4 * (outer_radius() ** 2 - inner_radius() ** 2), rel=1e-3
    )  # the land removes 0.02 % of it

    edit.add_support(p, "ring", [crown_cut.index], dofs=(1,), support_id="sym-x")
    edit.add_support(p, "ring", [side_cut.index], dofs=(2,), support_id="sym-y")
    edit.add_support(p, "ring", [f.index for f in mid_width], dofs=(3,), support_id="sym-z")
    edit.add_force_load(p, "ring", [land.index], (0.0, -FORCE / 4, 0.0), load_id="pinch")
    edit.set_mesh_size(p, mesh_size)

    workdir.mkdir(parents=True, exist_ok=True)
    result = solve(p, workdir=workdir)

    nodes = result.mesh.nodes
    u = result.solution.displacements
    tol = 1e-6 * MEAN_RADIUS
    for axis in range(3):
        on_cut = np.abs(nodes[:, axis]) < tol
        assert on_cut.sum() >= 3
        assert np.abs(u[on_cut, axis]).max() < 1e-12, f"cut plane {'xyz'[axis]} = 0 moves normal"

    def node(x: float, y: float) -> int:
        at = (
            (np.abs(nodes[:, 0] - x) < tol)
            & (np.abs(nodes[:, 1] - y) < tol)
            & (np.abs(nodes[:, 2]) < tol)
        )
        assert at.sum() == 1, f"no mesh node at ({x:g}, {y:g}, 0)"
        return int(np.flatnonzero(at)[0])

    delta_v = -2.0 * float(u[node(0.0, MEAN_RADIUS), 1])
    delta_h = 2.0 * float(u[node(MEAN_RADIUS, 0.0), 0])
    assert delta_v > 0.0, "ring did not shorten along the load: force sign or face"
    assert delta_h > 0.0, "ring did not widen across the load"
    return RingResult(
        delta_v=delta_v,
        delta_v_inner=-2.0 * float(u[node(0.0, inner_radius()), 1]),
        delta_h=delta_h,
        n_nodes=result.mesh.n_nodes,
    )


@pytest.fixture(scope="module")
def ring_step(tmp_path_factory: pytest.TempPathFactory) -> Path:
    directory = tmp_path_factory.mktemp("thin-ring-geometry")
    return write_ring_step(directory / "thin_ring_eighth.step")


@pytest.fixture(scope="module")
def ring_result(
    ring_step: Path, tmp_path_factory: pytest.TempPathFactory
) -> Callable[[float], RingResult]:
    """Cached ``mesh_size -> RingResult`` shared by the gate and the sweep."""
    cache: dict[float, RingResult] = {}

    def run(mesh_size: float) -> RingResult:
        if mesh_size not in cache:
            workdir = tmp_path_factory.mktemp(f"thin-ring-h{mesh_size:.4g}")
            cache[mesh_size] = _solve_ring(ring_step, mesh_size, workdir)
        return cache[mesh_size]

    return run


def _quarter_energy(m0: float, f: float, q: float = 0.0, n_points: int = 2001) -> float:
    """Strain energy (N mm) of the quarter ring, Simpson quadrature over phi.

    Loads: F at the crown in -y and Q at the side end (R, 0) in +x. The crown section carries
    the moment M0 and, by horizontal equilibrium of the quarter, the hoop force Q (its shear is
    zero by symmetry). On the free body from the crown to phi the rest of the ring exerts the
    force (Q, F); with the tangent t = (cos phi, -sin phi) and the outward radial direction
    (sin phi, cos phi) this gives N and V, and the moment of the crown loads about the section
    centroid gives M.
    """
    r, e = MEAN_RADIUS, neutral_axis_shift()
    a, modulus, shear_mod = section_area(), YOUNGS_MODULUS, shear_modulus()
    phi = np.linspace(0.0, math.pi / 2, n_points)
    s, c = np.sin(phi), np.cos(phi)
    moment = m0 - f * r * s - q * r * (1.0 - c)
    normal = q * c - f * s
    shear = q * s + f * c
    density = (
        moment**2 / (2 * a * modulus * e * r)
        + normal**2 / (2 * a * modulus)
        - moment * normal / (a * modulus * r)
        + SHEAR_FACTOR * shear**2 / (2 * a * shear_mod)
    )
    weights = np.ones(n_points)
    weights[1:-1:2], weights[2:-1:2] = 4.0, 2.0
    return float(r * (phi[1] - phi[0]) / 3.0 * np.dot(weights, density))


def _compatible_energy(f: float, q: float) -> float:
    """Quarter-ring energy with M0 from zero crown rotation, dU/dM0 = 0 (U is quadratic)."""
    step = f * MEAN_RADIUS
    u0 = _quarter_energy(0.0, f, q)
    up, um = _quarter_energy(step, f, q), _quarter_energy(-step, f, q)
    slope, curvature = (up - um) / (2 * step), (up - 2 * u0 + um) / step**2
    return _quarter_energy(-slope / curvature, f, q)


def test_closed_form_is_consistent() -> None:
    """Castigliano by quadrature and finite differences reproduces both closed forms."""
    r, f = MEAN_RADIUS, FORCE / 2
    e = neutral_axis_shift()

    # Vertical: M0 = 2 F (R - e) / pi makes the crown rotation dU/dM0 zero.
    m0 = 2.0 * f * (r - e) / math.pi
    dm = 1e-3 * m0
    dudm0 = (_quarter_energy(m0 + dm, f) - _quarter_energy(m0 - dm, f)) / (2 * dm)
    assert abs(dudm0) < 1e-9 * f * r**3 / (YOUNGS_MODULUS * second_moment())
    assert math.isclose(m0, FORCE * r / math.pi * (1 - e / r), rel_tol=1e-12)
    # Delta_V = dU_ring/dP = 2 dU_quarter/dF, with M0 eliminated (the total derivative).
    df = 1e-3 * f
    delta_quarter = (_compatible_energy(f + df, 0.0) - _compatible_energy(f - df, 0.0)) / (2 * df)
    assert math.isclose(2.0 * delta_quarter, diameter_shortening(), rel_tol=1e-8)

    # Roark's form with I* = A e R: K1 pi/4 - 2 K2^2 / pi.
    i_star = section_area() * e * r
    a_r = e / r
    beta = SHEAR_FACTOR * YOUNGS_MODULUS * i_star / (shear_modulus() * section_area() * r**2)
    k1, k2 = 1.0 - a_r + beta, 1.0 - a_r
    roark = FORCE * r**3 / (YOUNGS_MODULUS * i_star) * (k1 * math.pi / 4 - 2 * k2**2 / math.pi)
    assert math.isclose(roark, diameter_shortening(), rel_tol=1e-12)

    # Thin limit: the catalog coefficient, and bending = thin value with I -> A e R.
    assert math.isclose(math.pi / 4 - 2 / math.pi, 0.14878, abs_tol=5e-6)
    bending, hoop, shear = shortening_terms()
    assert math.isclose(bending * i_star / second_moment(), thin_ring_shortening(), rel_tol=1e-12)
    assert math.isclose(i_star / second_moment(), 1.0, rel_tol=1e-3)
    assert 0.0 < hoop / thin_ring_shortening() < 0.005
    assert 0.01 < shear / thin_ring_shortening() < 0.02
    assert MEAN_RADIUS / THICKNESS >= 10.0  # catalog: R/h >= 10

    # Horizontal: Delta_H / 2 = dU/dQ at Q = 0, with M0 eliminated for every Q.
    dq = 1e-3 * f
    delta_h_half = (_compatible_energy(f, dq) - _compatible_energy(f, -dq)) / (2 * dq)
    assert math.isclose(2.0 * delta_h_half, diameter_increase(), rel_tol=1e-8)
    thin_h = (2 / math.pi - 0.5) * FORCE * r**3 / (YOUNGS_MODULUS * second_moment())
    assert math.isclose(diameter_increase(), thin_h, rel_tol=0.01)


def test_reference_matches_formula(load_reference: Callable[[str], dict[str, Any]]) -> None:
    """The stored reference can never drift from the closed-form value."""
    ref = load_reference(REFERENCE)
    assert math.isclose(
        ref["reference_diameter_shortening_mm"], diameter_shortening(), rel_tol=1e-9
    )
    assert 0.0 < ref["rel_tolerance"] <= 0.03  # CATALOG target: within 3 %


def test_thin_ring_diameter_shortening(
    ring_result: Callable[[float], RingResult],
    load_reference: Callable[[str], dict[str, Any]],
) -> None:
    # Solve first: fail (not skip) if the solve breaks, before the reference check.
    delta_v = ring_result(GATE_MESH_SIZE).delta_v

    ref = load_reference(REFERENCE)
    assert delta_v == pytest.approx(
        ref["reference_diameter_shortening_mm"], rel=ref["rel_tolerance"]
    )


def test_thin_ring_mesh_convergence(
    ring_result: Callable[[float], RingResult],
    load_reference: Callable[[str], dict[str, Any]],
) -> None:
    # Solve first: fail (not skip) if the solve breaks, before the reference check.
    results = [ring_result(size) for size in CONVERGENCE_MESH_SIZES]

    ref = load_reference(REFERENCE)
    expected = ref["reference_diameter_shortening_mm"]
    tolerance = ref["rel_tolerance"]
    expected_h = diameter_increase()

    print(
        f"\nThin ring R/h = {MEAN_RADIUS / THICKNESS:g}, P = {FORCE:g} N: Delta_V ref "
        f"{expected:.6e} mm (thin ring {thin_ring_shortening():.6e}); Delta_H cross-check "
        f"{expected_h:.6e} mm"
    )
    print(
        f"{'h [mm]':>8} {'nodes':>8} {'Delta_V [mm]':>13} {'rel. error':>11}"
        f" {'vs thin':>9} {'inner':>9} {'Delta_H err':>12}"
    )
    for size, res in zip(CONVERGENCE_MESH_SIZES, results, strict=True):
        print(
            f"{size:>8.4g} {res.n_nodes:>8d} {res.delta_v:>13.6e}"
            f" {res.delta_v / expected - 1.0:>+11.4%}"
            f" {res.delta_v / thin_ring_shortening() - 1.0:>+9.3%}"
            f" {res.delta_v_inner / res.delta_v - 1.0:>+9.3%}"
            f" {res.delta_h / expected_h - 1.0:>+12.4%}"
        )
    print("(inner: crown inner-surface node vs the centroid node)")

    node_counts = [res.n_nodes for res in results]
    assert node_counts == sorted(node_counts) and len(set(node_counts)) == len(node_counts), (
        "mesh_size is not refining the mesh"
    )
    errors = [abs(res.delta_v - expected) / expected for res in results]
    for coarse, fine in itertools.pairwise(errors):
        assert fine <= coarse or fine < tolerance / 10
    assert errors[-1] < tolerance
