"""Simply supported circular plate under uniform pressure, centre deflection (TS-08 variant B).

Problem (mm, N, MPa per ADR-002)
--------------------------------
The variant A disc (``test_clamped_plate.py``): steel, radius a = 20, thickness h = 2
(a/h = 10), E = 210 000 MPa, nu = 0.30 (``steel-generic``), occupying 0 <= z <= h; uniform
pressure q = 1 MPa on the top face z = h (``PressureLoad``, positive = into the part, i.e. -z).
The plate rests on a ring: a line support (ADR-006) holds u_z = 0 on the bottom outer edge,
the edge between the bottom face z = 0 and the cylindrical face r = a, named by that face pair
in ``edit.add_edge_support(..., dofs=(3,))``. Only u_z is held, so the ring is free to move
radially and the plate is free to rotate about it. Linear static, quadratic tetrahedra (C3D10)
in CalculiX.

Modelled as a quarter disc (x, y >= 0) with symmetry rollers u_x = 0 on the x = 0 cut plane and
u_y = 0 on y = 0. Both rollers are exact for the axisymmetric problem. Together with the ring
they remove all rigid-body motion: the ring holds z and the two tilts, the rollers hold x, y
and the rotation about z.

Measured quantity
-----------------
w0: the centre deflection of the mid-surface, -u_z of the mesh node at x = y = 0, z = h/2,
which exists on every mesh (see variant A: the straight edge x = y = 0 is split uniformly and
quadratic elements put a node at its midpoint). The top and bottom centre nodes deflect less
than the mid-surface because of the Poisson through-thickness strain of bending (variant A
explains it). Here the predicted offset is nu / (1 - nu) * (k + 1) / (2 k) * (h/a)^2 of w_K
with k = (5 + nu) / (1 + nu), i.e. 0.264 % of w0. Measured: top 0.257 %, bottom 0.265 %. Only
the deflection is gated: the line support carries a finite force per length, so the stresses
at the ring are singular (ADR-006, ADR-007).

Reference: Kirchhoff bending + Mindlin transverse shear
-------------------------------------------------------
    w0 = w_K + w_S
    w_K = (5 + nu) / (1 + nu) * q a^4 / (64 D),   D = E h^3 / (12 (1 - nu^2))
    w_S = q a^2 / (4 kappa G h),                  G = E / (2 (1 + nu)), kappa = 5/6

    w_S / w_K = 16 (1 + nu) / (5 (5 + nu) (1 - nu)) (h/a)^2   (kappa = 5/6)

- w_K: thin-plate (Kirchhoff) solution of the simply supported circular plate under uniform
  load, w_K(r) = q (a^2 - r^2) / (64 D) * ((5 + nu) / (1 + nu) a^2 - r^2).
  S. P. Timoshenko & S. Woinowsky-Krieger, *Theory of Plates and Shells*, 2nd ed., McGraw-Hill
  1959, chapter 3, symmetrical bending of circular plates (section number **unverified**;
  centre value (5 + nu) q a^4 / (64 (1 + nu) D) **verified** only via public web-page
  excerpts, see CATALOG.md TS-08). The test re-checks it: D lap(lap(w)) = q, w(a) = 0 and
  M_r(a) = -D (w'' + nu w'/r)(a) = 0.
- w_S: **derived here**, not read in a source (status **unverified** against a book), with the
  same argument as variant A; only the theories are cited (R. D. Mindlin, J. Appl. Mech. 18
  (1951) 31-38; E. Reissner, J. Appl. Mech. 12 (1945) A69-A77). Why it carries over to the
  simple support: (1) Q_r = -q r / 2 follows from vertical equilibrium of the central disc of
  radius r alone, 2 pi r Q_r = -q pi r^2, whatever holds the edge (statically determinate).
  (2) The Mindlin rotation psi obeys the moment equilibrium d(r M_r)/dr - M_theta = r Q_r with
  M_r, M_theta given by psi exactly as Kirchhoff's are by the slope, the same Q_r, and the same
  boundary data psi(0) = 0 and M_r(a) = 0 (a free rotation, as for Kirchhoff). That ODE has a
  unique solution, so psi is the Kirchhoff slope and the bending part of w is w_K. In the
  axisymmetric case there is no twisting moment, so the "hard" and "soft" Mindlin simple
  supports coincide and there is no edge boundary layer. (3) The shear strain adds
  dw_S/dr = Q_r / (kappa G h) = -q r / (2 kappa G h), and w(a) = 0 gives
  w_S(r) = q (a^2 - r^2) / (4 kappa G h), the variant A term. The test checks (2) on the
  Kirchhoff moments and (3) by differentiation.

Numbers: D = 153 846.15 N mm, w_K = 0.066250000 mm, w_S = 0.000742857 mm (1.121 % of w_K),
w0 = 0.066992857 mm. w0 / h = 3.3 %, so linear theory holds. The shear term is relatively
smaller than in variant A (4.57 %) because the simply supported plate bends about four times as
much. A Kirchhoff-only reference would differ by 1.1 %, inside the gate, so unlike variant A
this case does not on its own prove transverse shear is present.

The 3D solid versus the plate reference: the support sits on the bottom edge
--------------------------------------------------------------------------------
Plate theory holds the mid-surface at r = a. The solid is held at the bottom corner, h/2 below
it. Effects of that, from analysis and measured:

- No edge moment and no membrane force. The reaction is vertical, so it has no moment about the
  mid-surface point above it. Since only u_z is held, the ring moves out freely. Measured
  radially at the gate mesh: u_r = 5.12e-3 mm against psi(a) h/2 = 5.07e-3 mm (Kirchhoff slope
  at the edge times h/2), so the support works as a roller, not a pin.
- Rotation about the corner instead of the mid-surface moves the mid-surface edge sideways by
  psi h/2. Its vertical part is O(psi^2 h), which linear analysis drops.
- The mid-surface at r = a deflects below the held bottom corner. Part of this is the Poisson
  through-thickness strain: M_r(a) = 0, but M_theta(a) is not zero, giving
  nu / (1 - nu) * h^2 / 8 * (-lap w_K)(a) = nu / (1 - nu) * (k - 3) / (2 k) (h/a)^2 w_K =
  0.057 % of w0. The larger part is the local compliance of the line load at a 90 degree
  corner. It grows like ln(1 / element size) and does not converge: a line force on an elastic
  wedge has a logarithmically singular displacement. Measured mid-surface deflection at r = a,
  in % of w0: 0.23 / 0.32 / 0.37 / 0.41 / 0.45 at element size 2 / 1 / 0.667 / 0.5 / 0.4,
  i.e. +0.13 to +0.18 % of w0 per unit of ln(1 / size). The ring reaction is q a / 2 = 10 N/mm.
  Refining from 0.4 to 0.1 mm would add only about another 0.2 %, so the effect matters little
  at any mesh a user can afford. It does mean the gate value creeps upward slowly under
  refinement instead of converging.
- Holding the ring in all directions (a pinned line, dofs (1, 2, 3)) instead stops the bottom
  fibres from stretching. Membrane action then halves the deflection (measured -52 % at the
  gate mesh). The case asserts that the ring moves outward, so a support that silently holds
  more than u_z fails here even before the tolerance check.

The same model with the ring on the mid-surface (the rim split at z = h/2 and the support on
that circle, measured once in a scratch run and not part of this case) gives -0.39 / -0.22 /
-0.20 % at 2 / 1 / 0.667 mm. That is the same O((h/a)^2) offset from the Mindlin reference as
variant A (-0.3 %; sigma_z neglected, kappa, the 3D edge). The bottom-edge support adds +0.41 %
at the gate mesh, almost all of it the mid-surface rim offset above.

Measured through Tessera's public API (gmsh mesher, CalculiX 2.21, October 2026):

    size [mm]  el. through h  nodes   w0 mid-surface   rim mid-surface   runtime
    2.0             1          1597      -0.08 %           +0.23 %       ~0.7 s
    1.0             2          6771      +0.14 %           +0.32 %       ~3 s
    0.6667          3         20193      +0.22 %           +0.37 %       ~12 s  (gate)
    0.5             4         ~40 k      +0.25 %           +0.41 %       ~31 s  (not in the sweep)
    0.4             5         ~75 k      +0.29 %           +0.45 %       ~92 s  (not in the sweep)

The coarsest mesh is below the reference (bending too stiff, as in variant A), and refinement
moves w0 up. Between 2 and 1 mm the shift is +0.22 %, and between 1 and 0.667 mm it is
+0.08 %. Past about 0.667 mm the bending error has largely converged, and what is left is the
slow logarithmic rise from the ring above. The mesh-convergence test therefore checks
that the steps shrink and every size is inside the gate, not that the error falls
monotonically.

Gate: 3 element layers (size h/3), as in variant A; one solve takes about 12 s.
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

REFERENCE = "simply_supported_plate"

# Geometry (mm), load (MPa), material (MPa): the variant A disc.
RADIUS = 20.0  # a
THICKNESS = 2.0  # h
PRESSURE = 1.0  # q, on the top face, acting in -z
YOUNGS_MODULUS = 210_000.0  # steel-generic
POISSON_RATIO = 0.30
SHEAR_CORRECTION = 5.0 / 6.0  # kappa (Mindlin/Reissner)

GATE_MESH_SIZE = THICKNESS / 3  # mm: 3 quadratic elements through the thickness
CONVERGENCE_MESH_SIZES = (THICKNESS, THICKNESS / 2, THICKNESS / 3)  # mm, coarse to fine


class PlateResult(NamedTuple):
    """Deflections (mm, positive in the load direction -z) from one solve."""

    centre_mid: float  # w0: mid-surface node at x = y = 0, z = h/2
    centre_top: float
    centre_bottom: float
    rim_mid: float  # mean over the mid-surface nodes at r = a (the ring itself is held at z = 0)
    ring_radial: float  # mean radial displacement of the held ring nodes
    n_nodes: int


def flexural_rigidity(
    thickness: float = THICKNESS,
    youngs_modulus: float = YOUNGS_MODULUS,
    poisson_ratio: float = POISSON_RATIO,
) -> float:
    """D (N mm) = E h^3 / (12 (1 - nu^2))."""
    return youngs_modulus * thickness**3 / (12.0 * (1.0 - poisson_ratio**2))


def kirchhoff_deflection(
    r: float,
    radius: float = RADIUS,
    pressure: float = PRESSURE,
    thickness: float = THICKNESS,
    poisson_ratio: float = POISSON_RATIO,
) -> float:
    """Thin-plate deflection (mm) of the simply supported plate.

    w_K(r) = q (a^2 - r^2) / (64 D) * ((5 + nu) / (1 + nu) a^2 - r^2).
    """
    k = (5.0 + poisson_ratio) / (1.0 + poisson_ratio)
    d = flexural_rigidity(thickness, poisson_ratio=poisson_ratio)
    return pressure * (radius**2 - r**2) * (k * radius**2 - r**2) / (64.0 * d)


def shear_deflection(
    r: float,
    radius: float = RADIUS,
    pressure: float = PRESSURE,
    thickness: float = THICKNESS,
    youngs_modulus: float = YOUNGS_MODULUS,
    poisson_ratio: float = POISSON_RATIO,
    kappa: float = SHEAR_CORRECTION,
) -> float:
    """Mindlin transverse-shear deflection (mm) w_S(r) = q (a^2 - r^2) / (4 kappa G h)."""
    shear_modulus = youngs_modulus / (2.0 * (1.0 + poisson_ratio))
    return pressure * (radius**2 - r**2) / (4.0 * kappa * shear_modulus * thickness)


def centre_deflection() -> float:
    """Reference w0 (mm) = w_K(0) + w_S(0)."""
    return kirchhoff_deflection(0.0) + shear_deflection(0.0)


def write_quarter_disc_step(path: Path) -> Path:
    """Quarter disc x, y >= 0, r <= a, 0 <= z <= h as STEP (gmsh OpenCASCADE).

    The same geometry as variant A. It is duplicated rather than imported because the case
    modules are not a package (pytest imports them with ``--import-mode=importlib``).
    """
    import gmsh

    owns_session = not gmsh.isInitialized()
    if owns_session:
        gmsh.initialize(interruptible=False)
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add(f"benchmark-{path.stem}")
    try:
        gmsh.model.occ.addCylinder(0.0, 0.0, 0.0, 0.0, 0.0, THICKNESS, RADIUS, angle=math.pi / 2)
        gmsh.model.occ.synchronize()
        gmsh.write(str(path))
    finally:
        gmsh.model.remove()
        if owns_session:
            gmsh.finalize()
    return path


def _face_with_normal(faces: list[Any], normal: tuple[float, float, float]) -> int:
    matches = [
        f.index for f in faces if f.normal is not None and np.allclose(f.normal, normal, atol=1e-6)
    ]
    assert len(matches) == 1, f"expected one face with normal {normal}, got {matches}"
    return matches[0]


def _curved_face_with_area(faces: list[Any], area: float) -> int:
    matches = [
        f.index for f in faces if f.normal is None and math.isclose(f.area, area, rel_tol=1e-6)
    ]
    assert len(matches) == 1, f"expected one curved face of area {area:.4f}, got {matches}"
    return matches[0]


def _solve_simply_supported_plate(step_path: Path, mesh_size: float, workdir: Path) -> PlateResult:
    """Run the quarter disc on a ring support through the public API."""
    from tessera.core import describe_edges, describe_faces, edit, new_project, solve
    from tessera.model import PartInstance

    p = new_project("simply-supported-plate")
    p.assembly.parts.append(
        PartInstance(id="disc", name="disc", geometry=step_path, material_id="steel-generic")
    )
    faces = describe_faces(p, "disc")
    top = _face_with_normal(faces, (0.0, 0.0, 1.0))
    bottom = _face_with_normal(faces, (0.0, 0.0, -1.0))
    sym_x = _face_with_normal(faces, (-1.0, 0.0, 0.0))  # x = 0 cut plane
    sym_y = _face_with_normal(faces, (0.0, -1.0, 0.0))  # y = 0 cut plane
    rim = _curved_face_with_area(faces, math.pi / 2 * RADIUS * THICKNESS)  # 2 pi a h / 4
    assert next(f for f in faces if f.index == rim).kind.lower() == "cylinder"
    assert next(f for f in faces if f.index == top).centroid[2] == pytest.approx(THICKNESS)

    ring = (min(bottom, rim), max(bottom, rim))
    (ring_info,) = [e for e in describe_edges(p, "disc") if e.faces == ring]
    assert ring_info.kind.lower() == "circle"
    assert ring_info.length == pytest.approx(math.pi / 2 * RADIUS)  # quarter of 2 pi a
    assert ring_info.point[2] == pytest.approx(0.0, abs=1e-9)  # the bottom edge, not the top

    edit.add_edge_support(p, "disc", [ring], dofs=(3,), support_id="ring")
    edit.add_support(p, "disc", [sym_x], dofs=(1,), support_id="sym-x")
    edit.add_support(p, "disc", [sym_y], dofs=(2,), support_id="sym-y")
    edit.add_pressure_load(p, "disc", [top], PRESSURE, load_id="q")
    edit.set_mesh_size(p, mesh_size)

    workdir.mkdir(parents=True, exist_ok=True)
    result = solve(p, workdir=workdir)

    nodes = result.mesh.nodes
    u = result.solution.displacements
    tol = 1e-6 * RADIUS
    r = np.hypot(nodes[:, 0], nodes[:, 1])
    on_rim = np.abs(r - RADIUS) < tol
    on_ring = on_rim & (np.abs(nodes[:, 2]) < tol)
    assert on_ring.sum() >= 3, "no mesh nodes on the bottom outer edge"
    assert np.abs(u[on_ring, 2]).max() < 1e-12, "ring support does not hold u_z"
    radial = (u[on_ring, 0] * nodes[on_ring, 0] + u[on_ring, 1] * nodes[on_ring, 1]) / RADIUS
    # Only u_z is held: the bottom fibres stretch, so the ring must move out (roller, not pin).
    assert radial.min() > 0.0, "ring support also holds the radial direction"
    for axis in (0, 1):
        on_cut = np.abs(nodes[:, axis]) < tol
        assert np.abs(u[on_cut, axis]).max() < 1e-12, f"cut plane {'xy'[axis]} = 0 moves normal"

    on_axis = (np.abs(nodes[:, 0]) < tol) & (np.abs(nodes[:, 1]) < tol)
    z = nodes[on_axis, 2]
    w = -u[on_axis, 2]

    def at_height(height: float) -> float:
        here = np.abs(z - height) < tol
        assert here.sum() == 1, f"no mesh node at the centre x = y = 0, z = {height:g}"
        return float(w[here][0])

    w_mid = at_height(THICKNESS / 2)
    assert w_mid > 0.0, "plate did not deflect in the load direction: pressure sign or face"
    rim_mid = on_rim & (np.abs(nodes[:, 2] - THICKNESS / 2) < tol)
    assert rim_mid.any(), "no mesh node on the mid-surface at r = a"
    return PlateResult(
        centre_mid=w_mid,
        centre_top=at_height(THICKNESS),
        centre_bottom=at_height(0.0),
        rim_mid=float(-u[rim_mid, 2].mean()),
        ring_radial=float(radial.mean()),
        n_nodes=result.mesh.n_nodes,
    )


@pytest.fixture(scope="module")
def disc_step(tmp_path_factory: pytest.TempPathFactory) -> Path:
    directory = tmp_path_factory.mktemp("simply-supported-plate-geometry")
    return write_quarter_disc_step(directory / "simply_supported_quarter_disc.step")


@pytest.fixture(scope="module")
def simply_supported_plate(
    disc_step: Path, tmp_path_factory: pytest.TempPathFactory
) -> Callable[[float], PlateResult]:
    """Cached ``mesh_size -> PlateResult`` shared by the gate and the sweep."""
    cache: dict[float, PlateResult] = {}

    def run(mesh_size: float) -> PlateResult:
        if mesh_size not in cache:
            workdir = tmp_path_factory.mktemp(f"simply-supported-plate-h{mesh_size:.4g}")
            cache[mesh_size] = _solve_simply_supported_plate(disc_step, mesh_size, workdir)
        return cache[mesh_size]

    return run


def test_closed_form_is_consistent() -> None:
    """Kirchhoff plate equation and simple-support conditions; Mindlin shear term."""
    a, q, nu = RADIUS, PRESSURE, POISSON_RATIO
    d = flexural_rigidity()
    k = (5.0 + nu) / (1.0 + nu)
    # w_K(r) = q / (64 D) (a^2 - r^2)(k a^2 - r^2) as a polynomial in r.
    w_k = np.polynomial.Polynomial([k * a**4, 0.0, -(k + 1.0) * a**2, 0.0, 1.0]) * (q / (64.0 * d))
    for r in (0.0, 0.3 * a, 0.7 * a, a):
        assert math.isclose(w_k(r), kirchhoff_deflection(r), rel_tol=1e-12, abs_tol=1e-15)

    def over_r(f: np.polynomial.Polynomial) -> np.polynomial.Polynomial:
        assert f.coef[0] == 0.0  # odd in r: f/r is again a polynomial
        return np.polynomial.Polynomial(f.coef[1:])

    def laplacian(f: np.polynomial.Polynomial) -> np.polynomial.Polynomial:
        """Axisymmetric Laplacian f'' + f'/r."""
        return f.deriv(2) + over_r(f.deriv())

    biharmonic = laplacian(laplacian(w_k))
    for r in (0.0, 0.5 * a, a):
        assert math.isclose(d * biharmonic(r), q, rel_tol=1e-12)
    assert math.isclose(w_k(a), 0.0, abs_tol=1e-15)

    # Kirchhoff moments (w positive in the load direction).
    m_r = -d * (w_k.deriv(2) + nu * over_r(w_k.deriv()))
    m_theta = -d * (over_r(w_k.deriv()) + nu * w_k.deriv(2))
    assert math.isclose(m_r(a), 0.0, abs_tol=1e-9 * q * a**2)  # simple support: no edge moment
    assert abs(m_theta(a)) > 0.1 * q * a**2 / 16  # the hoop moment at the edge is not zero
    # Moment equilibrium d(r M_r)/dr - M_theta = r Q_r with the statically determinate
    # Q_r = -q r / 2. Mindlin's rotation obeys this same equation with the same data
    # (psi(0) = 0, M_r(a) = 0), so it equals the Kirchhoff slope.
    r_poly = np.polynomial.Polynomial([0.0, 1.0])
    residual = (r_poly * m_r).deriv() - m_theta + r_poly * r_poly * (q / 2.0)
    assert np.allclose(residual.coef, 0.0, atol=1e-9 * q * a**2)
    # The same Q_r from the plate equation, Q_r = -D d(lap w)/dr.
    for r in (0.25 * a, 0.75 * a):
        assert math.isclose(-d * laplacian(w_k).deriv()(r), -q * r / 2.0, rel_tol=1e-12)

    # Shear part: kappa G h dw_S/dr = Q_r = -q r / 2, w_S(a) = 0.
    kappa_g_h = SHEAR_CORRECTION * YOUNGS_MODULUS / (2.0 * (1.0 + nu)) * THICKNESS
    step = 1e-4 * a
    for r in (0.25 * a, 0.5 * a, 0.75 * a):
        slope = (shear_deflection(r + step) - shear_deflection(r - step)) / (2.0 * step)
        assert math.isclose(kappa_g_h * slope, -q * r / 2.0, rel_tol=1e-6)
    assert math.isclose(shear_deflection(a), 0.0, abs_tol=1e-15)

    assert math.isclose(kirchhoff_deflection(0.0), k * q * a**4 / (64.0 * d), rel_tol=1e-12)
    ratio = shear_deflection(0.0) / kirchhoff_deflection(0.0)
    expected = 16.0 * (1.0 + nu) / (5.0 * (5.0 + nu) * (1.0 - nu)) * (THICKNESS / RADIUS) ** 2
    assert math.isclose(ratio, expected, rel_tol=1e-12)


def test_reference_matches_formula(load_reference: Callable[[str], dict[str, Any]]) -> None:
    """The stored reference can never drift from the closed-form value."""
    ref = load_reference(REFERENCE)
    assert math.isclose(ref["reference_centre_deflection_mm"], centre_deflection(), rel_tol=1e-9)
    assert 0.0 < ref["rel_tolerance"] <= 0.03  # ADR-006 / CATALOG target: within 3 %


def test_simply_supported_plate_centre_deflection(
    simply_supported_plate: Callable[[float], PlateResult],
    load_reference: Callable[[str], dict[str, Any]],
) -> None:
    # Solve first: fail (not skip) if the solve breaks, before the reference check.
    w0 = simply_supported_plate(GATE_MESH_SIZE).centre_mid

    ref = load_reference(REFERENCE)
    assert w0 == pytest.approx(ref["reference_centre_deflection_mm"], rel=ref["rel_tolerance"])


def test_simply_supported_plate_mesh_convergence(
    simply_supported_plate: Callable[[float], PlateResult],
    load_reference: Callable[[str], dict[str, Any]],
) -> None:
    # Solve first: fail (not skip) if the solve breaks, before the reference check.
    results = [simply_supported_plate(size) for size in CONVERGENCE_MESH_SIZES]

    ref = load_reference(REFERENCE)
    expected = ref["reference_centre_deflection_mm"]
    tolerance = ref["rel_tolerance"]

    print(
        f"\nSimply supported circular plate a/h = {RADIUS / THICKNESS:g}, q = {PRESSURE:g} MPa, "
        f"ring u_z = 0 on the bottom outer edge: w0 ref {expected:.6e} mm (Kirchhoff "
        f"{kirchhoff_deflection(0.0):.6e} + shear {shear_deflection(0.0):.6e})"
    )
    print(
        f"{'h [mm]':>8} {'nodes':>8} {'w0 mid [mm]':>13} {'rel. error':>11} {'top':>8}"
        f" {'bottom':>8} {'rim mid':>8} {'ring u_r':>10}"
    )
    for size, res in zip(CONVERGENCE_MESH_SIZES, results, strict=True):
        w = res.centre_mid
        print(
            f"{size:>8.4g} {res.n_nodes:>8d} {w:>13.6e} {w / expected - 1.0:>+11.4%}"
            f" {res.centre_top / w - 1.0:>+8.3%} {res.centre_bottom / w - 1.0:>+8.3%}"
            f" {res.rim_mid / expected:>+8.3%} {res.ring_radial:>10.3e}"
        )
    print("(top, bottom: vs the mid-surface; rim mid: mid-surface at r = a, % of the reference)")

    node_counts = [res.n_nodes for res in results]
    assert node_counts == sorted(node_counts) and len(set(node_counts)) == len(node_counts), (
        "mesh_size is not refining the mesh"
    )
    w0 = [res.centre_mid for res in results]
    # The line support's local compliance grows like ln(1 / size), so w0 has no strict limit
    # (see the module docstring); require that each refinement changes it less than the one
    # before, and that every mesh sits inside the gate.
    steps = [abs(fine - coarse) for coarse, fine in itertools.pairwise(w0)]
    assert all(later < earlier for earlier, later in itertools.pairwise(steps)), (
        f"refinement steps are not shrinking: {steps}"
    )
    errors = [abs(w - expected) / expected for w in w0]
    assert max(errors) < tolerance
