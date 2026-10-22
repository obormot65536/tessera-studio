"""Clamped circular plate under uniform pressure, centre deflection (TS-08 variant A; gate 2 %).

Problem (mm, N, MPa per ADR-002)
--------------------------------
Steel disc of radius a = 20 and thickness h = 2 (a/h = 10), E = 210 000 MPa, nu = 0.30
(``steel-generic``). The disc occupies 0 <= z <= h; uniform pressure q = 1 MPa acts on the top
face z = h (``PressureLoad``, positive = into the part, i.e. in -z). The outer cylindrical face
r = a is fully clamped (all three components of every node fixed). Linear static, quadratic
tetrahedra (C3D10) in CalculiX.

Modelled as a quarter disc (x, y >= 0) with symmetry rollers u_x = 0 on the x = 0 cut plane and
u_y = 0 on y = 0. The problem is axisymmetric, so both rollers are exact; the clamped edge
removes all rigid-body motion.

Variant B (simply supported: a ring line support u_z = 0 on the bottom outer edge, ADR-006) is
``test_simply_supported_plate.py``; it is not part of this case.

Measured quantity
-----------------
w0: the centre deflection of the mid-surface, -u_z of the mesh node at x = y = 0, z = h/2.
That node exists on every mesh: it lies on the straight geometry edge x = y = 0 (the
intersection of the two cut planes), which gmsh splits uniformly, and with quadratic elements
an edge split into n segments always has a node at its midpoint.

Why the mid-surface: plate theory describes the mid-surface. The top and bottom surfaces
deflect less than the mid-surface because of the Poisson through-thickness strain of bending,
eps_z = -nu / (1 - nu) (eps_r + eps_theta) (plane stress, sigma_z neglected), which is odd in
z and integrates over each half thickness to an offset of nu / (1 - nu) (h/a)^2 w_K = 0.43 % at
the centre (computed from the Kirchhoff curvature below). The top/bottom mean has the same
offset (both surfaces move the same way), so it is not a substitute for the mid-surface here.
Measured: 0.40 % (see the table below).

Reference: Kirchhoff bending + Mindlin transverse shear
-------------------------------------------------------
    w0 = w_K + w_S
    w_K = q a^4 / (64 D),       D = E h^3 / (12 (1 - nu^2))
    w_S = q a^2 / (4 kappa G h), G = E / (2 (1 + nu)), kappa = 5/6

    w_S / w_K = 16 D / (kappa G h a^2) = 16 / (5 (1 - nu)) (h/a)^2   (kappa = 5/6)

- w_K: thin-plate (Kirchhoff) solution of the clamped circular plate under uniform load.
  S. P. Timoshenko & S. Woinowsky-Krieger, *Theory of Plates and Shells*, 2nd ed., McGraw-Hill
  1959, chapter 3, symmetrical bending of circular plates (section number **unverified**;
  formula **verified** only via public web-page excerpts, see CATALOG.md TS-08). The test
  re-checks it: w_K(r) = q (a^2 - r^2)^2 / (64 D) satisfies D lap(lap(w)) = q,
  w(a) = w'(a) = 0.
- w_S: **derived here**, not read in a source (status **unverified** against a book). The
  plate theory is Mindlin's (R. D. Mindlin, J. Appl. Mech. 18 (1951) 31-38) or Reissner's
  (E. Reissner, J. Appl. Mech. 12 (1945) A69-A77); only the theories are cited, not this
  formula. Derivation: the axisymmetric shear force follows from vertical equilibrium of the
  central disc of radius r alone, 2 pi r Q_r = -q pi r^2, so Q_r = -q r / 2 (statically
  determinate). In Mindlin theory the normal rotation psi obeys the same moment equations as
  the Kirchhoff slope with the same data (Q_r, psi(0) = 0, clamp psi(a) = 0), so psi equals
  the Kirchhoff slope and the bending part of w is w_K. The shear strain adds
  dw_S/dr = Q_r / (kappa G h) = -q r / (2 kappa G h); with w_S(a) = 0 this integrates to
  w_S(r) = q (a^2 - r^2) / (4 kappa G h). The test checks both steps.

Numbers: D = 153 846.15 N mm, w_K = 0.016250000 mm, w_S = 0.000742857 mm (4.571 % of w_K),
w0 = 0.016992857 mm. w0 / h = 0.85 %, so linear theory holds.

Why a/h = 10 and a shear correction (rather than a/h >= 40)
----------------------------------------------------------
A 3D mesh with at least 3 quadratic elements through the thickness needs element size h/3
everywhere (Tessera meshes with one global size), so the node count grows as (a/h)^2. Measured
on this machine at size h/3: a/h = 10, 21 k nodes, ~9 s; a/h = 15, 45 k nodes, ~33 s;
a/h = 20, 78 k nodes, ~77 s. a/h = 40 would be far beyond the CI budget. At a/h = 10 the shear
term is 4.6 %, so a Kirchhoff-only reference would put the solid at +3.5 % and fail the 2 %
gate: the case therefore also checks that transverse shear deformation is present.

Remaining gap between the 3D solid and the Mindlin reference: plate theory clamps only the
mid-surface displacement and the rotation; fixing every node of the cylindrical face is a
stiffer, slightly different boundary condition, and Mindlin also neglects sigma_z and uses a
shear-correction factor. All of these are O((h/a)^2) effects. Measured, the solid converges
towards about -0.3 % of the Mindlin value at a/h = 10 (the same coefficient, w0 / w_K - 1 of
about 4.25 (h/a)^2 instead of 4.57 (h/a)^2, was seen at a/h = 20). Using kappa = pi^2/12
instead of 5/6 would change the reference by only +0.06 %.

Measured through Tessera's public API (gmsh mesher, CalculiX 2.21, October 2026):

    size [mm]  el. through h  nodes   w0 mid-surface   runtime
    2.0             1          1597      -1.89 %        ~0.3 s
    1.0             2          6771      -0.93 %        ~2 s
    0.6667          3         20193      -0.65 %        ~9 s   (gate)
    0.5             4         ~40 k      -0.57 %        ~31 s  (not in the sweep)
    0.4             5         ~75 k      -0.52 %        ~98 s  (not in the sweep)

(The last two rows come from the same problem scaled to a = 10, h = 1; the deflection error
depends only on a/h, nu and size/h.) The top/bottom mean sits 0.40 % below the mid-surface node
at every size.

Convergence is roughly first order in the element size (singular stresses at the corners of the
clamped face), from below.
"""

from __future__ import annotations

import itertools
import math
import shutil
from collections.abc import Callable
from pathlib import Path
from typing import Any

import numpy as np
import pytest

pytestmark = [
    pytest.mark.benchmark,
    pytest.mark.solver,
    pytest.mark.skipif(shutil.which("ccx") is None, reason="CalculiX (ccx) not on PATH"),
]

REFERENCE = "clamped_plate"

# Geometry (mm), load (MPa), material (MPa)
RADIUS = 20.0  # a
THICKNESS = 2.0  # h
PRESSURE = 1.0  # q, on the top face, acting in -z
YOUNGS_MODULUS = 210_000.0  # steel-generic
POISSON_RATIO = 0.30
SHEAR_CORRECTION = 5.0 / 6.0  # kappa (Mindlin/Reissner)

GATE_MESH_SIZE = THICKNESS / 3  # mm: 3 quadratic elements through the thickness
CONVERGENCE_MESH_SIZES = (THICKNESS, THICKNESS / 2, THICKNESS / 3)  # mm, coarse to fine
CONVERGENCE_SLACK = 0.002  # allowed relative-error increase between refinements (noise)


def flexural_rigidity(
    thickness: float = THICKNESS,
    youngs_modulus: float = YOUNGS_MODULUS,
    poisson_ratio: float = POISSON_RATIO,
) -> float:
    """D (N mm) = E h^3 / (12 (1 - nu^2))."""
    return youngs_modulus * thickness**3 / (12.0 * (1.0 - poisson_ratio**2))


def kirchhoff_deflection(
    r: float, radius: float = RADIUS, pressure: float = PRESSURE, thickness: float = THICKNESS
) -> float:
    """Thin-plate deflection (mm) w_K(r) = q (a^2 - r^2)^2 / (64 D) of the clamped plate."""
    return pressure * (radius**2 - r**2) ** 2 / (64.0 * flexural_rigidity(thickness))


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
    """Quarter disc x, y >= 0, r <= a, 0 <= z <= h as STEP (gmsh OpenCASCADE)."""
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


def _solve_clamped_plate(
    step_path: Path, mesh_size: float, workdir: Path
) -> tuple[float, float, int]:
    """Run the quarter disc through the public API.

    Returns (mid-surface centre deflection, top/bottom mean centre deflection, number of nodes),
    deflections positive in the load direction (-z).
    """
    return _centre_deflections(_solve_disc(step_path, mesh_size, workdir))


def _solve_disc(
    step_path: Path, mesh_size: float, workdir: Path, *, check_mesh: bool = False
) -> Any:
    """Run the quarter disc through the public API; return the ``StudyResult``."""
    from tessera.core import describe_faces, new_project, solve
    from tessera.model import FaceSet, PartInstance, PressureLoad, Support

    p = new_project("clamped-plate")
    p.assembly.parts.append(
        PartInstance(id="disc", name="disc", geometry=step_path, material_id="steel-generic")
    )
    faces = describe_faces(p, "disc")
    top = _face_with_normal(faces, (0.0, 0.0, 1.0))
    sym_x = _face_with_normal(faces, (-1.0, 0.0, 0.0))  # x = 0 cut plane
    sym_y = _face_with_normal(faces, (0.0, -1.0, 0.0))  # y = 0 cut plane
    rim = _curved_face_with_area(faces, math.pi / 2 * RADIUS * THICKNESS)  # 2 pi a h / 4
    assert next(f for f in faces if f.index == rim).kind.lower() == "cylinder"
    assert next(f for f in faces if f.index == top).centroid[2] == pytest.approx(THICKNESS)

    p.study.supports += [
        Support(id="clamp", faces=FaceSet(part_id="disc", faces=(rim,))),
        Support(id="sym-x", faces=FaceSet(part_id="disc", faces=(sym_x,)), dofs=(1,)),
        Support(id="sym-y", faces=FaceSet(part_id="disc", faces=(sym_y,)), dofs=(2,)),
    ]
    p.study.loads.append(
        PressureLoad(id="q", faces=FaceSet(part_id="disc", faces=(top,)), pressure=PRESSURE)
    )
    p.study.mesh_size = mesh_size

    workdir.mkdir(parents=True, exist_ok=True)
    return solve(p, workdir=workdir, check_mesh=check_mesh)


def _centre_deflections(result: Any) -> tuple[float, float, int]:
    """(mid-surface, top/bottom mean centre deflection, number of nodes) of a solved disc."""
    nodes = result.mesh.nodes
    u = result.solution.displacements
    tol = 1e-6 * RADIUS
    on_rim = np.abs(np.hypot(nodes[:, 0], nodes[:, 1]) - RADIUS) < tol
    assert np.abs(u[on_rim]).max() < 1e-12, "outer cylindrical face is not clamped"
    for axis in (0, 1):
        on_cut = np.abs(nodes[:, axis]) < tol
        assert np.abs(u[on_cut, axis]).max() < 1e-12, f"cut plane {'xy'[axis]} = 0 moves normal"

    on_axis = (np.abs(nodes[:, 0]) < tol) & (np.abs(nodes[:, 1]) < tol)
    z = nodes[on_axis, 2]
    w = -u[on_axis, 2]
    mid = np.abs(z - THICKNESS / 2) < tol
    assert mid.sum() == 1, "no mesh node at the mid-surface centre (x = y = 0, z = h/2)"
    surfaces = (np.abs(z) < tol) | (np.abs(z - THICKNESS) < tol)
    assert surfaces.sum() == 2, "missing top or bottom centre node"
    w_mid = float(w[mid][0])
    assert w_mid > 0.0, "plate did not deflect in the load direction: pressure sign or face"
    return w_mid, float(w[surfaces].mean()), result.mesh.n_nodes


@pytest.fixture(scope="module")
def disc_step(tmp_path_factory: pytest.TempPathFactory) -> Path:
    directory = tmp_path_factory.mktemp("clamped-plate-geometry")
    return write_quarter_disc_step(directory / "clamped_quarter_disc.step")


@pytest.fixture(scope="module")
def clamped_plate(
    disc_step: Path, tmp_path_factory: pytest.TempPathFactory
) -> Callable[[float], tuple[float, float, int]]:
    """Cached ``mesh_size -> (w0 mid, w0 top/bottom mean, n_nodes)`` shared by gate and sweep."""
    cache: dict[float, tuple[float, float, int]] = {}

    def run(mesh_size: float) -> tuple[float, float, int]:
        if mesh_size not in cache:
            workdir = tmp_path_factory.mktemp(f"clamped-plate-h{mesh_size:.4g}")
            cache[mesh_size] = _solve_clamped_plate(disc_step, mesh_size, workdir)
        return cache[mesh_size]

    return run


def test_closed_form_is_consistent() -> None:
    """Kirchhoff plate equation and clamp conditions; Mindlin shear term from equilibrium."""
    a, q = RADIUS, PRESSURE
    d = flexural_rigidity()
    # w_K(r) as a polynomial in r; axisymmetric Laplacian lap(f) = f'' + f'/r.
    w_k = np.polynomial.Polynomial([a**4, 0.0, -2.0 * a**2, 0.0, 1.0]) * (q / (64.0 * d))
    for r in (0.3 * a, 0.7 * a, a):
        assert math.isclose(w_k(r), kirchhoff_deflection(r), rel_tol=1e-12, abs_tol=1e-15)

    def laplacian(f: np.polynomial.Polynomial) -> np.polynomial.Polynomial:
        slope = f.deriv()
        assert slope.coef[0] == 0.0  # even in r: f'/r is again a polynomial
        return f.deriv(2) + np.polynomial.Polynomial(slope.coef[1:])

    biharmonic = laplacian(laplacian(w_k))
    for r in (0.0, 0.5 * a, a):
        assert math.isclose(d * biharmonic(r), q, rel_tol=1e-12)
    assert math.isclose(w_k(a), 0.0, abs_tol=1e-15)
    assert math.isclose(w_k.deriv()(a), 0.0, abs_tol=1e-15)

    # Shear part: kappa G h dw_S/dr = Q_r = -q r / 2 (equilibrium of the central disc), w_S(a) = 0.
    kappa_g_h = SHEAR_CORRECTION * YOUNGS_MODULUS / (2.0 * (1.0 + POISSON_RATIO)) * THICKNESS
    step = 1e-4 * a
    for r in (0.25 * a, 0.5 * a, 0.75 * a):
        slope = (shear_deflection(r + step) - shear_deflection(r - step)) / (2.0 * step)
        assert math.isclose(kappa_g_h * slope, -q * r / 2.0, rel_tol=1e-6)
    assert math.isclose(shear_deflection(a), 0.0, abs_tol=1e-15)

    ratio = shear_deflection(0.0) / kirchhoff_deflection(0.0)
    expected = 16.0 / (5.0 * (1.0 - POISSON_RATIO)) * (THICKNESS / RADIUS) ** 2
    assert math.isclose(ratio, expected, rel_tol=1e-12)


def test_reference_matches_formula(load_reference: Callable[[str], dict[str, Any]]) -> None:
    """The stored reference can never drift from the closed-form value."""
    ref = load_reference(REFERENCE)
    assert math.isclose(ref["reference_centre_deflection_mm"], centre_deflection(), rel_tol=1e-9)
    assert 0.0 < ref["rel_tolerance"] <= 0.02  # spec gate: within 2 %


def test_clamped_plate_centre_deflection(
    clamped_plate: Callable[[float], tuple[float, float, int]],
    load_reference: Callable[[str], dict[str, Any]],
) -> None:
    # Solve first: fail (not skip) if the solve breaks, before the reference check.
    w0, _, _ = clamped_plate(GATE_MESH_SIZE)

    ref = load_reference(REFERENCE)
    assert w0 == pytest.approx(ref["reference_centre_deflection_mm"], rel=ref["rel_tolerance"])


def test_clamped_plate_mesh_convergence(
    clamped_plate: Callable[[float], tuple[float, float, int]],
    load_reference: Callable[[str], dict[str, Any]],
) -> None:
    # Solve first: fail (not skip) if the solve breaks, before the reference check.
    for size in CONVERGENCE_MESH_SIZES:
        clamped_plate(size)

    ref = load_reference(REFERENCE)
    expected = ref["reference_centre_deflection_mm"]
    tolerance = ref["rel_tolerance"]

    rows = []
    for size in CONVERGENCE_MESH_SIZES:
        w_mid, w_surf, n_nodes = clamped_plate(size)
        rows.append((size, n_nodes, w_mid, abs(w_mid - expected) / expected, w_surf))

    print(
        f"\nClamped circular plate a/h = {RADIUS / THICKNESS:g}, q = {PRESSURE:g} MPa: "
        f"w0 ref {expected:.6e} mm (Kirchhoff {kirchhoff_deflection(0.0):.6e} + shear "
        f"{shear_deflection(0.0):.6e})"
    )
    print(
        f"{'h [mm]':>8} {'nodes':>8} {'w0 mid [mm]':>13} {'rel. error':>11}"
        f" {'top/bot mean':>13} {'vs mid':>8}"
    )
    for size, n_nodes, w_mid, error, w_surf in rows:
        print(
            f"{size:>8.4g} {n_nodes:>8d} {w_mid:>13.6e} {error:>11.4%}"
            f" {w_surf:>13.6e} {w_surf / w_mid - 1.0:>+8.3%}"
        )

    node_counts = [r[1] for r in rows]
    assert node_counts == sorted(node_counts) and len(set(node_counts)) == len(node_counts), (
        "mesh_size is not refining the mesh"
    )
    errors = [r[3] for r in rows]
    for coarse, fine in itertools.pairwise(errors):
        # Error must not grow under refinement, unless it already sits well inside the gate.
        assert fine <= coarse + CONVERGENCE_SLACK or fine < tolerance / 4
    assert errors[-1] < tolerance


THIN_WALL_MESH_SIZE = THICKNESS / 3  # mm: about 3 elements across the 2 mm plate


def test_clamped_plate_check_mesh_does_not_settle_the_clamp_on_a_thin_wall(
    disc_step: Path, tmp_path: Path
) -> None:
    """ADR-007 Amendment 2 (action item 11): at 2/3 mm the clamped rim's face peaks on the cut
    planes barely move (measured +4.1 % and +3.9 %, 3.0 coarse elements across the 2 mm plate)
    but keep rising at finer meshes (mesh-dependent from 0.47 mm). With fewer than 5 elements
    across, the guard reports them as not judged (thin-wall), never settled. The peaks sit
    0.23 mm from the rim, so the largest ball touching them there is 0.46 mm across (0.7 coarse
    elements), not the 2 mm thickness: a wall is never doubled across a symmetry plane, and the
    ball cannot cross the clamped rim. Measured: 37 s for both solves."""
    from tessera.model import Verdict

    result = _solve_disc(disc_step, THIN_WALL_MESH_SIZE, tmp_path, check_mesh=True)
    check = result.mesh_check
    assert check is not None
    judged = {Verdict.SETTLED, Verdict.MESH_DEPENDENT, Verdict.THIN_WALL}
    flagged = [s for s in result.face_peaks if s.flagged and s.face in check.face_peaks]
    assert flagged, "the clamp flags the face peaks at the rim"
    thin = []
    for spot in flagged:
        verdict = check.face_peaks[spot.face]
        assert spot.wall is not None and spot.wall.thickness <= THICKNESS * (1 + 1e-6)
        assert verdict.verdict is not Verdict.SETTLED or (
            verdict.elements_across is not None
            and verdict.elements_across >= check.wall_elements_min
        )
        assert verdict.verdict in judged | {Verdict.NOT_REFINED}
        if verdict.verdict is Verdict.THIN_WALL:
            assert verdict.change <= check.threshold
            thin.append(spot.face)
    assert thin, "a clamp peak that barely moved on this thin wall reads thin-wall"
    # Face 4, the y = 0 cut plane: its peak at the rim read settled (+4.1 %) without the guard.
    assert check.face_peaks[4].verdict is Verdict.THIN_WALL
