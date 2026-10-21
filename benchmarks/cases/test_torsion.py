"""Torsion of a circular shaft and a square bar (TS-09, Saint-Venant; gates 2 % and 3 %).

Problem (mm, N, MPa per ADR-002)
--------------------------------
Two steel bars along z, E = 210 000 MPa, nu = 0.30 (``steel-generic``),
G = E / (2 (1 + nu)) = 80 769.23 MPa, length L = 60:

(a) circular shaft, radius R = 5 (diameter 10);
(b) square bar, side a = 10, centred on the z axis (-a/2 <= x, y <= a/2).

In both, the root face z = 0 is fully clamped (``Support`` with the default dofs), the end face
z = L carries a torque T = 1e4 N mm about +z (``MomentLoad``, moment (0, 0, T)), and the lateral
surface is free. Linear static, quadratic tetrahedra (C3D10) in CalculiX. Whole bars, no
symmetry: torsion is antisymmetric under the mirror planes, so a reduced model would need
cyclic or antisymmetry constraints that Tessera does not offer.

``MomentLoad`` spreads the torque like a distributing coupling: nodal forces
F_i = w_i (omega x r_i), w_i = integral of N_i dA, r_i relative to the weighted face centroid,
omega = I^-1 M with I = sum w_i (|r_i|^2 1 - r_i r_i^T). On these flat end faces r_z = 0, so
I_xz = I_yz = 0 and omega = (0, 0, T / I_zz): nodal forces statically equivalent to the
traction (T / I_zz) e_z x r (lumped by area share, as CalculiX's *DISTRIBUTING does), with zero
net force.

Measured quantities
-------------------
- theta: twist per unit length. Every node gets the rotation about the axis
  phi = (x u_y - y u_x) / (x^2 + y^2) (nodes with x^2 + y^2 < (0.2 c)^2 are skipped, c = R or
  a/2). theta is the slope of a straight-line least-squares fit of phi against z over all nodes
  in the middle window 2 D <= z <= L - 2 D (D = 10 = diameter or side), i.e. 20 <= z <= 40. In
  Saint-Venant torsion the in-plane displacements are exactly u_x = -theta z y,
  u_y = theta z x for any section (warping only moves nodes along z), so phi = theta z + const
  at every node and the fit has no model bias.
- tau_max, from the nodal stresses (CalculiX extrapolates to nodes and averages) of surface
  nodes in the same window, as the tangential shear tau = -n_y tau_zx + n_x tau_yz with
  (n_x, n_y) the in-plane outward normal:
  (a) the mean over all lateral-surface nodes (the exact tau = T r / J is uniform there);
  (b) the four side faces pooled, with s the coordinate along the face from its midside: a
      least-squares fit tau = c0 + c2 s^2 over |s| <= a/4, tau_max = c0. The exact tau(s) is
      even in s with its maximum at s = 0, so no node needs to sit on the midside line. Fit
      bias on the exact series profile: +0.11 % (tau(a/4) = 0.833 tau_max).

Reference (a): circular shaft
-----------------------------
    J = pi R^4 / 2,   theta = T / (G J),   tau_max = T R / J

Source: Timoshenko, *Strength of Materials* Part I (torsion of a circular shaft; section and
edition **unverified**, as in CATALOG.md TS-09). Formula re-derived from the Saint-Venant field
u = -theta z y, v = theta z x, w = 0 and recomputed in ``test_reference_matches_formula``.
Numbers: J = 981.7477 mm^4, theta = 1.261113e-4 rad/mm, tau_max = 50.92958 MPa.

This case is exact everywhere, not only in the window: the Saint-Venant field vanishes at
z = 0 (the clamp is exact, a circle does not warp), and its end traction G theta e_z x r is
statically equivalent to the distributing coupling's nodal forces (lumped, not consistent, so
the end differs only within an element or so of z = L). The field is linear in x, y, z times
z (quadratic), which C3D10 holds exactly; only the faceted-curved boundary is approximate.

Reference (b): square bar
-------------------------
    J = k1 a^4,   theta = T / (G J),   tau_max = T / (k2 a^3)

    k1 = (1/3) [1 - (192 / pi^5) sum_{n odd} tanh(n pi / 2) / n^5]         = 0.1405770
    k  = 1 - (8 / pi^2) sum_{n odd} 1 / (n^2 cosh(n pi / 2)),  k2 = k1 / k  = 0.2081653

(tau_max = G theta a k at the midside.) Source: the Saint-Venant series for a rectangular bar,
Timoshenko & Goodier, *Theory of Elasticity*, ch. 10 (sec. 109 per Karami's course notes;
**unverified**, CATALOG.md TS-09). k1 = 0.1406 and k2 = 0.2082 are **computed** from the series
(CATALOG.md); the tabulated 0.141 / 0.208 are **verified** in the Atlantis Press and CliffsNotes
excerpts (CATALOG.md). ``test_saint_venant_constants`` recomputes both and checks them against
the catalog to its stated digits. Numbers: J = 1405.770 mm^4, theta = 8.807238e-5 rad/mm,
tau_max = 48.03876 MPa.

Why these choices
-----------------
- Ends (b). The clamp suppresses warping and the linear end traction is not the Saint-Venant
  shear distribution of a square, so both ends carry a self-equilibrated disturbance (zero net
  force and moment). By Saint-Venant's principle it decays along the bar, for a compact section
  within about one side length. The independent check below measured it: the local twist rate
  in 3 mm bins (h = 1.25, L = 60) is -5.5 %, -0.55 %, -0.07 % off the plateau's -0.05 % in the
  first three bins from the root, and +3.2 %, +0.6 %, +0.007 % in the last three at the loaded
  end. The disturbance is below 0.03 % beyond about 0.9 a from either end. The window keeps 2 a
  from both ends, and L/a = 6 leaves it 2 a long. With L = 100 the plateau is the same: theta
  -0.0474 % against -0.0475 % at h = 1.25.
- Ends (a). Nothing to avoid (exact everywhere); the same window is used for symmetry of the
  code. The twist-rate profile in 3 mm bins is flat to 0.1 % from end to end at h = 2.5.
- tau_max (b) uses surface stresses at the midside, as the catalog asks. They converge more
  slowly than theta (-7.1 % at h = 5, -0.78 % at h = 2.5), hence the catalog's looser 3 %.
- Torque size. T = 1e4 N mm gives tau about 50 MPa and a total twist below 0.01 rad, so the
  rotations are small and the measure needs no geometric correction.
- Gate mesh h = 2.5 (2 elements per radius, 4 across the square). It passes with margin and
  keeps both gate solves to a few seconds; the convergence sweep adds h = 5 and h = 1.25.

Check (Tessera's mesher and deck writer, the distributing-coupling forces
F_i = w_i omega x r_i computed by hand from ``tri6_load_weights`` and added as ``NodalLoads``;
not independent of those weights, see the *DISTRIBUTING run below for that,
CalculiX 2.21, October 2026), measured with ``_measure`` from this file, errors against the
references:

    case        h [mm]   nodes   theta      tau_max
    (a) circle  5          583   -0.11 %    +0.76 %
                2.5       3139   +0.004 %   +0.22 %
                1.25     18602   +0.0004 %  +0.04 %    (~8 s mesh + ccx)
    (b) square  5          663   -3.50 %    -7.08 %
                2.5       4046   -0.42 %    -0.78 %
                1.25     23881   -0.048 %   -0.10 %    (~10 s mesh + ccx)

At h = 2.5 and 1.25 a second run with CalculiX's own ``*COUPLING`` / ``*DISTRIBUTING`` on the
end face and ``*CLOAD`` moment 6 on a reference node gave the same theta to within 0.0001 %
(circle +0.0039 %, +0.0004 %; square -0.4241 %, -0.0477 %). Net force of the hand forces:
below 1e-12 N; their moment: T to machine precision.
"""

from __future__ import annotations

import itertools
import math
import shutil
from collections.abc import Callable
from pathlib import Path
from typing import Any, Literal

import numpy as np
import pytest
from numpy.typing import NDArray

pytestmark = [
    pytest.mark.benchmark,
    pytest.mark.solver,
    pytest.mark.skipif(shutil.which("ccx") is None, reason="CalculiX (ccx) not on PATH"),
]

REFERENCE = "torsion"

Shape = Literal["circle", "square"]
SHAPES: tuple[Shape, ...] = ("circle", "square")

# Geometry (mm), load (N mm), material (MPa)
LENGTH = 60.0  # L, along z; root (clamped) at z = 0, torque on z = L
RADIUS = 5.0  # R, circular shaft
SIDE = 10.0  # a, square bar
TORQUE = 1.0e4  # T, about +z
YOUNGS_MODULUS = 210_000.0  # steel-generic
POISSON_RATIO = 0.30
SHEAR_MODULUS = YOUNGS_MODULUS / (2.0 * (1.0 + POISSON_RATIO))

# Catalog values (CATALOG.md TS-09), to the digits it states.
CATALOG_K1 = 0.1406
CATALOG_K2 = 0.2082
SERIES_TERMS = 2000  # odd n = 1 .. 3999; the k1 tail (~1/n^4) is below 1e-15 relative

WINDOW_MARGIN = 2.0 * SIDE  # middle window 2 D <= z <= L - 2 D (D = 2 R = a = 10)
AXIS_EXCLUSION = 0.2  # skip nodes closer to the axis than 0.2 R (or 0.2 a/2) for phi
MIDSIDE_BAND = SIDE / 4  # |s| <= a/4 for the square's midside tau fit

GATE_MESH_SIZE = 2.5  # mm: 2 elements per radius, 4 across the square
CONVERGENCE_MESH_SIZES = (5.0, 2.5, 1.25)  # mm, coarse to fine
CONVERGENCE_SLACK = 0.002  # allowed relative-error increase between refinements (noise)


def _odd(terms: int = SERIES_TERMS) -> range:
    return range(1, 2 * terms, 2)


def square_torsion_k1(terms: int = SERIES_TERMS) -> float:
    """Saint-Venant J / a^4 of a square: (1/3)[1 - 192/pi^5 sum tanh(n pi/2)/n^5], n odd."""
    total = sum(math.tanh(n * math.pi / 2) / n**5 for n in _odd(terms))
    return (1.0 - 192.0 / math.pi**5 * total) / 3.0


def square_torsion_k2(terms: int = SERIES_TERMS) -> float:
    """Saint-Venant T / (tau_max a^3) of a square: k1 / k, k = tau_max / (G theta a)."""
    # sech(x) = 2 e^-x / (1 + e^-2x): math.cosh overflows for n above about 450.
    total = sum(
        2.0 * math.exp(-n * math.pi / 2) / (1.0 + math.exp(-n * math.pi)) / n**2
        for n in _odd(terms)
    )
    return square_torsion_k1(terms) / (1.0 - 8.0 / math.pi**2 * total)


def torsion_constant(shape: Shape) -> float:
    """J (mm^4): pi R^4 / 2 for the circle, k1 a^4 for the square."""
    if shape == "circle":
        return math.pi * RADIUS**4 / 2.0
    return square_torsion_k1() * SIDE**4


def twist_rate(shape: Shape, torque: float = TORQUE, shear_modulus: float = SHEAR_MODULUS) -> float:
    """theta (rad/mm) = T / (G J)."""
    return torque / (shear_modulus * torsion_constant(shape))


def max_shear_stress(shape: Shape, torque: float = TORQUE) -> float:
    """tau_max (MPa): T R / J for the circle, T / (k2 a^3) at the square's midsides."""
    if shape == "circle":
        return torque * RADIUS / torsion_constant(shape)
    return torque / (square_torsion_k2() * SIDE**3)


def write_bar_step(path: Path, shape: Shape) -> Path:
    """Circular (radius R) or square (side a, centred) bar 0 <= z <= L as STEP (gmsh OCC)."""
    import gmsh

    owns_session = not gmsh.isInitialized()
    if owns_session:
        gmsh.initialize(interruptible=False)
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add(f"benchmark-{path.stem}")
    try:
        if shape == "circle":
            gmsh.model.occ.addCylinder(0.0, 0.0, 0.0, 0.0, 0.0, LENGTH, RADIUS)
        else:
            gmsh.model.occ.addBox(-SIDE / 2, -SIDE / 2, 0.0, SIDE, SIDE, LENGTH)
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


def _measure(
    shape: Shape,
    nodes: NDArray[np.float64],
    u: NDArray[np.float64],
    s: NDArray[np.float64],
) -> tuple[float, float]:
    """Return (theta, tau_max) from nodal displacements and stresses (see module docstring).

    ``s`` columns: sxx, syy, szz, sxy, syz, szx.
    """
    x, y, z = nodes[:, 0], nodes[:, 1], nodes[:, 2]
    half = RADIUS if shape == "circle" else SIDE / 2
    window = (z >= WINDOW_MARGIN) & (z <= LENGTH - WINDOW_MARGIN)
    rho2 = x**2 + y**2

    use = window & (rho2 > (AXIS_EXCLUSION * half) ** 2)
    assert np.unique(np.round(z[use], 6)).size >= 3, "too few node levels in the window"
    phi = (x[use] * u[use, 1] - y[use] * u[use, 0]) / rho2[use]
    theta = float(np.polyfit(z[use], phi, 1)[0])

    tau_zx, tau_yz = s[:, 5], s[:, 4]
    tol = 1e-6 * half
    if shape == "circle":
        surface = window & (np.abs(np.sqrt(rho2) - RADIUS) < tol)
        assert surface.sum() >= 10, "no lateral-surface nodes found in the window"
        r = np.sqrt(rho2[surface])
        tau = (-y[surface] * tau_zx[surface] + x[surface] * tau_yz[surface]) / r
        return theta, float(tau.mean())

    along: list[NDArray[np.float64]] = []
    values: list[NDArray[np.float64]] = []
    for n_x, n_y in ((1.0, 0.0), (-1.0, 0.0), (0.0, 1.0), (0.0, -1.0)):
        on_face = np.abs(n_x * x + n_y * y - SIDE / 2) < tol
        coord = -n_y * x + n_x * y  # along the face, 0 at its midside
        band = window & on_face & (np.abs(coord) <= MIDSIDE_BAND)
        assert band.sum() >= 6, f"too few nodes near the midside of face n = ({n_x}, {n_y})"
        along.append(coord[band])
        values.append(-n_y * tau_zx[band] + n_x * tau_yz[band])
    c2, c0 = np.polyfit(np.concatenate(along) ** 2, np.concatenate(values), 1)
    assert c2 < 0.0, "midside shear is not a maximum: wrong torque direction or no torsion?"
    return theta, float(c0)


def _solve_torsion(
    shape: Shape, step_path: Path, mesh_size: float, workdir: Path
) -> tuple[float, float, int]:
    """Run one clamped bar under an end torque through the public API.

    Returns (theta, tau_max, number of nodes).
    """
    from tessera.core import describe_faces, new_project, solve
    from tessera.model import FaceSet, MomentLoad, PartInstance, Support

    p = new_project(f"torsion-{shape}")
    p.assembly.parts.append(
        PartInstance(id="bar", name="bar", geometry=step_path, material_id="steel-generic")
    )
    faces = describe_faces(p, "bar")
    root = _face_with_normal(faces, (0.0, 0.0, -1.0))
    tip = _face_with_normal(faces, (0.0, 0.0, 1.0))
    tip_info = next(f for f in faces if f.index == tip)
    area = math.pi * RADIUS**2 if shape == "circle" else SIDE**2
    assert tip_info.area == pytest.approx(area, rel=1e-6)
    assert tip_info.centroid == pytest.approx((0.0, 0.0, LENGTH), abs=1e-6 * LENGTH)

    p.study.supports.append(Support(id="root", faces=FaceSet(part_id="bar", faces=(root,))))
    p.study.loads.append(
        MomentLoad(id="torque", faces=FaceSet(part_id="bar", faces=(tip,)), moment=(0, 0, TORQUE))
    )
    p.study.mesh_size = mesh_size

    workdir.mkdir(parents=True, exist_ok=True)
    result = solve(p, workdir=workdir)

    nodes = result.mesh.nodes
    u = result.solution.displacements
    theta, tau_max = _measure(shape, nodes, u, result.solution.stresses)
    # First, so that the failure is clear if the moment is dropped or applied the wrong way.
    assert theta > 0.0, f"bar did not twist about +z (theta = {theta:.3e} rad/mm)"
    on_root = np.abs(nodes[:, 2]) < 1e-6 * LENGTH
    assert np.abs(u[on_root]).max() < 1e-12, "root face is not clamped"
    return theta, tau_max, result.mesh.n_nodes


@pytest.fixture(scope="module")
def bar_steps(tmp_path_factory: pytest.TempPathFactory) -> dict[Shape, Path]:
    directory = tmp_path_factory.mktemp("torsion-geometry")
    return {shape: write_bar_step(directory / f"torsion_{shape}.step", shape) for shape in SHAPES}


@pytest.fixture(scope="module")
def torsion(
    bar_steps: dict[Shape, Path], tmp_path_factory: pytest.TempPathFactory
) -> Callable[[Shape, float], tuple[float, float, int]]:
    """Cached ``(shape, mesh_size) -> (theta, tau_max, n_nodes)`` shared by gate and sweep."""
    cache: dict[tuple[Shape, float], tuple[float, float, int]] = {}

    def run(shape: Shape, mesh_size: float) -> tuple[float, float, int]:
        key = (shape, mesh_size)
        if key not in cache:
            workdir = tmp_path_factory.mktemp(f"torsion-{shape}-h{mesh_size:g}")
            cache[key] = _solve_torsion(shape, bar_steps[shape], mesh_size, workdir)
        return cache[key]

    return run


def test_saint_venant_constants() -> None:
    """The series reproduce the catalog's k1 and k2 to its stated (4) digits, and the
    tabulated 0.141 / 0.208 to 3; the truncated series has converged."""
    k1, k2 = square_torsion_k1(), square_torsion_k2()
    assert round(k1, 4) == CATALOG_K1 and round(k1, 3) == 0.141
    assert round(k2, 4) == CATALOG_K2 and round(k2, 3) == 0.208
    assert math.isclose(k1, square_torsion_k1(4 * SERIES_TERMS), rel_tol=1e-12)
    assert math.isclose(k2, square_torsion_k2(4 * SERIES_TERMS), rel_tol=1e-12)


def test_reference_matches_formula(load_reference: Callable[[str], dict[str, Any]]) -> None:
    """The stored reference can never drift from the closed-form value."""
    ref = load_reference(REFERENCE)
    for shape in SHAPES:
        assert math.isclose(ref[f"reference_twist_rate_{shape}"], twist_rate(shape), rel_tol=1e-9)
        assert math.isclose(
            ref[f"reference_max_shear_{shape}_mpa"], max_shear_stress(shape), rel_tol=1e-9
        )
    assert 0.0 < ref["rel_tolerance_twist_rate"] <= 0.02  # spec gate: within 2 %
    assert 0.0 < ref["rel_tolerance_max_shear"] <= 0.03  # spec gate: within 3 %


@pytest.mark.parametrize("shape", SHAPES)
def test_torsion(
    shape: Shape,
    torsion: Callable[[Shape, float], tuple[float, float, int]],
    load_reference: Callable[[str], dict[str, Any]],
) -> None:
    # Solve first: fail (not skip) while moment loads are missing, before the reference check.
    theta, tau_max, _ = torsion(shape, GATE_MESH_SIZE)

    ref = load_reference(REFERENCE)
    assert theta == pytest.approx(
        ref[f"reference_twist_rate_{shape}"], rel=ref["rel_tolerance_twist_rate"]
    )
    assert tau_max == pytest.approx(
        ref[f"reference_max_shear_{shape}_mpa"], rel=ref["rel_tolerance_max_shear"]
    )


@pytest.mark.parametrize("shape", SHAPES)
def test_torsion_mesh_convergence(
    shape: Shape,
    torsion: Callable[[Shape, float], tuple[float, float, int]],
    load_reference: Callable[[str], dict[str, Any]],
) -> None:
    # Solve first: fail (not skip) while moment loads are missing, before the reference check.
    for size in CONVERGENCE_MESH_SIZES:
        torsion(shape, size)

    ref = load_reference(REFERENCE)
    theta_ref = ref[f"reference_twist_rate_{shape}"]
    tau_ref = ref[f"reference_max_shear_{shape}_mpa"]

    rows = []
    for size in CONVERGENCE_MESH_SIZES:
        theta, tau_max, n_nodes = torsion(shape, size)
        rows.append(
            (
                size,
                n_nodes,
                theta,
                abs(theta - theta_ref) / theta_ref,
                tau_max,
                abs(tau_max - tau_ref) / tau_ref,
            )
        )

    print(
        f"\nTorsion, {shape} bar, T = {TORQUE:g} N mm: theta ref {theta_ref:.6e} rad/mm, "
        f"tau_max ref {tau_ref:.5f} MPa"
    )
    print(
        f"{'h [mm]':>8} {'nodes':>8} {'theta':>13} {'rel. error':>11}"
        f" {'tau_max':>10} {'rel. error':>11}"
    )
    for size, n_nodes, theta, e_theta, tau_max, e_tau in rows:
        print(
            f"{size:>8g} {n_nodes:>8d} {theta:>13.6e} {e_theta:>11.4%}"
            f" {tau_max:>10.5f} {e_tau:>11.4%}"
        )

    node_counts = [r[1] for r in rows]
    assert node_counts == sorted(node_counts) and len(set(node_counts)) == len(node_counts), (
        "mesh_size is not refining the mesh"
    )
    for column, tolerance in (
        (3, ref["rel_tolerance_twist_rate"]),
        (5, ref["rel_tolerance_max_shear"]),
    ):
        errors = [r[column] for r in rows]
        for coarse, fine in itertools.pairwise(errors):
            # Error must not grow under refinement, unless it already sits well inside the gate.
            assert fine <= coarse + CONVERGENCE_SLACK or fine < tolerance / 4
        assert errors[-1] < tolerance
