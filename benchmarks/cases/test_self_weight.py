"""Bar hanging under its own weight: elongation and axial stress (TS-04, v0.1 gate: 2 %).

Problem (mm, N, MPa, tonne, s per ADR-002)
------------------------------------------
Steel prism along z, length L = 200, square cross-section 20 x 20 (half-width c = 10),
E = 210 000 MPa, nu = 0.30, rho = 7.85e-9 t/mm^3 (material ``steel-generic``). It hangs from its
top face z = L under gravity g = 9810 mm/s^2 in -z (``Study.gravity = (0, 0, -9810)``); the
bottom face z = 0 is free. Linear static, quadratic tetrahedra (C3D10) in CalculiX.

Modelled as one quarter (0 <= x, y <= c) with symmetry rollers u_x = 0 on x = 0 and u_y = 0 on
y = 0, and a roller u_z = 0 on the top face (see "Why a roller top face" below). No other
load: the only load is the self-weight.

Measured quantities
-------------------
- delta: elongation = -(mean u_z of the bottom-face nodes), since the top face has u_z = 0.
- sigma_mid: sigma_zz at mid-length z = L/2, from a straight-line least-squares fit of the
  nodal sigma_zz against z over the nodes with |z - L/2| <= L/20 (no node is guaranteed to sit
  exactly at z = L/2; the exact field is linear in z, so the fit adds no bias).

Reference: bar under its own weight
-----------------------------------
With z measured up from the free end, the weight below a section gives

    sigma_zz(z) = rho g z,     delta = integral_0^L sigma_zz / E dz = rho g L^2 / (2 E)

    sigma_mid = rho g L / 2

Source: elementary strength of materials, S. P. Timoshenko, *Strength of Materials* Part I,
3rd ed., Van Nostrand 1955, chapter on tension and compression (section **unverified**, as in
CATALOG.md TS-04). Verification status: **unverified** in the book; the formula is re-derived
above and recomputed in ``test_reference_matches_formula``.

Numbers: rho g = 7.700850e-5 N/mm^3, delta = 7.334143e-6 mm, sigma_mid = 7.700850e-3 MPa.

The values are small but well resolved: the CalculiX ``.frd`` stores six significant digits
in floating point (relative error <= 5e-6), so there is no absolute noise floor. The
independent check below reproduces both values to 3e-5 relative. Physical g and density are
kept on purpose: the case exists to catch unit slips in the density/gravity chain (ADR-002).
A slip shows up as a ratio delta_FE / delta_ref of about
1e3 (kg/mm^3), 1e9 (g/cm^3 read as t/mm^3), 1e12 (kg/m^3 read as t/mm^3) or 1e-3 (g in m/s^2),
all far outside any tolerance.

Why a roller top face, not a clamped one
----------------------------------------
The 3D elasticity solution with a uniform support traction on the top face is exact and simple
(stresses linear in z, so compatibility holds; it can be checked by substitution):

    sigma_zz = rho g z,  all other stresses 0
    u = -nu rho g x z / E,  v = -nu rho g y z / E,
    w = rho g z^2 / (2 E) + nu rho g (x^2 + y^2) / (2 E) + const

so every fibre stretches by exactly rho g L^2 / (2 E), whatever its (x, y). The field is
quadratic, which C3D10 represents exactly, and u = 0 on x = 0, v = 0 on y = 0 (the symmetry
rollers are exact). The only departure is at the top face: the roller forces it flat, while the
exact solution bows it by nu rho g c^2 / (2 E). That perturbation is self-equilibrated and
decays within about c of the top (Saint-Venant). Its effect on delta is of order
nu (c/L)^2 < 0.1 %.

A fully clamped top face (as worded in CATALOG.md) also suppresses the lateral Poisson
contraction u = -nu rho g x L / E at the top. That gives a mesh-independent bias of order
nu^2 c / L. The independent check measured -0.61 % on this geometry (-0.68 % at h = 10,
-0.61 % at h = 2.5), against -0.003 % for the roller model. The roller model therefore isolates
what this case tests, the gravity load and the density units, from a support artefact.

Independent check (Tessera's mesher and deck writer, gravity added by hand as CalculiX
``*DLOAD GRAV``, CalculiX 2.21, October 2026):
h = 10: delta -0.008 %, sigma_mid -0.001 % (898 nodes);
h = 5: -0.004 %, -0.001 % (1839 nodes); h = 2.5: -0.003 %, -0.001 % (13181 nodes, ~3 s).
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

pytestmark = [
    pytest.mark.benchmark,
    pytest.mark.solver,
    pytest.mark.skipif(shutil.which("ccx") is None, reason="CalculiX (ccx) not on PATH"),
]

REFERENCE = "self_weight_bar"

# Geometry (mm), load (mm/s^2), material (MPa, t/mm^3)
LENGTH = 200.0  # along z; top (supported) face at z = L, free end at z = 0
HALF_WIDTH = 10.0  # c: quarter model 0 <= x, y <= c of a 20 x 20 section
GRAVITY = 9810.0  # magnitude, applied in -z
YOUNGS_MODULUS = 210_000.0  # steel-generic
POISSON_RATIO = 0.30
DENSITY = 7.85e-9  # steel-generic, t/mm^3

GATE_MESH_SIZE = 5.0  # mm: 2 elements across the quarter section
CONVERGENCE_MESH_SIZES = (10.0, 5.0, 2.5)  # mm, coarse to fine
CONVERGENCE_SLACK = 0.002  # allowed relative-error increase between refinements (noise)
MID_BAND = LENGTH / 20  # half-width of the z band used for the mid-length stress fit


def self_weight_elongation(
    length: float = LENGTH,
    density: float = DENSITY,
    gravity: float = GRAVITY,
    youngs_modulus: float = YOUNGS_MODULUS,
) -> float:
    """Elongation (mm) of a prismatic bar hanging under its own weight: rho g L^2 / (2 E)."""
    return density * gravity * length**2 / (2.0 * youngs_modulus)


def self_weight_mid_stress(
    length: float = LENGTH, density: float = DENSITY, gravity: float = GRAVITY
) -> float:
    """Axial stress (MPa) at mid-length: rho g L / 2."""
    return density * gravity * length / 2.0


def _face_with_normal(faces: list[Any], normal: tuple[float, float, float]) -> int:
    matches = [
        f.index for f in faces if f.normal is not None and np.allclose(f.normal, normal, atol=1e-6)
    ]
    assert len(matches) == 1, f"expected one face with normal {normal}, got {matches}"
    return matches[0]


def _solve_self_weight(
    step_path: Path, mesh_size: float, workdir: Path
) -> tuple[float, float, int]:
    """Run the quarter bar through the public API.

    Returns (elongation delta, sigma_zz at mid-length, number of nodes).
    """
    from tessera.core import describe_faces, new_project, solve
    from tessera.model import FaceSet, PartInstance, Support

    p = new_project("self-weight")
    p.assembly.parts.append(
        PartInstance(id="bar", name="bar", geometry=step_path, material_id="steel-generic")
    )
    faces = describe_faces(p, "bar")
    top = _face_with_normal(faces, (0.0, 0.0, 1.0))
    sym_x = _face_with_normal(faces, (-1.0, 0.0, 0.0))
    sym_y = _face_with_normal(faces, (0.0, -1.0, 0.0))
    top_info = next(f for f in faces if f.index == top)
    assert top_info.centroid[2] == pytest.approx(LENGTH)
    assert top_info.area == pytest.approx(HALF_WIDTH * HALF_WIDTH)

    rollers: tuple[tuple[str, int, Literal[1, 2, 3]], ...] = (
        ("top", top, 3),
        ("sym-x", sym_x, 1),
        ("sym-y", sym_y, 2),
    )
    for name, face, dof in rollers:
        p.study.supports.append(
            Support(id=name, faces=FaceSet(part_id="bar", faces=(face,)), dofs=(dof,))
        )
    p.study.gravity = (0.0, 0.0, -GRAVITY)
    p.study.mesh_size = mesh_size

    workdir.mkdir(parents=True, exist_ok=True)
    result = solve(p, workdir=workdir)

    nodes = result.mesh.nodes
    u = result.solution.displacements
    tol = 1e-6 * LENGTH
    on_bottom = np.abs(nodes[:, 2]) < tol
    assert on_bottom.sum() >= 6, "no mesh nodes found on the free bottom face"
    delta = -float(u[on_bottom, 2].mean())
    # First, so that the failure is clear while gravity is still ignored by the physics.
    assert delta > 0.0, f"free end did not move down (mean u_z = {-delta:.3e} mm): no gravity?"

    on_top = np.abs(nodes[:, 2] - LENGTH) < tol
    assert np.abs(u[on_top, 2]).max() < 1e-12, "top face moves along z"
    # Roller, not clamp: the top face must still contract laterally (Poisson).
    assert np.abs(u[on_top, :2]).max() > 1e-3 * delta, "top face is clamped, expected a roller"
    assert np.abs(u[np.abs(nodes[:, 0]) < tol, 0]).max() < 1e-12, "x = 0 plane moves along x"
    assert np.abs(u[np.abs(nodes[:, 1]) < tol, 1]).max() < 1e-12, "y = 0 plane moves along y"

    band = np.abs(nodes[:, 2] - LENGTH / 2) <= MID_BAND
    assert np.unique(np.round(nodes[band, 2], 9)).size >= 2, "too few node levels near L/2"
    slope, intercept = np.polyfit(nodes[band, 2], result.solution.stresses[band, 2], 1)
    sigma_mid = float(slope * LENGTH / 2 + intercept)
    return delta, sigma_mid, result.mesh.n_nodes


@pytest.fixture(scope="module")
def bar_step(box_step: Callable[..., Path]) -> Path:
    return box_step("self_weight_quarter", (0.0, 0.0, 0.0), (HALF_WIDTH, HALF_WIDTH, LENGTH))


@pytest.fixture(scope="module")
def self_weight(
    bar_step: Path, tmp_path_factory: pytest.TempPathFactory
) -> Callable[[float], tuple[float, float, int]]:
    """Cached ``mesh_size -> (delta, sigma_mid, n_nodes)`` so the gate and the sweep share runs."""
    cache: dict[float, tuple[float, float, int]] = {}

    def run(mesh_size: float) -> tuple[float, float, int]:
        if mesh_size not in cache:
            workdir = tmp_path_factory.mktemp(f"self-weight-h{mesh_size:g}")
            cache[mesh_size] = _solve_self_weight(bar_step, mesh_size, workdir)
        return cache[mesh_size]

    return run


def test_reference_matches_formula(load_reference: Callable[[str], dict[str, Any]]) -> None:
    """The stored reference can never drift from the closed-form value."""
    ref = load_reference(REFERENCE)
    assert math.isclose(ref["reference_elongation_mm"], self_weight_elongation(), rel_tol=1e-9)
    assert math.isclose(ref["reference_mid_stress_mpa"], self_weight_mid_stress(), rel_tol=1e-9)
    assert 0.0 < ref["rel_tolerance"] <= 0.02  # spec gate: within 2 %


def test_self_weight(
    self_weight: Callable[[float], tuple[float, float, int]],
    load_reference: Callable[[str], dict[str, Any]],
) -> None:
    # Solve first: fail (not skip) while gravity is ignored, before the reference check can skip.
    delta, sigma_mid, _ = self_weight(GATE_MESH_SIZE)

    ref = load_reference(REFERENCE)
    ratio = delta / ref["reference_elongation_mm"]
    assert delta == pytest.approx(ref["reference_elongation_mm"], rel=ref["rel_tolerance"]), (
        f"delta_FE / delta_ref = {ratio:.6g} (about 1e3, 1e9, 1e12 or 1e-3 means a unit slip)"
    )
    assert sigma_mid == pytest.approx(ref["reference_mid_stress_mpa"], rel=ref["rel_tolerance"])


def test_self_weight_mesh_convergence(
    self_weight: Callable[[float], tuple[float, float, int]],
    load_reference: Callable[[str], dict[str, Any]],
) -> None:
    # Solve first: fail (not skip) while gravity is ignored, before the reference check can skip.
    for size in CONVERGENCE_MESH_SIZES:
        self_weight(size)

    ref = load_reference(REFERENCE)
    delta_ref = ref["reference_elongation_mm"]
    sigma_ref = ref["reference_mid_stress_mpa"]
    tolerance = ref["rel_tolerance"]

    rows = []
    for size in CONVERGENCE_MESH_SIZES:
        delta, sigma_mid, n_nodes = self_weight(size)
        rows.append(
            (
                size,
                n_nodes,
                delta,
                abs(delta - delta_ref) / delta_ref,
                sigma_mid,
                abs(sigma_mid - sigma_ref) / sigma_ref,
            )
        )

    print(
        f"\nBar under its own weight, delta_ref {delta_ref:.6e} mm, "
        f"sigma_ref(L/2) {sigma_ref:.6e} MPa"
    )
    print(
        f"{'h [mm]':>8} {'nodes':>8} {'delta [mm]':>13} {'rel. error':>11}"
        f" {'sigma_mid':>13} {'rel. error':>11}"
    )
    for size, n_nodes, delta, e_delta, sigma_mid, e_sigma in rows:
        print(
            f"{size:>8g} {n_nodes:>8d} {delta:>13.6e} {e_delta:>11.4%}"
            f" {sigma_mid:>13.6e} {e_sigma:>11.4%}"
        )

    node_counts = [r[1] for r in rows]
    assert node_counts == sorted(node_counts) and len(set(node_counts)) == len(node_counts), (
        "mesh_size is not refining the mesh"
    )
    for column in (3, 5):
        errors = [r[column] for r in rows]
        for coarse, fine in itertools.pairwise(errors):
            # C3D10 holds the exact quadratic field, so the error is already ~0.01 % on the
            # coarsest mesh. It must not grow, unless it sits well inside the gate.
            assert fine <= coarse + CONVERGENCE_SLACK or fine < tolerance / 4
        assert errors[-1] < tolerance
