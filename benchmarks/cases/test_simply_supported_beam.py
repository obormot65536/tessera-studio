"""Simply supported beam under uniform pressure, midspan deflection (TS-05; gate 2 %).

Problem (mm, N, MPa per ADR-002)
--------------------------------
Prismatic steel beam along +x: span L = 200, width b = 10 (y), height h = 10 (z), L/h = 20,
E = 210 000 MPa, nu = 0.30 (``steel-generic``). The beam occupies 0 <= x <= L,
-b/2 <= y <= b/2, 0 <= z <= h. Uniform pressure p = 0.5 MPa acts on the top face z = h
(``PressureLoad``, positive = into the part, i.e. in -z), a line load q = p b = 5 N/mm.
Linear static, quadratic tetrahedra (C3D10) in CalculiX.

Supports are line supports (ADR-006) on the bottom edges of the two end faces, i.e. on the
lines x = 0, z = 0 and x = L, z = 0:

- ``pin``: the edge between the end face x = 0 and the bottom face, all three directions held;
- ``roller``: the edge between the end face x = L and the bottom face, u_z only.

Modelled as a half beam (0 <= y <= b/2) with a symmetry roller u_y = 0 on the cut plane y = 0.
Geometry, load and both supports are symmetric about y = 0, so the roller is exact. A quarter
model (also cutting at midspan x = L/2 with u_x = 0) is *not* exact with a pinned end: the pin
and the midspan plane would both hold u_x and stop the bottom fibre from stretching (arch
action). It would be exact only with rollers at both ends, which is not the ADR-006 set-up.

Rigid-body modes: the pinned line holds every translation and every rotation except the one
about the line itself (the y axis); the roller line, at a different x, removes that rotation.
So even the full-width beam needs no extra restraint, and the symmetry roller is there only to
halve the model. The pinned line also holds u_y along its whole length, which stops the lateral
Poisson strain of the bottom fibre at the pin, a slight local over-constraint; measured, holding
only (u_x, u_z) on the pinned line changes the midspan deflection by 0.001 % (2.5 and 1.67 mm
meshes), so the decided xyz pin is kept.

Measured quantity
-----------------
w_mid: -u_z of the mesh node at the centroid of the midspan section, x = L/2, y = 0, z = h/2,
i.e. on the neutral axis in the symmetry plane. Beam theory describes this centroidal line.

That node exists on every mesh because the geometry puts an edge there: the symmetry face
y = 0 is split along the neutral axis z = h/2 into two faces (a CAD split line, imprinted with
gmsh's ``fragment``; the part is still one solid). The edge between the two halves is a
straight geometry edge of length L, which gmsh splits uniformly; with quadratic elements an
edge split into n segments always has a node at its midpoint (a corner node for even n, a
midside node for odd n). Both halves carry the symmetry roller.

The bottom (z = 0) and top (z = h) nodes of the midspan section on y = 0 sit 0.09 % below the
neutral-axis node at every mesh (through-thickness Poisson strain of bending), so they are not
used for the gate.

Reference: Euler-Bernoulli bending + Timoshenko shear
-----------------------------------------------------
    w_mid = w_B + w_S
    w_B = 5 q L^4 / (384 E I),       I = b h^3 / 12
    w_S = q L^2 / (8 kappa G A),     A = b h, G = E / (2 (1 + nu))
    kappa = 10 (1 + nu) / (12 + 11 nu)   (Cowper, rectangle; 5/6 at nu = 0)

    w_S / w_B = 8 (1 + nu) / (5 kappa) (h/L)^2

- w_B: simply supported beam under uniform load. S. P. Timoshenko & J. M. Gere, *Mechanics of
  Materials*, table of beam deflections (appendix; edition and table number **unverified**);
  Roark's *Formulas for Stress and Strain*, beam table (**unverified**, not opened). Formula
  status: not read in a source for this case; the test re-derives it: w_B(x) =
  q x (L^3 - 2 L x^2 + x^3) / (24 E I) satisfies E I w'''' = q, w = 0 and w'' = 0 (zero
  moment) at both ends, and -E I w'' = M(x) = q x (L - x) / 2.
- w_S: **derived here** from Timoshenko beam theory (S. P. Timoshenko, *Phil. Mag.* 41 (1921)
  744-746, theory only; this formula **unverified** against a book). The beam is statically
  determinate, so the shear force follows from equilibrium alone, V(x) = q (L/2 - x) = M'(x),
  and the shear strain adds dw_S/dx = V / (kappa G A); with w_S(0) = 0 this integrates to
  w_S(x) = q x (L - x) / (2 kappa G A), so w_S(L/2) = q L^2 / (8 kappa G A). The test checks
  both steps. Being statically determinate, the bending part is unchanged by shear
  flexibility (the rotation still obeys E I psi' = -M with the same end conditions).
- kappa: G. R. Cowper, "The shear coefficient in Timoshenko's beam theory", *J. Appl. Mech.*
  33 (1966) 335-340, Table 1 (rectangle), the same choice as TS-01. kappa = 5/6 instead would
  raise the reference by 0.012 %.

Numbers: q = 5 N/mm, I = 833.33 mm^4, w_B = 0.595238095 mm, w_S = 0.003642857 mm (0.612 % of
w_B), w_mid = 0.598880952 mm (L/334, so linear theory holds).

Cross-check, not gated: the plane-stress elasticity solution of the uniformly loaded simply
supported beam (Timoshenko & Goodier, *Theory of Elasticity*, 3rd ed., art. 22; formula quoted
from memory, **unverified**) gives w_B [1 + 12/5 (h/L)^2 (4/5 + nu/2)] = 0.598631 mm, 0.042 %
below the Timoshenko-beam value. Both corrections are O((h/L)^2) and well inside the gate.

Why L/h = 20, and the support singularity
-----------------------------------------
A line support carrying a finite force per length is a stress singularity (ADR-006/007), and
in the continuum even the displacement next to the line is logarithmically singular (Flamant):
the support line indents the beam by an amount that grows like ln(1/mesh size). Relative to the
bending deflection this indentation scales like (h/L)^3, so at L/h = 20 it is tiny but visible:
the measured midspan deflection drifts *up* by about 0.06 % per unit of ln(1/size) (roughly 2.5
times the half-space Flamant rate, plausible for a load at the apex of a 90-degree wedge;
inferred, not derived). It would take an absurdly fine mesh for this drift to reach the gate,
so the case gates on deflection only and never on a stress value. The peak stress at the
support lines is checked only for being *flagged* (ADR-007).

Measured through Tessera's public API (gmsh mesher, CalculiX 2.21, October 2026; 4 cores):

    size [mm]  el. through h  nodes   w_mid vs reference  runtime
    5.0             2          2243      +0.086 %         ~0.5 s
    2.5             4          7315      +0.136 %         ~2 s    (gate)
    1.667           6         22630      +0.154 %         ~6 s
    1.25            8         45138      +0.172 %         ~20 s   (not in the sweep)

Without the shear term the errors would be +0.70, +0.75, +0.77, +0.79 %: the 2 % gate does not
by itself detect a missing shear term at this slenderness (TS-08 does that for plates). The
error grows slowly under refinement because of the support indentation above; the convergence
test therefore allows growth inside a quarter of the tolerance, like TS-01.
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

REFERENCE = "simply_supported_beam"

# Geometry (mm), load (MPa), material (MPa)
SPAN = 200.0  # L
WIDTH = 10.0  # b, along y (the model is the half 0 <= y <= b/2)
HEIGHT = 10.0  # h, along z (load direction)
PRESSURE = 0.5  # p on the top face, acting in -z; line load q = p b
YOUNGS_MODULUS = 210_000.0  # steel-generic
POISSON_RATIO = 0.30

GATE_MESH_SIZE = HEIGHT / 4  # mm: 4 quadratic elements through the depth
CONVERGENCE_MESH_SIZES = (HEIGHT / 2, HEIGHT / 4, HEIGHT / 6)  # mm, coarse to fine
CONVERGENCE_SLACK = 0.002  # allowed relative-error increase between refinements (noise)


def line_load(pressure: float = PRESSURE, width: float = WIDTH) -> float:
    """q (N/mm) = p b."""
    return pressure * width


def cowper_kappa(poisson_ratio: float = POISSON_RATIO) -> float:
    """Cowper's (1966) shear coefficient of a rectangle, 10 (1 + nu) / (12 + 11 nu)."""
    return 10.0 * (1.0 + poisson_ratio) / (12.0 + 11.0 * poisson_ratio)


def bending_deflection(x: float, span: float = SPAN) -> float:
    """Euler-Bernoulli deflection (mm) w_B(x) = q x (L^3 - 2 L x^2 + x^3) / (24 E I)."""
    inertia = WIDTH * HEIGHT**3 / 12.0
    return (
        line_load() * x * (span**3 - 2.0 * span * x**2 + x**3) / (24.0 * YOUNGS_MODULUS * inertia)
    )


def shear_deflection(x: float, span: float = SPAN, kappa: float | None = None) -> float:
    """Timoshenko shear deflection (mm) w_S(x) = q x (L - x) / (2 kappa G A)."""
    kappa = cowper_kappa() if kappa is None else kappa
    shear_modulus = YOUNGS_MODULUS / (2.0 * (1.0 + POISSON_RATIO))
    return line_load() * x * (span - x) / (2.0 * kappa * shear_modulus * WIDTH * HEIGHT)


def midspan_deflection() -> float:
    """Reference w_mid (mm) = 5 q L^4 / (384 E I) + q L^2 / (8 kappa G A)."""
    return bending_deflection(SPAN / 2) + shear_deflection(SPAN / 2)


def write_half_beam_step(path: Path) -> Path:
    """Half beam 0 <= x <= L, 0 <= y <= b/2, 0 <= z <= h as STEP (gmsh OpenCASCADE).

    The face y = 0 is split along the neutral axis z = h/2 (a CAD split line), so the midspan
    centroid x = L/2, y = 0, z = h/2 lies on a geometry edge and is a mesh node on every mesh.
    """
    import gmsh

    owns_session = not gmsh.isInitialized()
    if owns_session:
        gmsh.initialize(interruptible=False)
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add(f"benchmark-{path.stem}")
    try:
        occ = gmsh.model.occ
        box = occ.addBox(0.0, 0.0, 0.0, SPAN, WIDTH / 2, HEIGHT)
        # Lower half of the y = 0 face: a rectangle in the xy plane turned into the xz plane.
        lower = occ.addRectangle(0.0, 0.0, 0.0, SPAN, HEIGHT / 2)
        occ.rotate([(2, lower)], 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, math.pi / 2)
        occ.fragment([(3, box)], [(2, lower)])  # imprint: still one solid, y = 0 face split
        occ.synchronize()
        assert len(gmsh.model.getEntities(3)) == 1
        gmsh.write(str(path))
    finally:
        gmsh.model.remove()
        if owns_session:
            gmsh.finalize()
    return path


def _face_with(faces: list[Any], normal: tuple[float, float, float], z: float | None = None) -> int:
    matches = [
        f.index
        for f in faces
        if f.normal is not None
        and np.allclose(f.normal, normal, atol=1e-6)
        and (z is None or math.isclose(f.centroid[2], z, abs_tol=1e-6 * HEIGHT))
    ]
    assert len(matches) == 1, f"expected one face with normal {normal} (z = {z}), got {matches}"
    return matches[0]


def _edge(edges: list[Any], a: int, b: int, length: float) -> tuple[int, int]:
    pair = (min(a, b), max(a, b))
    match = [e for e in edges if e.faces == pair]
    assert len(match) == 1, f"faces {a} and {b} do not share an edge"
    assert match[0].kind == "Line" and match[0].curves == 1
    assert match[0].length == pytest.approx(length)
    return pair


def _solve_beam(step_path: Path, mesh_size: float, workdir: Path) -> Any:
    """Run the half beam through the public API; return the ``StudyResult``."""
    from tessera.core import describe_edges, describe_faces, edit, new_project, solve

    p = new_project("simply-supported-beam")
    edit.add_part(p, step_path, part_id="beam", material_id="steel-generic")
    faces = describe_faces(p, "beam")
    end_pin = _face_with(faces, (-1.0, 0.0, 0.0))  # x = 0
    end_roller = _face_with(faces, (1.0, 0.0, 0.0))  # x = L
    bottom = _face_with(faces, (0.0, 0.0, -1.0))
    top = _face_with(faces, (0.0, 0.0, 1.0))
    sym_lower = _face_with(faces, (0.0, -1.0, 0.0), z=HEIGHT / 4)  # y = 0, below the axis
    sym_upper = _face_with(faces, (0.0, -1.0, 0.0), z=3 * HEIGHT / 4)  # y = 0, above it
    assert len(faces) == 7
    assert next(f for f in faces if f.index == top).area == pytest.approx(SPAN * WIDTH / 2)

    edges = describe_edges(p, "beam")
    pin_line = _edge(edges, end_pin, bottom, WIDTH / 2)
    roller_line = _edge(edges, end_roller, bottom, WIDTH / 2)
    neutral_axis = _edge(edges, sym_lower, sym_upper, SPAN)
    axis_info = next(e for e in edges if e.faces == neutral_axis)
    assert axis_info.point[2] == pytest.approx(HEIGHT / 2)

    edit.add_edge_support(p, "beam", [pin_line], dofs=(1, 2, 3), support_id="pin")
    edit.add_edge_support(p, "beam", [roller_line], dofs=(3,), support_id="roller")
    edit.add_support(p, "beam", [sym_lower, sym_upper], dofs=(2,), support_id="symmetry")
    edit.add_pressure_load(p, "beam", [top], PRESSURE, load_id="q")
    edit.set_mesh_size(p, mesh_size)

    workdir.mkdir(parents=True, exist_ok=True)
    result = solve(p, workdir=workdir)

    nodes = result.mesh.nodes
    u = result.solution.displacements
    tol = 1e-6 * SPAN
    x, y, z = nodes.T
    on_pin = (np.abs(x) < tol) & (np.abs(z) < tol)
    on_roller = (np.abs(x - SPAN) < tol) & (np.abs(z) < tol)
    assert on_pin.sum() >= 2 and on_roller.sum() >= 2, "no nodes on a support line"
    assert np.abs(u[on_pin]).max() < 1e-12, "pinned line moves"
    assert np.abs(u[on_roller, 2]).max() < 1e-12, "roller line moves vertically"
    assert np.abs(u[np.abs(y) < tol, 1]).max() < 1e-12, "symmetry plane y = 0 moves normal"
    # Line supports, not face supports: the ends rotate and the roller slides.
    assert np.abs(u[(np.abs(x) < tol) & (np.abs(z - HEIGHT) < tol), 0]).min() > 1e-4
    assert u[on_roller, 0].min() > 1e-4, "roller line does not slide out as the beam bends"
    return result


def _midspan_deflection(result: Any) -> tuple[float, int]:
    """(-u_z of the midspan centroid node x = L/2, y = 0, z = h/2, number of nodes)."""
    nodes = result.mesh.nodes
    tol = 1e-6 * SPAN
    at = (
        (np.abs(nodes[:, 0] - SPAN / 2) < tol)
        & (np.abs(nodes[:, 1]) < tol)
        & (np.abs(nodes[:, 2] - HEIGHT / 2) < tol)
    )
    assert at.sum() == 1, "no mesh node at the midspan centroid (x = L/2, y = 0, z = h/2)"
    w = -float(result.solution.displacements[at, 2][0])
    assert w > 0.0, "beam did not deflect in the load direction: pressure sign or face"
    return w, result.mesh.n_nodes


@pytest.fixture(scope="module")
def beam_step(tmp_path_factory: pytest.TempPathFactory) -> Path:
    directory = tmp_path_factory.mktemp("simply-supported-beam-geometry")
    return write_half_beam_step(directory / "simply_supported_half_beam.step")


@pytest.fixture(scope="module")
def beam_result(
    beam_step: Path, tmp_path_factory: pytest.TempPathFactory
) -> Callable[[float], Any]:
    """Cached ``mesh_size -> StudyResult`` shared by the gate, the sweep and the flag check."""
    cache: dict[float, Any] = {}

    def run(mesh_size: float) -> Any:
        if mesh_size not in cache:
            workdir = tmp_path_factory.mktemp(f"simply-supported-beam-h{mesh_size:.4g}")
            cache[mesh_size] = _solve_beam(beam_step, mesh_size, workdir)
        return cache[mesh_size]

    return run


def test_closed_form_is_consistent() -> None:
    """Beam equation, simple-support conditions and the shear term from equilibrium."""
    span, q = SPAN, line_load()
    ei = YOUNGS_MODULUS * WIDTH * HEIGHT**3 / 12.0
    # w_B(x) as a polynomial in x.
    w_b = np.polynomial.Polynomial([0.0, span**3, 0.0, -2.0 * span, 1.0]) * (q / (24.0 * ei))
    for x in (0.1 * span, 0.37 * span, 0.5 * span, span):
        assert math.isclose(w_b(x), bending_deflection(x), rel_tol=1e-12, abs_tol=1e-15)
    for x in (0.0, 0.3 * span, span):
        assert math.isclose(ei * w_b.deriv(4)(x), q, rel_tol=1e-12)
    for end in (0.0, span):
        assert math.isclose(w_b(end), 0.0, abs_tol=1e-12)
        assert math.isclose(w_b.deriv(2)(end), 0.0, abs_tol=1e-12)
    moment = np.polynomial.Polynomial([0.0, q * span / 2.0, -q / 2.0])  # M = q x (L - x) / 2
    for x in (0.2 * span, 0.5 * span, 0.9 * span):
        assert math.isclose(-ei * w_b.deriv(2)(x), moment(x), rel_tol=1e-12)
    assert math.isclose(
        bending_deflection(span / 2), 5.0 * q * span**4 / (384.0 * ei), rel_tol=1e-12
    )

    # Shear part: kappa G A dw_S/dx = V = M' = q (L/2 - x), w_S(0) = w_S(L) = 0.
    kappa = cowper_kappa()
    kappa_g_a = kappa * YOUNGS_MODULUS / (2.0 * (1.0 + POISSON_RATIO)) * WIDTH * HEIGHT
    step = 1e-4 * span
    for x in (0.1 * span, 0.3 * span, 0.8 * span):
        slope = (shear_deflection(x + step) - shear_deflection(x - step)) / (2.0 * step)
        assert math.isclose(kappa_g_a * slope, moment.deriv()(x), rel_tol=1e-6)
        assert math.isclose(moment.deriv()(x), q * (span / 2.0 - x), rel_tol=1e-12)
    assert math.isclose(shear_deflection(0.0), 0.0, abs_tol=1e-15)
    assert math.isclose(shear_deflection(span), 0.0, abs_tol=1e-15)
    assert math.isclose(shear_deflection(span / 2), q * span**2 / (8.0 * kappa_g_a), rel_tol=1e-12)

    ratio = shear_deflection(span / 2) / bending_deflection(span / 2)
    expected = 8.0 * (1.0 + POISSON_RATIO) / (5.0 * kappa) * (HEIGHT / span) ** 2
    assert math.isclose(ratio, expected, rel_tol=1e-12)
    assert math.isclose(cowper_kappa(0.0), 5.0 / 6.0, rel_tol=1e-15)
    assert SPAN / HEIGHT >= 20.0  # TS-05: slender, support indentation negligible


def test_reference_matches_formula(load_reference: Callable[[str], dict[str, Any]]) -> None:
    """The stored reference can never drift from the closed-form value."""
    ref = load_reference(REFERENCE)
    assert math.isclose(ref["reference_midspan_deflection_mm"], midspan_deflection(), rel_tol=1e-9)
    assert 0.0 < ref["rel_tolerance"] <= 0.02  # spec gate: within 2 %


def test_simply_supported_beam_midspan_deflection(
    beam_result: Callable[[float], Any],
    load_reference: Callable[[str], dict[str, Any]],
) -> None:
    # Solve first: fail (not skip) if the solve breaks, before the reference check.
    w_mid, _ = _midspan_deflection(beam_result(GATE_MESH_SIZE))

    ref = load_reference(REFERENCE)
    assert w_mid == pytest.approx(ref["reference_midspan_deflection_mm"], rel=ref["rel_tolerance"])


def test_simply_supported_beam_mesh_convergence(
    beam_result: Callable[[float], Any],
    load_reference: Callable[[str], dict[str, Any]],
) -> None:
    # Solve first: fail (not skip) if the solve breaks, before the reference check.
    for size in CONVERGENCE_MESH_SIZES:
        beam_result(size)

    ref = load_reference(REFERENCE)
    expected = ref["reference_midspan_deflection_mm"]
    tolerance = ref["rel_tolerance"]

    rows = []
    for size in CONVERGENCE_MESH_SIZES:
        w_mid, n_nodes = _midspan_deflection(beam_result(size))
        rows.append((size, n_nodes, w_mid, (w_mid - expected) / expected))

    print(
        f"\nSimply supported beam L/h = {SPAN / HEIGHT:g}, p = {PRESSURE:g} MPa: w_mid ref "
        f"{expected:.6e} mm (bending {bending_deflection(SPAN / 2):.6e} + shear "
        f"{shear_deflection(SPAN / 2):.6e})"
    )
    print(f"{'h [mm]':>8} {'nodes':>8} {'w_mid [mm]':>13} {'rel. error':>11}")
    for size, n_nodes, w_mid, error in rows:
        print(f"{size:>8.4g} {n_nodes:>8d} {w_mid:>13.6e} {error:>+11.4%}")

    node_counts = [r[1] for r in rows]
    assert node_counts == sorted(node_counts) and len(set(node_counts)) == len(node_counts), (
        "mesh_size is not refining the mesh"
    )
    errors = [abs(r[3]) for r in rows]
    for coarse, fine in itertools.pairwise(errors):
        # Error must not grow under refinement, unless it already sits well inside the gate (the
        # line supports indent the beam by an amount that grows like ln(1/size)).
        assert fine <= coarse + CONVERGENCE_SLACK or fine < tolerance / 4
    assert errors[-1] < tolerance


def test_simply_supported_beam_line_supports_are_flagged(
    beam_result: Callable[[float], Any],
) -> None:
    """ADR-007: both line supports are flagged, the one-direction symmetry roller is not; the
    raw peak sits at a support line and is flagged, the clear peak is the midspan bending
    stress. No reference value: this checks where and why, not how much."""
    from tessera.model import FlagKind, StressFlag

    result = beam_result(GATE_MESH_SIZE)
    pin = StressFlag(FlagKind.SUPPORT, support="pin")
    roller = StressFlag(FlagKind.SUPPORT, support="roller")
    assert set(result.node_flags.masks) == {pin, roller}

    nodes = result.mesh.nodes
    tol = 1e-6 * SPAN
    bottom_ends = (np.abs(nodes[:, 2]) < tol) & (
        (np.abs(nodes[:, 0]) < tol) | (np.abs(nodes[:, 0] - SPAN) < tol)
    )
    assert result.node_flags.flagged[bottom_ends].all()

    peak = result.peak_von_mises
    assert len(peak.flags) == 1 and peak.flags[0] in (pin, roller)
    assert min(peak.location[0], SPAN - peak.location[0]) <= result.node_flags.distance

    clear = result.clear_peak_von_mises
    assert clear is not None and not clear.flagged
    assert clear.value < peak.value
    assert abs(clear.location[0] - SPAN / 2) <= 0.1 * SPAN, "clear peak is not near midspan"
    assert min(clear.location[2], HEIGHT - clear.location[2]) < tol, "clear peak not on a fibre"
