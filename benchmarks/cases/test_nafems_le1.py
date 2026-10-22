"""NAFEMS LE1 elliptic membrane: tangential stress at the inner ellipse end D (TS-11).

Problem (described in our own words; the NAFEMS figure is not reproduced)
-------------------------------------------------------------------------
A flat membrane of uniform thickness occupies, in its own plane, the region between two
concentric, coaxial ellipses: an inner ellipse with semi-axes 2 m (x) by 1 m (y) and an outer
ellipse with semi-axes 3.25 m (x) by 2.75 m (y). The thickness is 0.1 m. The inner boundary is
free; the outer boundary carries a uniform outward (tensile) normal pressure of 10 MPa. Steel,
E = 210 GPa, nu = 0.3. The problem is symmetric about both axes, so one quarter (x, y >= 0) is
analysed, with roller (symmetry) supports on the two cut planes. Linear static.

The four corner points of the quarter are named as in the published benchmark: D = (2, 0) is
the end of the inner ellipse on the x axis, where the boundary curvature is largest and the
tangential stress peaks.

Units (ADR-002): the benchmark is stated in metres and MPa. Tessera works in mm, N, MPa, so
every length is multiplied by 1000 (2 m -> 2000 mm, 3.25 m -> 3250 mm, 0.1 m -> 100 mm, D =
(2000, 0) mm), E = 210 GPa = 210 000 MPa, and the pressure stays 10 MPa. In linear elasticity
the stress field of geometrically similar bodies under the same pressure does not depend on
the length scale, so the target stress is unchanged in MPa.

Model (mm, N, MPa)
------------------
The thin plate is a 3D solid meshed with quadratic tetrahedra (C3D10, CalculiX). Besides the
two in-plane symmetry planes it is also symmetric about its mid-plane, so half the thickness is
modelled: 0 <= z <= 50 with a roller u_z = 0 on z = 0. That roller also removes the rigid-body
motion along z, which the in-plane rollers leave free; it is exact for a plate of uniform
thickness under in-plane load. So, in total: u_x = 0 on x = 0, u_y = 0 on y = 0, u_z = 0 on
z = 0 (three rollers, each holding only the normal component), and a pressure of -10 MPa on
the outer elliptic face (Tessera's positive pressure pushes into the part, so a negative one
pulls the face outward). Material ``steel-generic`` (E = 210 000 MPa, nu = 0.30), which is
exactly the LE1 material.

Measured quantity
-----------------
sigma_yy (the tangential stress there: the inner boundary at D is normal to x) at the mesh
node nearest D at mid-thickness, (2000, 0, 0). That point is a vertex of the quarter-model
geometry, so a node always sits exactly on it. Tessera has no point probe in ``tessera.post``
yet, so the nodal value is read straight from the solver's nodal stresses (CalculiX's
extrapolated and averaged values, the same field every Tessera peak is taken from) with the
nearest-node rule ``tessera.post.local_edge_length`` uses (a KD-tree query). The peak
conventions of ``tessera.post`` (raw peak, clear peak, per-face peaks, node flags; ADR-007) are
checked separately: the largest von Mises stress must be unflagged and lie on the D edge.

Reference: published target value, not a formula
------------------------------------------------
Target sigma_yy(D) = 92.7 MPa. NAFEMS, *The Standard NAFEMS Benchmarks*, TNSB Rev. 3, October
1990 (NAFEMS R0015), test LE1, as reproduced in the Abaqus Benchmarks Guide (LE1), Altair
OptiStruct verification OS-V 0010 and the COMSOL "elliptic membrane" and Cast3M benchmark
pages. Status per ``benchmarks/CATALOG.md``: the 92.7 MPa target is ``verified`` (seen in
excerpts of the COMSOL, Cast3M and Abaqus pages, which call it analytical); the outer semi-axes
3.25 x 2.75 m are ``unverified`` (they match the 1.25 / 1.75 m wall widths in the COMSOL
figure). No closed-form expression for this value is given in those sources, so this module
has no ``test_reference_matches_formula``. Instead ``test_reference_matches_published_target``
pins the reference file to the published number.

Published fine-mesh quadratic results bracket the target from above: Abaqus CPS8 92.81 MPa
(+0.12 %), per the catalog.

Measured through Tessera's public API (gmsh, CalculiX 2.21, October 2026), sigma_yy at D:

    size [mm]  nodes    mid-thickness       surface z = 50   runtime
    200         2 475   92.029 (-0.72 %)    91.791          ~1 s
    100         9 231   91.999 (-0.76 %)    91.520          ~2.4 s  (sweep)
     70.71     17 247   92.481 (-0.24 %)    91.965          ~5 s    (sweep)
     50        33 931   92.647 (-0.06 %)    92.255          ~11 s   (sweep, gate)
     35.36     65 945   92.744 (+0.05 %)    92.331          ~26 s
     25       174 389   92.785 (+0.09 %)    92.386          ~139 s

The mid-thickness value converges from below to about 92.8 MPa, +0.1 % above the target and
in line with the Abaqus plane-stress 92.81. The surface value sits about 0.4 % lower: in 3D the
stress at a curved free edge is higher at mid-thickness than at the faces (the same effect as
in TS-02), here with full thickness / radius of curvature at D = 100 / 500 = 0.2.

Gate: 50 mm, i.e. one element through the modelled half thickness (two through the full
plate, as the catalog asks) and about ten element edges per 500 mm radius of curvature at D.
The sweep refines by sqrt(2) per step: 100, 70.71, 50 mm.
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
from scipy.spatial import KDTree

pytestmark = [
    pytest.mark.benchmark,
    pytest.mark.solver,
    pytest.mark.skipif(shutil.which("ccx") is None, reason="CalculiX (ccx) not on PATH"),
]

REFERENCE = "nafems_le1"

M_TO_MM = 1000.0  # ADR-002: the benchmark is stated in metres

# Geometry (mm), load (MPa), material (MPa)
INNER_SEMI_AXES = (2.0 * M_TO_MM, 1.0 * M_TO_MM)  # (x, y)
OUTER_SEMI_AXES = (3.25 * M_TO_MM, 2.75 * M_TO_MM)  # (x, y), unverified (see docstring)
THICKNESS = 0.1 * M_TO_MM  # full plate
HALF_THICKNESS = THICKNESS / 2  # modelled, with the mid-plane z = 0 as a symmetry plane
PRESSURE = 10.0  # outward (tensile) on the outer edge
YOUNGS_MODULUS = 210_000.0  # 210 GPa, steel-generic
POISSON_RATIO = 0.30
POINT_D = (INNER_SEMI_AXES[0], 0.0, 0.0)  # (2, 0) m at mid-thickness

NAFEMS_TARGET_SIGMA_YY = 92.7
"""Published LE1 target (MPa), NAFEMS TNSB Rev. 3 (1990); ``verified`` per CATALOG.md."""

GATE_MESH_SIZE = 50.0  # mm: two elements through the full thickness
CONVERGENCE_MESH_SIZES = (100.0, 100.0 / math.sqrt(2.0), 50.0)  # mm, coarse to fine
CONVERGENCE_SLACK = 0.001  # allowed relative-error increase between refinements (noise)


class LE1Result(NamedTuple):
    sigma_yy_mid: float  # at the node on D, mid-thickness (the gated quantity)
    sigma_yy_surface: float  # at the node above D on the free face z = t/2
    n_nodes: int
    result: Any  # the StudyResult


def quarter_ellipse_perimeter(a: float, b: float) -> float:
    """Quarter of an ellipse's perimeter, Ramanujan's second approximation (rel. error < 1e-6
    for b/a >= 0.5)."""
    h = ((a - b) / (a + b)) ** 2
    return math.pi * (a + b) * (1.0 + 3.0 * h / (10.0 + math.sqrt(4.0 - 3.0 * h))) / 4.0


def write_le1_step(path: Path) -> Path:
    """Quarter membrane x, y >= 0 between the ellipses, 0 <= z <= t/2, as STEP (gmsh OCC)."""
    import gmsh

    owns_session = not gmsh.isInitialized()
    if owns_session:
        gmsh.initialize(interruptible=False)
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add(f"benchmark-{path.stem}")
    try:
        occ = gmsh.model.occ
        outer = occ.addDisk(0.0, 0.0, 0.0, *OUTER_SEMI_AXES)
        inner = occ.addDisk(0.0, 0.0, 0.0, *INNER_SEMI_AXES)
        ring, _ = occ.cut([(2, outer)], [(2, inner)])
        span = 2.0 * max(OUTER_SEMI_AXES)
        quadrant = occ.addRectangle(0.0, 0.0, 0.0, span, span)
        quarter, _ = occ.intersect(ring, [(2, quadrant)])
        occ.extrude(quarter, 0.0, 0.0, HALF_THICKNESS)
        occ.synchronize()
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
        f.index for f in faces if f.normal is None and math.isclose(f.area, area, rel_tol=1e-4)
    ]
    assert len(matches) == 1, f"expected one curved face of area {area:.1f}, got {matches}"
    return matches[0]


def _solve_le1(step_path: Path, mesh_size: float, workdir: Path) -> LE1Result:
    """Run the quarter membrane through the public API."""
    from tessera.core import describe_faces, edit, new_project, solve
    from tessera.model import PartInstance

    p = new_project("nafems-le1")
    p.assembly.parts.append(
        PartInstance(
            id="membrane", name="membrane", geometry=step_path, material_id="steel-generic"
        )
    )
    faces = describe_faces(p, "membrane")
    sym_x = _face_with_normal(faces, (-1.0, 0.0, 0.0))  # x = 0 cut plane
    sym_y = _face_with_normal(faces, (0.0, -1.0, 0.0))  # y = 0 cut plane
    sym_z = _face_with_normal(faces, (0.0, 0.0, -1.0))  # mid-plane z = 0
    outer = _curved_face_with_area(
        faces, quarter_ellipse_perimeter(*OUTER_SEMI_AXES) * HALF_THICKNESS
    )
    _curved_face_with_area(faces, quarter_ellipse_perimeter(*INNER_SEMI_AXES) * HALF_THICKNESS)
    assert len(faces) == 6, f"expected 6 faces on the quarter membrane, got {len(faces)}"

    edit.add_support(p, "membrane", [sym_x], dofs=(1,), support_id="sym-x")
    edit.add_support(p, "membrane", [sym_y], dofs=(2,), support_id="sym-y")
    edit.add_support(p, "membrane", [sym_z], dofs=(3,), support_id="sym-z")
    # Positive pressure pushes into the part; the LE1 edge load pulls outward.
    edit.add_pressure_load(p, "membrane", [outer], -PRESSURE, load_id="edge-tension")
    edit.set_mesh_size(p, mesh_size)

    workdir.mkdir(parents=True, exist_ok=True)
    result = solve(p, workdir=workdir)

    nodes = result.mesh.nodes
    u = result.solution.displacements
    tol = 1e-6 * OUTER_SEMI_AXES[0]
    # Rollers, not clamps: normal displacement zero, in-plane motion free.
    on_x0 = np.abs(nodes[:, 0]) < tol
    assert np.abs(u[on_x0, 0]).max() < 1e-12, "x = 0 symmetry plane moves along x"
    assert np.abs(u[on_x0, 1]).max() > 1e-6, "x = 0 face is clamped in y, expected a roller"
    on_y0 = np.abs(nodes[:, 1]) < tol
    assert np.abs(u[on_y0, 1]).max() < 1e-12, "y = 0 symmetry plane moves along y"
    assert np.abs(u[on_y0, 0]).max() > 1e-6, "y = 0 face is clamped in x, expected a roller"
    on_z0 = np.abs(nodes[:, 2]) < tol
    assert np.abs(u[on_z0, 2]).max() < 1e-12, "mid-plane z = 0 moves along z"

    # The load does positive work (twice the strain energy), so with an outward traction the
    # outer face moves outward on average. (Single points need not: the x-axis end moves in,
    # because the hole stretches in y and narrows in x.)
    on_outer = result.mesh.node_sets[f"face:{outer}"]
    a, b = OUTER_SEMI_AXES
    gradient = nodes[on_outer, :2] / np.array([a**2, b**2])
    normal = gradient / np.linalg.norm(gradient, axis=1, keepdims=True)
    outward = float(np.mean(np.sum(u[on_outer, :2] * normal, axis=1)))
    assert outward > 0.0, "outer edge moved inward: pressure sign or face is wrong"

    tree = KDTree(nodes)

    distance, at_d = tree.query(POINT_D)
    assert distance < tol, f"no mesh node at D, nearest is {distance:g} mm away"
    surface_distance, above_d = tree.query((POINT_D[0], POINT_D[1], HALF_THICKNESS))
    assert surface_distance < tol, "no mesh node above D on the free face"
    sigma_yy = result.solution.stresses[:, 1]
    return LE1Result(
        sigma_yy_mid=float(sigma_yy[at_d]),
        sigma_yy_surface=float(sigma_yy[above_d]),
        n_nodes=result.mesh.n_nodes,
        result=result,
    )


@pytest.fixture(scope="module")
def le1_step(tmp_path_factory: pytest.TempPathFactory) -> Path:
    directory = tmp_path_factory.mktemp("nafems-le1-geometry")
    return write_le1_step(directory / "nafems_le1_quarter.step")


@pytest.fixture(scope="module")
def le1_runs(
    le1_step: Path, tmp_path_factory: pytest.TempPathFactory
) -> Callable[[float], LE1Result]:
    """Cached ``mesh_size -> LE1Result`` shared by the gate, the sweep and the peak check."""
    cache: dict[float, LE1Result] = {}

    def run(mesh_size: float) -> LE1Result:
        if mesh_size not in cache:
            workdir = tmp_path_factory.mktemp(f"nafems-le1-h{mesh_size:.4g}")
            cache[mesh_size] = _solve_le1(le1_step, mesh_size, workdir)
        return cache[mesh_size]

    return run


def test_unit_conversion_and_geometry() -> None:
    """ADR-002 conversion and the quarter-perimeter helper used to find the faces."""
    assert POINT_D == (2000.0, 0.0, 0.0)
    assert OUTER_SEMI_AXES == (3250.0, 2750.0)
    assert THICKNESS == 100.0
    assert YOUNGS_MODULUS == 210.0e9 / 1.0e6  # 210 GPa in MPa (N/mm^2)
    # Circle: exact. Ellipse 2 x 1: series value 2.422112 (complete elliptic integral E(k^2=3/4)).
    assert math.isclose(quarter_ellipse_perimeter(1.0, 1.0), math.pi / 2, rel_tol=1e-15)
    assert math.isclose(quarter_ellipse_perimeter(2.0, 1.0), 2.4221121, rel_tol=1e-6)


def test_reference_matches_published_target(
    load_reference: Callable[[str], dict[str, Any]],
) -> None:
    """The stored reference can never drift from the published NAFEMS LE1 target."""
    ref = load_reference(REFERENCE)
    assert math.isclose(ref["reference_sigma_yy_mpa"], NAFEMS_TARGET_SIGMA_YY, rel_tol=1e-9)
    assert 0.0 < ref["rel_tolerance"] <= 0.02  # v0.1 gate: within 2 %


def test_nafems_le1_sigma_yy_at_d(
    le1_runs: Callable[[float], LE1Result],
    load_reference: Callable[[str], dict[str, Any]],
) -> None:
    # Solve first: fail (not skip) if the solve breaks, before the reference check.
    sigma_yy = le1_runs(GATE_MESH_SIZE).sigma_yy_mid

    ref = load_reference(REFERENCE)
    assert sigma_yy == pytest.approx(ref["reference_sigma_yy_mpa"], rel=ref["rel_tolerance"])


def test_nafems_le1_mesh_convergence(
    le1_runs: Callable[[float], LE1Result],
    load_reference: Callable[[str], dict[str, Any]],
) -> None:
    # Solve first: fail (not skip) if the solve breaks, before the reference check.
    results = [le1_runs(size) for size in CONVERGENCE_MESH_SIZES]

    ref = load_reference(REFERENCE)
    expected = ref["reference_sigma_yy_mpa"]
    tolerance = ref["rel_tolerance"]

    print(f"\nNAFEMS LE1 elliptic membrane, sigma_yy at D: reference {expected:g} MPa")
    print(f"{'h [mm]':>8} {'nodes':>8} {'mid [MPa]':>10} {'rel. error':>11} {'surf [MPa]':>11}")
    for size, res in zip(CONVERGENCE_MESH_SIZES, results, strict=True):
        print(
            f"{size:>8.4g} {res.n_nodes:>8d} {res.sigma_yy_mid:>10.4f}"
            f" {res.sigma_yy_mid / expected - 1.0:>+11.4%} {res.sigma_yy_surface:>11.4f}"
        )

    node_counts = [res.n_nodes for res in results]
    assert node_counts == sorted(node_counts) and len(set(node_counts)) == len(node_counts), (
        "mesh_size is not refining the mesh"
    )
    errors = [abs(res.sigma_yy_mid - expected) / expected for res in results]
    for coarse, fine in itertools.pairwise(errors):
        # The error must not grow under refinement, unless it already sits well inside the
        # gate (the converged 3D mid-thickness value is about 0.1 % above the target).
        assert fine <= coarse + CONVERGENCE_SLACK or fine < tolerance / 4
    assert errors[-1] < tolerance


def test_nafems_le1_peak_stress_is_clear_on_the_d_edge(
    le1_runs: Callable[[float], LE1Result],
) -> None:
    """ADR-007 peak conventions: three single-direction rollers and no inside edges, so no node
    is flagged; the raw and the clear von Mises peak are the same node, on the edge through D
    (x = 2000, y = 0), where the inner ellipse meets the y = 0 symmetry plane. No reference
    value: where, not how much."""
    result = le1_runs(GATE_MESH_SIZE).result

    assert not result.node_flags.flagged.any()
    peak = result.peak_von_mises
    assert not peak.flagged
    assert result.clear_peak_von_mises == peak
    x, y, z = peak.location
    tol = 1e-6 * OUTER_SEMI_AXES[0]
    assert x == pytest.approx(POINT_D[0], abs=tol) and abs(y) < tol, "peak not on the D edge"
    assert 0.0 <= z <= HALF_THICKNESS + tol
    assert peak.value == pytest.approx(result.max_von_mises)
