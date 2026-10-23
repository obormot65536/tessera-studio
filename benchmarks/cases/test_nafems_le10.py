"""NAFEMS LE10, thick elliptic plate under pressure, sigma_yy at point D (TS-12).

Problem (metres in the original; here mm, N, MPa per ADR-002)
-------------------------------------------------------------
A quarter of a thick plate whose planform is the NAFEMS LE1 one (TS-11): the region of the
first quadrant (x, y >= 0) between an inner ellipse x^2/2000^2 + y^2/1000^2 = 1 and an outer
ellipse x^2/3250^2 + y^2/2750^2 = 1. The plate is 600 mm thick and symmetric about its
mid-plane, -300 <= z <= 300. Steel, E = 210 000 MPa, nu = 0.30 (``steel-generic``). A uniform
pressure of 1 MPa acts on the upper face z = +300 (``PressureLoad``, positive = into the part,
i.e. -z). Linear static, quadratic tetrahedra (C3D10) in CalculiX.

Corner points of the planform, as named by the benchmark: A = (0, 1000) and B = (0, 2750) on the
y axis, C = (3250, 0) and D = (2000, 0) on the x axis; a primed letter is the same point on the
other face. Boundary conditions, as stated in the catalog (CATALOG.md TS-12):

- face DCD'C' (the y = 0 cut): u_y = 0 (symmetry roller);
- face ABA'B' (the x = 0 cut): u_x = 0 (symmetry roller);
- the outer curved face BCB'C': u_x = u_y = 0;
- u_z = 0 **only on the line EE'**, where the outer curved face meets the mid-plane z = 0.

The plate therefore bends as if simply supported on its outer rim at mid-thickness, with the
rim held in-plane.

How the mid-plane line support is modelled (ADR-006)
----------------------------------------------------
EE' is a line inside the outer curved face, and ADR-006 (decision 7) does not support lines
inside a face: a line support names an edge by the two faces it separates. The geometry is
therefore split at the mid-plane. It is built as two half-thickness slabs (-300..0 and 0..300)
fused with gmsh's ``Geometry.OCCUnionUnify = 0``, so the result is one solid with no internal
face, but every side face is cut in two along z = 0. The outer curved face becomes a lower and an
upper face whose shared edge is exactly EE', and ``edit.add_edge_support(..., [(lower, upper)],
dofs=(3,))`` holds u_z on that edge only. The two outer halves together get the face support
u_x = u_y = 0, and both halves of each cut plane get its roller. The split of the inner face and
the cut planes is harmless: they carry the same condition (or none) on both halves, and the mesh
is conforming anyway; it only forces mesh nodes onto the z = 0 lines.

The alternatives were rejected. A face support u_z = 0 on the whole outer face is a different
problem: it clamps the rim's rotation and moves sigma_yy at D to about -5.24 MPa (measured below,
-2.6 % at 75 mm), the value some later transcriptions of LE10 quote (see "Revision issue"). A
narrow strip face around z = 0 would also hold rotation partly (ADR-006, option B). A node set by
coordinates is not part of the public API (ADR-006, option C).

The case asserts that the support does what LE10 asks: every node on EE' has u_z = 0, and the
rest of the outer face moves in z (the rim rotates about EE'), with u_x = u_y = 0 there.

Measured quantity
-----------------
sigma_yy at D on the upper surface, i.e. the nodal sigma_yy (CalculiX extrapolates from the
integration points and averages over the elements at the node) of the mesh node at
(2000, 0, 300). D is a geometry vertex (where the upper face, the inner curved face and the
y = 0 cut meet), so this node exists on every mesh. It is about 1250 mm from the line support
and not on any concave edge, so ADR-007 flags nothing there; the case checks that through
``StudyResult.node_flags`` (``tessera.post.flag_nodes``) and reports the local element size at
D with ``tessera.post.local_edge_length``.

Reference (published target, not a closed form)
-----------------------------------------------
    sigma_yy(D, top) = -5.38 MPa

NAFEMS, *The Standard NAFEMS Benchmarks*, TNSB Rev. 3, October 1990 (NAFEMS R0015), test LE10.
Status **verified** per CATALOG.md TS-12: seen in web-search excerpts of the Abaqus Benchmarks
Guide (LE10), Altair OptiStruct OS-V 0060 and FeenoX's parametric LE10 (about -5.380 MPa); the
book itself was not opened. Abaqus and Altair quote the magnitude, 5.38; the stress is
compressive because the upper face is the compressed side of the bending. The outer-ellipse
semi-axes 3250 x 2750 come from TS-11 and are **unverified** there (they match the published
wall widths 1.25 m and 1.75 m). The problem is described here in our own words; no NAFEMS text or
figure is copied (licence note, CATALOG.md).

Revision issue: the original TNSB Rev. 3 holds u_z only on EE', and -5.38 MPa belongs to that
support. Some later guides hold u_z over the whole outer face; their converged value is about
-5.24 to -5.25 MPa (**unverified**, ESRD guide and andeplane/fem-lab#183 seen as search excerpts
only). Tessera, run once with that whole-face support in a scratch script (not part of this
case), gives -5.2138 / -5.2280 / -5.2387 MPa at 120 / 100 / 75 mm (-3.09 / -2.83 / -2.63 %),
which agrees with the -5.24 to -5.25 report.

Measured through Tessera's public API (gmsh 4.15, CalculiX 2.21, October 2026, 4 cores). The
sweep rows are from this module; the other rows are from a scratch run of the same solid built
with a smaller construction rectangle, which meshes slightly differently (e.g. 16 182 instead of
16 132 nodes at 120 mm).

    size [mm]  h at D [mm]  nodes    sigma_yy(D) [MPa]  error vs -5.38  runtime
    300           208        2 074       -5.5185          +2.57 %         0.5 s
    250           217        2 669       -5.4033          +0.43 %         0.7 s
    200           178        4 826       -5.5120          +2.45 %         1.6 s   (sweep)
    175           172        5 975       -5.4868          +1.99 %         2.1 s
    150           139        8 781       -5.4438          +1.19 %         3.7 s   (sweep)
    140           112       11 037       -5.3784          -0.03 %         4.0 s
    130           103       13 314       -5.3728          -0.13 %         5.5 s
    120           107       16 132       -5.3721          -0.15 %         7.2 s   (sweep, gate)
    110           103       21 016       -5.3807          +0.01 %        11.4 s
    100            99       27 170       -5.3730          -0.13 %        ~16 s   (sweep)
     75            79       57 755       -5.3765          -0.06 %        82.2 s
     60            61      107 174       -5.3755          -0.08 %       244.0 s

("h at D" is ``local_edge_length`` at D: the median corner-to-corner edge of the elements
there. Error = sigma_yy / (-5.38) - 1, positive = more compressive.)

LE10 is known to converge slowly at D. Here the coarse meshes (local element size at D above
about 140 mm, i.e. fewer than about 4.5 elements through the 600 mm thickness) scatter between
+0.4 % and +2.6 % depending on how the tetrahedra happen to fall at D, so the error is not
monotonic in the global size. Once the local size at D is at or below about 112 mm, every mesh
tried (7 sizes from 140 to 60 mm) gives -5.372 to -5.381 MPa, within -0.15 % / +0.01 % of the
target, and the fine meshes settle at -5.376 MPa (-0.07 %). The remaining difference is about
the rounding of the published three-digit target (+-0.09 %).

Gate: 120 mm (5 element layers through the thickness at D), about 7 s per solve. The sweep is
200, 150, 120, 100 mm (about 25 s in total). It checks that the node count and the element size
at D shrink, that the two finest meshes are inside the gate and within 0.5 % of each other, and
that the finest mesh is closer to the target than the coarsest. It does not require a monotonic
error, for the reason above.
"""

from __future__ import annotations

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

REFERENCE = "nafems_le10"

# Geometry (mm): LE1 planform (TS-11), LE10 thickness. The original is in metres.
INNER_A, INNER_B = 2000.0, 1000.0  # inner ellipse semi-axes along x, y
OUTER_A, OUTER_B = 3250.0, 2750.0  # outer ellipse semi-axes along x, y
THICKNESS = 600.0  # -300 <= z <= 300
PRESSURE = 1.0  # MPa on the upper face z = +300, acting in -z
MATERIAL_ID = "steel-generic"  # E = 210 000 MPa, nu = 0.30, as LE10

POINT_D = (INNER_A, 0.0, THICKNESS / 2)  # upper surface, inner ellipse on the x axis

# NAFEMS TNSB Rev. 3 (1990), LE10: sigma_yy at D on the upper surface (verified per CATALOG.md).
PUBLISHED_SIGMA_YY_D = -5.38  # MPa

GATE_MESH_SIZE = 120.0  # mm
CONVERGENCE_MESH_SIZES = (200.0, 150.0, 120.0, 100.0)  # mm, coarse to fine
SETTLED_CHANGE = 0.005  # the two finest sweep meshes agree to within 0.5 % of the target


class Le10Result(NamedTuple):
    """sigma_yy at D (MPa) and mesh data from one solve."""

    sigma_yy_d: float
    edge_length_d: float  # local element edge length at D (tessera.post), mm
    n_nodes: int


def quarter_ellipse_length(a: float, b: float, n: int = 20_001) -> float:
    """Arc length (mm) of the quarter ellipse with semi-axes a, b (trapezoidal rule)."""
    t = np.linspace(0.0, math.pi / 2, n)
    speed = np.hypot(a * np.sin(t), b * np.cos(t))
    return float(np.sum((speed[1:] + speed[:-1]) / 2.0 * np.diff(t)))


def planform_area() -> float:
    """Area (mm^2) between the quarter ellipses: pi / 4 (a_o b_o - a_i b_i)."""
    return math.pi / 4.0 * (OUTER_A * OUTER_B - INNER_A * INNER_B)


def write_le10_step(path: Path) -> Path:
    """The quarter plate as one solid whose side faces are split at the mid-plane z = 0.

    Two half-thickness slabs are fused without unifying faces (``Geometry.OCCUnionUnify = 0``),
    so the internal face disappears but the outer curved face stays split into a lower and an
    upper face that meet on EE'. The LE1 planform is built here rather than shared with TS-11
    because the case modules are not a package (pytest imports them with
    ``--import-mode=importlib``).
    """
    import gmsh

    owns_session = not gmsh.isInitialized()
    if owns_session:
        gmsh.initialize(interruptible=False)
    gmsh.option.setNumber("General.Terminal", 0)
    unify = gmsh.option.getNumber("Geometry.OCCUnionUnify")
    gmsh.option.setNumber("Geometry.OCCUnionUnify", 0)
    gmsh.model.add(f"benchmark-{path.stem}")
    occ = gmsh.model.occ
    try:
        slabs: list[tuple[int, int]] = []
        for z0 in (-THICKNESS / 2, 0.0):
            outer = occ.addDisk(0.0, 0.0, z0, OUTER_A, OUTER_B)
            inner = occ.addDisk(0.0, 0.0, z0, INNER_A, INNER_B)
            quadrant = occ.addRectangle(0.0, 0.0, z0, 2 * OUTER_A, 2 * OUTER_B)
            ring, _ = occ.cut([(2, outer)], [(2, inner)])
            planform, _ = occ.intersect(ring, [(2, quadrant)])
            slabs += [e for e in occ.extrude(planform, 0.0, 0.0, THICKNESS / 2) if e[0] == 3]
        assert len(slabs) == 2
        occ.fuse([slabs[0]], [slabs[1]])
        occ.synchronize()
        assert len(gmsh.model.getEntities(3)) == 1, "the two slabs did not fuse into one solid"
        gmsh.write(str(path))
    finally:
        gmsh.model.remove()
        gmsh.option.setNumber("Geometry.OCCUnionUnify", unify)
        if owns_session:
            gmsh.finalize()
    return path


def _faces_with_normal(faces: list[Any], normal: tuple[float, float, float]) -> list[int]:
    return [
        f.index for f in faces if f.normal is not None and np.allclose(f.normal, normal, atol=1e-6)
    ]


def _curved_faces_with_area(faces: list[Any], area: float) -> list[int]:
    return [f.index for f in faces if f.normal is None and math.isclose(f.area, area, rel_tol=1e-4)]


def _lower_upper(faces: list[Any], indices: list[int]) -> tuple[int, int]:
    """The two halves of a split side face, ordered (z < 0, z > 0)."""
    assert len(indices) == 2, f"expected a side face split in two, got {indices}"
    by_index = {f.index: f for f in faces}
    lower, upper = sorted(indices, key=lambda i: by_index[i].centroid[2])
    assert by_index[lower].centroid[2] < 0.0 < by_index[upper].centroid[2]
    return lower, upper


def _solve_le10(step_path: Path, mesh_size: float, workdir: Path) -> Le10Result:
    """Run LE10 through the public API and read sigma_yy at D."""
    from tessera.core import describe_edges, describe_faces, edit, new_project, solve
    from tessera.model import PartInstance
    from tessera.post import local_edge_length

    p = new_project("nafems-le10")
    p.assembly.parts.append(
        PartInstance(id="plate", name="plate", geometry=step_path, material_id=MATERIAL_ID)
    )
    faces = describe_faces(p, "plate")
    assert len(faces) == 10, "expected top, bottom and 4 side faces split in two"
    (top,) = _faces_with_normal(faces, (0.0, 0.0, 1.0))
    (bottom,) = _faces_with_normal(faces, (0.0, 0.0, -1.0))
    for index in (top, bottom):
        info = next(f for f in faces if f.index == index)
        assert info.area == pytest.approx(planform_area(), rel=1e-6)
    assert next(f for f in faces if f.index == top).centroid[2] == pytest.approx(THICKNESS / 2)
    cut_y = _faces_with_normal(faces, (0.0, -1.0, 0.0))  # DCD'C', y = 0
    cut_x = _faces_with_normal(faces, (-1.0, 0.0, 0.0))  # ABA'B', x = 0
    outer_area = quarter_ellipse_length(OUTER_A, OUTER_B) * THICKNESS / 2
    outer = _curved_faces_with_area(faces, outer_area)  # BCB'C', split at z = 0
    _lower_upper(faces, cut_y)
    _lower_upper(faces, cut_x)
    outer_lower, outer_upper = _lower_upper(faces, outer)

    line_ee = (min(outer_lower, outer_upper), max(outer_lower, outer_upper))
    (ee_info,) = [e for e in describe_edges(p, "plate") if e.faces == line_ee]
    assert ee_info.kind.lower() == "ellipse"
    assert ee_info.length == pytest.approx(quarter_ellipse_length(OUTER_A, OUTER_B), rel=1e-6)
    assert ee_info.point[2] == pytest.approx(0.0, abs=1e-6)  # the mid-plane line EE'

    edit.add_support(p, "plate", cut_y, dofs=(2,), support_id="sym-y")
    edit.add_support(p, "plate", cut_x, dofs=(1,), support_id="sym-x")
    edit.add_support(p, "plate", outer, dofs=(1, 2), support_id="outer-xy")
    edit.add_edge_support(p, "plate", [line_ee], dofs=(3,), support_id="ee-z")
    edit.add_pressure_load(p, "plate", [top], PRESSURE, load_id="q")
    edit.set_mesh_size(p, mesh_size)

    workdir.mkdir(parents=True, exist_ok=True)
    result = solve(p, workdir=workdir)

    nodes = result.mesh.nodes
    u = result.solution.displacements
    tol = 1e-6 * OUTER_A
    on_outer = np.abs((nodes[:, 0] / OUTER_A) ** 2 + (nodes[:, 1] / OUTER_B) ** 2 - 1.0) < 1e-6
    on_ee = on_outer & (np.abs(nodes[:, 2]) < tol)
    assert on_ee.sum() >= 3, "no mesh nodes on the mid-plane line EE'"
    assert np.abs(u[on_outer, :2]).max() < 1e-12, "outer face moves in x or y"
    assert np.abs(u[on_ee, 2]).max() < 1e-12, "EE' does not hold u_z"
    # Only EE' holds u_z: the rim rotates about it, so the rest of the outer face moves in z.
    # (A whole-face u_z support, the later-transcription variant, would make this zero.)
    rim_uz = np.abs(u[on_outer & ~on_ee, 2]).max()
    assert rim_uz > 1e-2 * np.abs(u[:, 2]).max(), "u_z is held over the outer face, not a line"
    assert np.abs(u[np.abs(nodes[:, 0]) < tol, 0]).max() < 1e-12, "x = 0 cut moves in x"
    assert np.abs(u[np.abs(nodes[:, 1]) < tol, 1]).max() < 1e-12, "y = 0 cut moves in y"

    (d_node,) = np.flatnonzero(np.linalg.norm(nodes - np.array(POINT_D), axis=1) < tol)
    assert u[d_node, 2] < 0.0, "plate did not deflect in -z: pressure sign or face"
    # D is far from the line support and on no concave edge: ADR-007 must not flag it.
    assert result.node_flags.at(int(d_node)) == (), "point D is flagged (ADR-007)"
    return Le10Result(
        sigma_yy_d=float(result.solution.stresses[d_node, 1]),
        edge_length_d=local_edge_length(result.mesh, POINT_D),
        n_nodes=result.mesh.n_nodes,
    )


@pytest.fixture(scope="module")
def le10_step(tmp_path_factory: pytest.TempPathFactory) -> Path:
    directory = tmp_path_factory.mktemp("nafems-le10-geometry")
    return write_le10_step(directory / "nafems_le10_quarter_plate.step")


@pytest.fixture(scope="module")
def nafems_le10(
    le10_step: Path, tmp_path_factory: pytest.TempPathFactory
) -> Callable[[float], Le10Result]:
    """Cached ``mesh_size -> Le10Result`` shared by the gate and the sweep."""
    cache: dict[float, Le10Result] = {}

    def run(mesh_size: float) -> Le10Result:
        if mesh_size not in cache:
            workdir = tmp_path_factory.mktemp(f"nafems-le10-h{mesh_size:.4g}")
            cache[mesh_size] = _solve_le10(le10_step, mesh_size, workdir)
        return cache[mesh_size]

    return run


def test_planform_helpers() -> None:
    """The arc length and area helpers used to identify faces are right."""
    # Circle: quarter perimeter pi r / 2.
    assert math.isclose(quarter_ellipse_length(1000.0, 1000.0), math.pi * 500.0, rel_tol=1e-9)
    # Ramanujan's second approximation for the full ellipse perimeter, error ~1e-10 here.
    a, b = OUTER_A, OUTER_B
    h = ((a - b) / (a + b)) ** 2
    ramanujan = math.pi * (a + b) * (1.0 + 3.0 * h / (10.0 + math.sqrt(4.0 - 3.0 * h)))
    assert math.isclose(4.0 * quarter_ellipse_length(a, b), ramanujan, rel_tol=1e-8)
    assert math.isclose(planform_area(), math.pi / 4.0 * (3250 * 2750 - 2000 * 1000))


def test_reference_matches_published_value(
    load_reference: Callable[[str], dict[str, Any]],
) -> None:
    """The stored reference can never drift from the published NAFEMS target."""
    ref = load_reference(REFERENCE)
    assert math.isclose(ref["reference_sigma_yy_d_mpa"], PUBLISHED_SIGMA_YY_D, rel_tol=1e-9)
    assert 0.0 < ref["rel_tolerance"] <= 0.03  # CATALOG TS-12 proposal: within 3 %


def test_nafems_le10_sigma_yy_at_d(
    nafems_le10: Callable[[float], Le10Result],
    load_reference: Callable[[str], dict[str, Any]],
) -> None:
    # Solve first: fail (not skip) if the solve breaks, before the reference check.
    sigma = nafems_le10(GATE_MESH_SIZE).sigma_yy_d

    ref = load_reference(REFERENCE)
    assert sigma == pytest.approx(ref["reference_sigma_yy_d_mpa"], rel=ref["rel_tolerance"])


def test_nafems_le10_mesh_convergence(
    nafems_le10: Callable[[float], Le10Result],
    load_reference: Callable[[str], dict[str, Any]],
) -> None:
    # Solve first: fail (not skip) if the solve breaks, before the reference check.
    results = [nafems_le10(size) for size in CONVERGENCE_MESH_SIZES]

    ref = load_reference(REFERENCE)
    expected = ref["reference_sigma_yy_d_mpa"]
    tolerance = ref["rel_tolerance"]

    print(
        f"\nNAFEMS LE10, t = {THICKNESS:g} mm, q = {PRESSURE:g} MPa, u_z = 0 on EE' only: "
        f"sigma_yy(D, top) ref {expected:.4f} MPa"
    )
    print(f"{'size [mm]':>10} {'h at D':>8} {'nodes':>8} {'sigma_yy D':>11} {'rel. error':>11}")
    for size, res in zip(CONVERGENCE_MESH_SIZES, results, strict=True):
        print(
            f"{size:>10.4g} {res.edge_length_d:>8.1f} {res.n_nodes:>8d} {res.sigma_yy_d:>11.4f}"
            f" {res.sigma_yy_d / expected - 1.0:>+11.3%}"
        )

    node_counts = [res.n_nodes for res in results]
    assert node_counts == sorted(node_counts) and len(set(node_counts)) == len(node_counts), (
        "mesh_size is not refining the mesh"
    )
    assert results[-1].edge_length_d < results[0].edge_length_d, "the mesh at D did not refine"
    errors = [abs(res.sigma_yy_d - expected) / abs(expected) for res in results]
    # Coarse meshes scatter with the tetrahedra at D (see the module docstring), so the error
    # is not required to fall monotonically; the fine end must be settled and inside the gate.
    assert errors[-1] < errors[0], f"refinement did not approach the target: {errors}"
    assert max(errors[-2:]) < tolerance
    settled = abs(results[-1].sigma_yy_d - results[-2].sigma_yy_d) / abs(expected)
    assert settled < SETTLED_CHANGE, f"sigma_yy at D not settled: {settled:.3%} change"
