"""Thick-walled cylinder under internal pressure, open ends (TS-06, Lamé; gates 2 % and 1 %).

Problem (mm, N, MPa per ADR-002)
--------------------------------
Steel tube along z, bore radius a = 10, outer radius b = 20 (b/a = 2), E = 210 000 MPa,
nu = 0.30 (``steel-generic``), uniform internal pressure p = 10 MPa on the bore
(``PressureLoad``, positive = into the part), outer surface free, ends open (no axial load).
Linear static, quadratic tetrahedra (C3D10) in CalculiX.

Modelled as a quarter-annulus slice 0 <= theta <= 90 deg of length L = 2.5:
- symmetry rollers u_x = 0 on the x = 0 cut plane and u_y = 0 on the y = 0 cut plane;
- roller u_z = 0 on z = 0, the z = L end free.

With sigma_z = 0 the open-end Lamé field is an exact 3D solution: the axial strain
eps_z = -2 nu A / E is uniform, so w = eps_z z is consistent with the u_z = 0 roller at z = 0
and the z = L end stays plane and traction-free. The slice length is therefore arbitrary; it is
kept short for speed.

The bore is curved, and Tessera meshes STEP geometry with curved quadratic elements. The
pressure must act along the local surface normal (a single mean direction for the whole
quarter bore would be grossly wrong).

Measured quantities
-------------------
- sigma_theta(a): on a symmetry plane, sigma_theta is a Cartesian component. It is sigma_yy on
  y = 0 (theta = 0) and sigma_xx on x = 0 (theta = 90 deg). The measure is the mean over all
  bore nodes on the two cut planes (they lie on the straight edges r = a of the slice, so nodes
  always exist there, including the four geometry vertices).
- u_r(a): mean over all bore nodes of the radial displacement u . e_r.

Reference: Lamé
---------------
With A = p a^2 / (b^2 - a^2) and B = p a^2 b^2 / (b^2 - a^2):

    sigma_r = A - B / r^2,   sigma_theta = A + B / r^2,   sigma_z = 0 (open ends)

    sigma_theta(a) = p (a^2 + b^2) / (b^2 - a^2) = 5 p / 3 for b/a = 2

Source: Lamé's solution as given in S. P. Timoshenko, *Strength of Materials* Part II, 3rd ed.,
Van Nostrand 1956, art. 42 (art. number **verified** only via the Ansys VM25 public page
excerpt, CATALOG.md TS-06); formula **verified** on the Wikipedia "Cylinder stress" page
excerpt; 5p/3 is **computed**. Also in Timoshenko & Goodier, *Theory of Elasticity* (section
**unverified**).

Radial displacement (open ends, sigma_z = 0), **computed** here from those stresses with
Hooke's law, u_r = r eps_theta = r (sigma_theta - nu sigma_r) / E:

    u_r(a) = (1 - nu) / E * A a + (1 + nu) / E * B / a

Numbers (p = 10 MPa): A = 3.333333 MPa, B = 1333.333 N, sigma_theta(a) = 16.666667 MPa,
u_r(a) = 9.365079e-4 mm.

Independent check (Tessera's mesher and deck writer, the bore pressure integrated by hand
along the local normal of each curved tri6 into consistent nodal forces, CalculiX 2.21,
October 2026), L = 2.5:

    h [mm]   nodes   sigma_theta(a)   u_r(a)
    2.5        911      +0.17 %       -0.008 %
    1.25      3554      +0.08 %       -0.003 %
    0.625    20761      +0.01 %       -0.0002 %   (~7 s)

The individual bore-node hoop stresses spread by about +-0.5 % around the mean at h = 1.25
(nodal stress extrapolation), which is why the measure averages them.
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

REFERENCE = "lame_cylinder"

# Geometry (mm), load (MPa), material (MPa)
INNER_RADIUS = 10.0  # a
OUTER_RADIUS = 20.0  # b
SLICE_LENGTH = 2.5  # L, along z
PRESSURE = 10.0  # p, internal, on the bore
YOUNGS_MODULUS = 210_000.0  # steel-generic
POISSON_RATIO = 0.30

GATE_MESH_SIZE = 1.25  # mm: 8 elements through the wall
CONVERGENCE_MESH_SIZES = (2.5, 1.25, 0.625)  # mm, coarse to fine
CONVERGENCE_SLACK = 0.002  # allowed relative-error increase between refinements (noise)


def _lame_constants(
    pressure: float, inner_radius: float, outer_radius: float
) -> tuple[float, float]:
    a2, b2 = inner_radius**2, outer_radius**2
    return pressure * a2 / (b2 - a2), pressure * a2 * b2 / (b2 - a2)


def lame_hoop_stress_at_bore(
    pressure: float = PRESSURE,
    inner_radius: float = INNER_RADIUS,
    outer_radius: float = OUTER_RADIUS,
) -> float:
    """sigma_theta(a) (MPa) = p (a^2 + b^2) / (b^2 - a^2)."""
    a2, b2 = inner_radius**2, outer_radius**2
    return pressure * (a2 + b2) / (b2 - a2)


def lame_radial_displacement_at_bore(
    pressure: float = PRESSURE,
    inner_radius: float = INNER_RADIUS,
    outer_radius: float = OUTER_RADIUS,
    youngs_modulus: float = YOUNGS_MODULUS,
    poisson_ratio: float = POISSON_RATIO,
) -> float:
    """u_r(a) (mm), open ends: (1 - nu)/E A a + (1 + nu)/E B / a."""
    a_const, b_const = _lame_constants(pressure, inner_radius, outer_radius)
    return (1.0 - poisson_ratio) / youngs_modulus * a_const * inner_radius + (
        1.0 + poisson_ratio
    ) / youngs_modulus * b_const / inner_radius


def write_quarter_tube_step(path: Path) -> Path:
    """Quarter-annulus slice (0 <= theta <= 90 deg, 0 <= z <= L) as STEP (gmsh OpenCASCADE)."""
    import gmsh

    owns_session = not gmsh.isInitialized()
    if owns_session:
        gmsh.initialize(interruptible=False)
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add(f"benchmark-{path.stem}")
    try:
        occ = gmsh.model.occ
        outer = occ.addCylinder(
            0.0, 0.0, 0.0, 0.0, 0.0, SLICE_LENGTH, OUTER_RADIUS, angle=math.pi / 2
        )
        inner = occ.addCylinder(
            0.0, 0.0, 0.0, 0.0, 0.0, SLICE_LENGTH, INNER_RADIUS, angle=math.pi / 2
        )
        occ.cut([(3, outer)], [(3, inner)])
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
        f.index for f in faces if f.normal is None and math.isclose(f.area, area, rel_tol=1e-6)
    ]
    assert len(matches) == 1, f"expected one curved face of area {area:.4f}, got {matches}"
    return matches[0]


def _solve_lame(step_path: Path, mesh_size: float, workdir: Path) -> tuple[float, float, int]:
    """Run the quarter slice through the public API.

    Returns (mean sigma_theta at r = a on the cut planes, mean u_r at r = a, number of nodes).
    """
    from tessera.core import describe_faces, new_project, solve
    from tessera.model import FaceSet, PartInstance, PressureLoad, Support

    p = new_project("lame-cylinder")
    p.assembly.parts.append(
        PartInstance(id="tube", name="tube", geometry=step_path, material_id="steel-generic")
    )
    faces = describe_faces(p, "tube")
    bore = _curved_face_with_area(faces, math.pi / 2 * INNER_RADIUS * SLICE_LENGTH)
    sym_x = _face_with_normal(faces, (-1.0, 0.0, 0.0))  # x = 0 cut plane
    sym_y = _face_with_normal(faces, (0.0, -1.0, 0.0))  # y = 0 cut plane
    end_0 = _face_with_normal(faces, (0.0, 0.0, -1.0))  # z = 0
    bore_info = next(f for f in faces if f.index == bore)
    assert bore_info.kind.lower() == "cylinder"

    rollers: tuple[tuple[str, int, Literal[1, 2, 3]], ...] = (
        ("sym-x", sym_x, 1),
        ("sym-y", sym_y, 2),
        ("end-z0", end_0, 3),
    )
    for name, face, dof in rollers:
        p.study.supports.append(
            Support(id=name, faces=FaceSet(part_id="tube", faces=(face,)), dofs=(dof,))
        )
    p.study.loads.append(
        PressureLoad(id="bore", faces=FaceSet(part_id="tube", faces=(bore,)), pressure=PRESSURE)
    )
    p.study.mesh_size = mesh_size

    workdir.mkdir(parents=True, exist_ok=True)
    result = solve(p, workdir=workdir)

    nodes = result.mesh.nodes
    u = result.solution.displacements
    s = result.solution.stresses  # sxx, syy, szz, sxy, syz, szx
    tol = 1e-6 * OUTER_RADIUS
    radius = np.hypot(nodes[:, 0], nodes[:, 1])
    on_bore = np.abs(radius - INNER_RADIUS) < 1e-6 * INNER_RADIUS
    on_x0 = np.abs(nodes[:, 0]) < tol
    on_y0 = np.abs(nodes[:, 1]) < tol
    on_z0 = np.abs(nodes[:, 2]) < tol
    on_zl = np.abs(nodes[:, 2] - SLICE_LENGTH) < tol

    u_r = (u[on_bore, 0] * nodes[on_bore, 0] + u[on_bore, 1] * nodes[on_bore, 1]) / radius[on_bore]
    assert u_r.mean() > 0.0, "bore did not expand: pressure missing or acting the wrong way"
    assert np.abs(u[on_x0, 0]).max() < 1e-12, "x = 0 cut plane moves along x"
    assert np.abs(u[on_y0, 1]).max() < 1e-12, "y = 0 cut plane moves along y"
    assert np.abs(u[on_z0, 2]).max() < 1e-12, "z = 0 end moves along z"
    # Rollers, not clamps; open end: the z = L face shortens (eps_z = -2 nu A / E < 0).
    assert np.abs(u[on_x0, 1]).max() > 1e-3 * u_r.mean(), "x = 0 plane is clamped in y"
    assert u[on_zl, 2].mean() < 0.0, "z = L end does not contract axially: end is not open"

    hoop_y0 = s[on_bore & on_y0, 1]  # theta = 0: sigma_theta = sigma_yy
    hoop_x0 = s[on_bore & on_x0, 0]  # theta = 90 deg: sigma_theta = sigma_xx
    assert hoop_y0.size >= 3 and hoop_x0.size >= 3, "no bore nodes found on the cut planes"
    sigma_theta = float(np.concatenate([hoop_y0, hoop_x0]).mean())
    return sigma_theta, float(u_r.mean()), result.mesh.n_nodes


@pytest.fixture(scope="module")
def tube_step(tmp_path_factory: pytest.TempPathFactory) -> Path:
    directory = tmp_path_factory.mktemp("lame-geometry")
    return write_quarter_tube_step(directory / "lame_quarter_tube.step")


@pytest.fixture(scope="module")
def lame(
    tube_step: Path, tmp_path_factory: pytest.TempPathFactory
) -> Callable[[float], tuple[float, float, int]]:
    """Cached ``mesh_size -> (sigma_theta(a), u_r(a), n_nodes)`` shared by gate and sweep."""
    cache: dict[float, tuple[float, float, int]] = {}

    def run(mesh_size: float) -> tuple[float, float, int]:
        if mesh_size not in cache:
            workdir = tmp_path_factory.mktemp(f"lame-h{mesh_size:g}")
            cache[mesh_size] = _solve_lame(tube_step, mesh_size, workdir)
        return cache[mesh_size]

    return run


def test_reference_matches_formula(load_reference: Callable[[str], dict[str, Any]]) -> None:
    """The stored reference can never drift from the closed-form value."""
    ref = load_reference(REFERENCE)
    assert math.isclose(
        ref["reference_hoop_stress_bore_mpa"], lame_hoop_stress_at_bore(), rel_tol=1e-9
    )
    assert math.isclose(
        ref["reference_radial_displacement_bore_mm"],
        lame_radial_displacement_at_bore(),
        rel_tol=1e-9,
    )
    assert 0.0 < ref["rel_tolerance_hoop_stress"] <= 0.02  # spec gate: within 2 %
    assert 0.0 < ref["rel_tolerance_radial_displacement"] <= 0.01  # spec gate: within 1 %


def test_lame_cylinder(
    lame: Callable[[float], tuple[float, float, int]],
    load_reference: Callable[[str], dict[str, Any]],
) -> None:
    # Solve first: fail (not skip) while pressure loads are missing, before the reference check.
    sigma_theta, u_r, _ = lame(GATE_MESH_SIZE)

    ref = load_reference(REFERENCE)
    assert sigma_theta == pytest.approx(
        ref["reference_hoop_stress_bore_mpa"], rel=ref["rel_tolerance_hoop_stress"]
    )
    assert u_r == pytest.approx(
        ref["reference_radial_displacement_bore_mm"], rel=ref["rel_tolerance_radial_displacement"]
    )


def test_lame_cylinder_mesh_convergence(
    lame: Callable[[float], tuple[float, float, int]],
    load_reference: Callable[[str], dict[str, Any]],
) -> None:
    # Solve first: fail (not skip) while pressure loads are missing, before the reference check.
    for size in CONVERGENCE_MESH_SIZES:
        lame(size)

    ref = load_reference(REFERENCE)
    s_ref = ref["reference_hoop_stress_bore_mpa"]
    u_ref = ref["reference_radial_displacement_bore_mm"]

    rows = []
    for size in CONVERGENCE_MESH_SIZES:
        sigma_theta, u_r, n_nodes = lame(size)
        rows.append(
            (
                size,
                n_nodes,
                sigma_theta,
                abs(sigma_theta - s_ref) / s_ref,
                u_r,
                abs(u_r - u_ref) / u_ref,
            )
        )

    print(
        f"\nLamé cylinder b/a = 2, p = {PRESSURE:g} MPa: sigma_theta(a) ref {s_ref:.6f} MPa, "
        f"u_r(a) ref {u_ref:.6e} mm"
    )
    print(
        f"{'h [mm]':>8} {'nodes':>8} {'sigma_th(a)':>12} {'rel. error':>11}"
        f" {'u_r(a) [mm]':>13} {'rel. error':>11}"
    )
    for size, n_nodes, sigma_theta, e_s, u_r, e_u in rows:
        print(
            f"{size:>8g} {n_nodes:>8d} {sigma_theta:>12.6f} {e_s:>11.4%} {u_r:>13.6e} {e_u:>11.4%}"
        )

    node_counts = [r[1] for r in rows]
    assert node_counts == sorted(node_counts) and len(set(node_counts)) == len(node_counts), (
        "mesh_size is not refining the mesh"
    )
    for column, tolerance in (
        (3, ref["rel_tolerance_hoop_stress"]),
        (5, ref["rel_tolerance_radial_displacement"]),
    ):
        errors = [r[column] for r in rows]
        for coarse, fine in itertools.pairwise(errors):
            # Error must not grow under refinement, unless it already sits well inside the gate.
            assert fine <= coarse + CONVERGENCE_SLACK or fine < tolerance / 4
        assert errors[-1] < tolerance
