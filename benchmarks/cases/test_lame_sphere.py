"""Thick-walled sphere under internal pressure (TS-07, Lamé; gates 2 % and 1 %).

Problem (mm, N, MPa per ADR-002)
--------------------------------
Hollow steel sphere centred at the origin, inner radius a = 10, outer radius b = 20 (b/a = 2),
E = 210 000 MPa, nu = 0.30 (``steel-generic``), uniform internal pressure p = 10 MPa on the
inner surface (``PressureLoad``, positive = into the part), outer surface free. Linear static,
quadratic tetrahedra (C3D10) in CalculiX.

Modelled as one octant (x, y, z >= 0) with symmetry rollers u_x = 0 on the x = 0 cut plane,
u_y = 0 on y = 0 and u_z = 0 on z = 0. The Lamé field is point-symmetric, so the displacement
normal to each cut plane vanishes and the shear tractions on it are zero: the three rollers are
exact and also remove all rigid-body motion.

The inner surface is doubly curved, and Tessera meshes STEP geometry with curved quadratic
elements. The pressure must act along the local surface normal.

Measured quantities
-------------------
- sigma_t(a): the stress state is equi-biaxial in the tangent plane of the sphere, so every
  direction normal to e_r carries the hoop stress. On a symmetry plane the plane normal is such a
  direction, so sigma_t is a Cartesian component: sigma_xx on x = 0, sigma_yy on y = 0 and
  sigma_zz on z = 0. The measure is the mean over all inner-surface nodes on the three cut planes
  (they lie on the three quarter-circle edges r = a of the octant, so nodes always exist there,
  including the three geometry vertices).
- u_r(a): mean over all inner-surface nodes of the radial displacement u . e_r.

Reference: Lamé
---------------
    sigma_r = -p (b^3/r^3 - 1) / (b^3/a^3 - 1)
    sigma_t = p (b^3/(2 r^3) + 1) / (b^3/a^3 - 1)

    sigma_t(a) = p (b^3 + 2 a^3) / (2 (b^3 - a^3)) = 5 p / 7 for b/a = 2,   sigma_r(a) = -p

    u_r(r) = p a^3 r / (E (b^3 - a^3)) [(1 - 2 nu) + (1 + nu) b^3 / (2 r^3)]

Source: Lamé's solution for a thick spherical shell, S. P. Timoshenko & J. N. Goodier, *Theory of
Elasticity*, 3rd ed., McGraw-Hill 1970, spherical vessel under internal pressure (section
**unverified**); also A. C. Ugural & S. K. Fenster, *Advanced Strength and Applied Elasticity*,
thick-walled spheres (section **unverified**). Formulas written from those standard results and
**not seen** in a source in this session. They are cross-checked here, which needs no book:
sigma_r(a) = -p, sigma_r(b) = 0, and sigma_r + 2 sigma_t is constant (the field is the sum of a
uniform hydrostatic part and a 1/r^3 deviatoric part, which is the general radially symmetric
equilibrium solution). u_r(a) agrees with u_r = r (sigma_t - nu (sigma_r + sigma_t)) / E
(Hooke's law, eps_t = u_r / r), which the test also checks. 5p/7 and the numbers below are
**computed**.

Numbers (p = 10 MPa): sigma_t(a) = 7.142857 MPa, u_r(a) = 3.809524e-4 mm.

Measured through Tessera's public API (gmsh mesher, CalculiX 2.21, October 2026). This is not an
independent check of the pressure load; TS-06 already covers that for a curved face:

    h [mm]   nodes   sigma_t(a)   u_r(a)
    2.5       2561     +1.39 %     -0.084 %    (~1 s)
    1.667     7066     +0.66 %     -0.023 %    (~4 s)
    1.25     15483     +0.48 %     -0.009 %    (~14 s, gate)
    1.0      27785     +0.23 %     -0.005 %    (~34 s, not in the sweep)
    0.833    45103     +0.15 %     -0.003 %    (~89 s, not in the sweep)

h = 5 (1355 nodes, +1.59 %, -0.18 %) is hardly coarser than h = 2.5 on the inner surface (the
curvature limit of the mesher governs), so the sweep starts at 2.5. At h = 2.5 the hoop error
sits at 70 % of the gate, so the gate mesh is 1.25 (8 elements through the wall, as in TS-06).

The individual inner-surface hoop stresses spread from about -0.5 % to +1.6 % at h = 1.25 (nodal
stress extrapolation on a curved surface with a steep 1/r^3 gradient; up to +3.6 % at h = 2.5),
which is why the measure averages them. The mean is biased high and converges from above. The
nodal sigma_r on the inner surface is -9.85 MPa at h = 1.25 (-1.5 % against -p); it is not
gated.
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

REFERENCE = "lame_sphere"

# Geometry (mm), load (MPa), material (MPa)
INNER_RADIUS = 10.0  # a
OUTER_RADIUS = 20.0  # b
PRESSURE = 10.0  # p, internal
YOUNGS_MODULUS = 210_000.0  # steel-generic
POISSON_RATIO = 0.30

GATE_MESH_SIZE = 1.25  # mm: 8 elements through the wall
CONVERGENCE_MESH_SIZES = (2.5, 1.667, 1.25)  # mm, coarse to fine
CONVERGENCE_SLACK = 0.002  # allowed relative-error increase between refinements (noise)


def lame_sphere_radial_stress(
    r: float,
    pressure: float = PRESSURE,
    inner_radius: float = INNER_RADIUS,
    outer_radius: float = OUTER_RADIUS,
) -> float:
    """sigma_r(r) (MPa) = -p (b^3/r^3 - 1) / (b^3/a^3 - 1)."""
    return -pressure * (outer_radius**3 / r**3 - 1.0) / (outer_radius**3 / inner_radius**3 - 1.0)


def lame_sphere_hoop_stress(
    r: float,
    pressure: float = PRESSURE,
    inner_radius: float = INNER_RADIUS,
    outer_radius: float = OUTER_RADIUS,
) -> float:
    """sigma_t(r) (MPa) = p (b^3 / (2 r^3) + 1) / (b^3/a^3 - 1)."""
    return (
        pressure
        * (outer_radius**3 / (2.0 * r**3) + 1.0)
        / (outer_radius**3 / inner_radius**3 - 1.0)
    )


def lame_sphere_radial_displacement(
    r: float,
    pressure: float = PRESSURE,
    inner_radius: float = INNER_RADIUS,
    outer_radius: float = OUTER_RADIUS,
    youngs_modulus: float = YOUNGS_MODULUS,
    poisson_ratio: float = POISSON_RATIO,
) -> float:
    """u_r(r) (mm) = p a^3 r / (E (b^3 - a^3)) [(1 - 2 nu) + (1 + nu) b^3 / (2 r^3)]."""
    a3, b3 = inner_radius**3, outer_radius**3
    return (
        pressure
        * a3
        * r
        / (youngs_modulus * (b3 - a3))
        * ((1.0 - 2.0 * poisson_ratio) + (1.0 + poisson_ratio) * b3 / (2.0 * r**3))
    )


def write_octant_shell_step(path: Path) -> Path:
    """Octant (x, y, z >= 0) of the hollow sphere a <= r <= b as STEP (gmsh OpenCASCADE)."""
    import gmsh

    owns_session = not gmsh.isInitialized()
    if owns_session:
        gmsh.initialize(interruptible=False)
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add(f"benchmark-{path.stem}")
    try:
        occ = gmsh.model.occ
        # angle1/angle2: latitude range (equator to pole, z >= 0); angle3: azimuth (x, y >= 0).
        outer, inner = (
            occ.addSphere(0.0, 0.0, 0.0, radius, angle1=0.0, angle2=math.pi / 2, angle3=math.pi / 2)
            for radius in (OUTER_RADIUS, INNER_RADIUS)
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


def _solve_lame_sphere(
    step_path: Path, mesh_size: float, workdir: Path
) -> tuple[float, float, int]:
    """Run the octant through the public API.

    Returns (mean sigma_t at r = a on the cut planes, mean u_r at r = a, number of nodes).
    """
    from tessera.core import describe_faces, new_project, solve
    from tessera.model import FaceSet, PartInstance, PressureLoad, Support

    p = new_project("lame-sphere")
    p.assembly.parts.append(
        PartInstance(id="shell", name="shell", geometry=step_path, material_id="steel-generic")
    )
    faces = describe_faces(p, "shell")
    inner_face = _curved_face_with_area(faces, math.pi / 2 * INNER_RADIUS**2)  # 4 pi a^2 / 8
    sym_x = _face_with_normal(faces, (-1.0, 0.0, 0.0))  # x = 0 cut plane
    sym_y = _face_with_normal(faces, (0.0, -1.0, 0.0))  # y = 0 cut plane
    sym_z = _face_with_normal(faces, (0.0, 0.0, -1.0))  # z = 0 cut plane
    inner_info = next(f for f in faces if f.index == inner_face)
    assert inner_info.kind.lower() == "sphere"

    rollers: tuple[tuple[str, int, Literal[1, 2, 3]], ...] = (
        ("sym-x", sym_x, 1),
        ("sym-y", sym_y, 2),
        ("sym-z", sym_z, 3),
    )
    for name, face, dof in rollers:
        p.study.supports.append(
            Support(id=name, faces=FaceSet(part_id="shell", faces=(face,)), dofs=(dof,))
        )
    p.study.loads.append(
        PressureLoad(
            id="inner", faces=FaceSet(part_id="shell", faces=(inner_face,)), pressure=PRESSURE
        )
    )
    p.study.mesh_size = mesh_size

    workdir.mkdir(parents=True, exist_ok=True)
    result = solve(p, workdir=workdir)

    nodes = result.mesh.nodes
    u = result.solution.displacements
    s = result.solution.stresses  # sxx, syy, szz, sxy, syz, szx
    tol = 1e-6 * OUTER_RADIUS
    radius = np.linalg.norm(nodes, axis=1)
    on_inner = np.abs(radius - INNER_RADIUS) < 1e-6 * INNER_RADIUS
    on_plane = [np.abs(nodes[:, axis]) < tol for axis in range(3)]  # x = 0, y = 0, z = 0

    u_r = np.einsum("ij,ij->i", u[on_inner], nodes[on_inner]) / radius[on_inner]
    assert u_r.mean() > 0.0, (
        "inner surface did not expand: pressure missing or acting the wrong way"
    )
    for axis, on_cut in enumerate(on_plane):
        assert np.abs(u[on_cut, axis]).max() < 1e-12, f"cut plane {'xyz'[axis]} = 0 moves normal"
        # Rollers, not clamps: the cut plane slides in its own plane.
        in_plane = (axis + 1) % 3
        assert np.abs(u[on_cut, in_plane]).max() > 1e-3 * u_r.mean(), (
            f"cut plane {'xyz'[axis]} = 0 is clamped in {'xyz'[in_plane]}"
        )

    # On the plane normal to axis k, sigma_t = sigma_kk (the plane normal is tangent to the sphere).
    hoop = [s[on_inner & on_cut, axis] for axis, on_cut in enumerate(on_plane)]
    assert all(h.size >= 3 for h in hoop), "no inner-surface nodes found on the cut planes"
    sigma_t = float(np.concatenate(hoop).mean())
    return sigma_t, float(u_r.mean()), result.mesh.n_nodes


@pytest.fixture(scope="module")
def shell_step(tmp_path_factory: pytest.TempPathFactory) -> Path:
    directory = tmp_path_factory.mktemp("lame-sphere-geometry")
    return write_octant_shell_step(directory / "lame_octant_shell.step")


@pytest.fixture(scope="module")
def lame_sphere(
    shell_step: Path, tmp_path_factory: pytest.TempPathFactory
) -> Callable[[float], tuple[float, float, int]]:
    """Cached ``mesh_size -> (sigma_t(a), u_r(a), n_nodes)`` shared by gate and sweep."""
    cache: dict[float, tuple[float, float, int]] = {}

    def run(mesh_size: float) -> tuple[float, float, int]:
        if mesh_size not in cache:
            workdir = tmp_path_factory.mktemp(f"lame-sphere-h{mesh_size:g}")
            cache[mesh_size] = _solve_lame_sphere(shell_step, mesh_size, workdir)
        return cache[mesh_size]

    return run


def test_closed_form_is_consistent() -> None:
    """Boundary conditions, equilibrium invariant and Hooke's law of the Lamé sphere field."""
    a, b, p = INNER_RADIUS, OUTER_RADIUS, PRESSURE
    assert math.isclose(lame_sphere_radial_stress(a), -p, rel_tol=1e-12)
    assert math.isclose(lame_sphere_radial_stress(b), 0.0, abs_tol=1e-12)
    assert math.isclose(lame_sphere_hoop_stress(a), 5.0 * p / 7.0, rel_tol=1e-12)  # b/a = 2
    mean_stress = [lame_sphere_radial_stress(r) + 2.0 * lame_sphere_hoop_stress(r) for r in (a, b)]
    assert math.isclose(mean_stress[0], mean_stress[1], rel_tol=1e-12)
    for r in (a, 0.5 * (a + b), b):
        s_r, s_t = lame_sphere_radial_stress(r), lame_sphere_hoop_stress(r)
        hooke = r * (s_t - POISSON_RATIO * (s_r + s_t)) / YOUNGS_MODULUS
        assert math.isclose(lame_sphere_radial_displacement(r), hooke, rel_tol=1e-12)


def test_reference_matches_formula(load_reference: Callable[[str], dict[str, Any]]) -> None:
    """The stored reference can never drift from the closed-form value."""
    ref = load_reference(REFERENCE)
    assert math.isclose(
        ref["reference_hoop_stress_inner_mpa"],
        lame_sphere_hoop_stress(INNER_RADIUS),
        rel_tol=1e-9,
    )
    assert math.isclose(
        ref["reference_radial_displacement_inner_mm"],
        lame_sphere_radial_displacement(INNER_RADIUS),
        rel_tol=1e-9,
    )
    assert 0.0 < ref["rel_tolerance_hoop_stress"] <= 0.02  # spec gate: within 2 %
    assert 0.0 < ref["rel_tolerance_radial_displacement"] <= 0.01  # spec gate: within 1 %


def test_lame_sphere(
    lame_sphere: Callable[[float], tuple[float, float, int]],
    load_reference: Callable[[str], dict[str, Any]],
) -> None:
    # Solve first: fail (not skip) if the solve breaks, before the reference check.
    sigma_t, u_r, _ = lame_sphere(GATE_MESH_SIZE)

    ref = load_reference(REFERENCE)
    assert sigma_t == pytest.approx(
        ref["reference_hoop_stress_inner_mpa"], rel=ref["rel_tolerance_hoop_stress"]
    )
    assert u_r == pytest.approx(
        ref["reference_radial_displacement_inner_mm"],
        rel=ref["rel_tolerance_radial_displacement"],
    )


def test_lame_sphere_mesh_convergence(
    lame_sphere: Callable[[float], tuple[float, float, int]],
    load_reference: Callable[[str], dict[str, Any]],
) -> None:
    # Solve first: fail (not skip) if the solve breaks, before the reference check.
    for size in CONVERGENCE_MESH_SIZES:
        lame_sphere(size)

    ref = load_reference(REFERENCE)
    s_ref = ref["reference_hoop_stress_inner_mpa"]
    u_ref = ref["reference_radial_displacement_inner_mm"]

    rows = []
    for size in CONVERGENCE_MESH_SIZES:
        sigma_t, u_r, n_nodes = lame_sphere(size)
        rows.append(
            (
                size,
                n_nodes,
                sigma_t,
                abs(sigma_t - s_ref) / s_ref,
                u_r,
                abs(u_r - u_ref) / u_ref,
            )
        )

    print(
        f"\nLamé sphere b/a = 2, p = {PRESSURE:g} MPa: sigma_t(a) ref {s_ref:.6f} MPa, "
        f"u_r(a) ref {u_ref:.6e} mm"
    )
    print(
        f"{'h [mm]':>8} {'nodes':>8} {'sigma_t(a)':>12} {'rel. error':>11}"
        f" {'u_r(a) [mm]':>13} {'rel. error':>11}"
    )
    for size, n_nodes, sigma_t, e_s, u_r, e_u in rows:
        print(f"{size:>8g} {n_nodes:>8d} {sigma_t:>12.6f} {e_s:>11.4%} {u_r:>13.6e} {e_u:>11.4%}")

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
