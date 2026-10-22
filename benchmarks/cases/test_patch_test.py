"""Uniform tension patch test on an unstructured mesh (TS-03, exact up to output precision).

Problem (mm, N, MPa per ADR-002)
--------------------------------
Steel cube 0 <= x, y, z <= L with L = 10, E = 210 000 MPa, nu = 0.30 (``steel-generic``).
Rollers on three orthogonal faces: u_x = 0 on x = 0, u_y = 0 on y = 0, u_z = 0 on z = 0
(``Support.dofs``). Uniform traction sigma = 10 MPa in +x on the face x = L, applied as a
``ForceLoad`` of sigma * L^2 = 1000 N, which Tessera spreads uniformly per unit area with
consistent (shape-function weighted) nodal forces. All other faces are free. Linear static,
quadratic tetrahedra (C3D10) in CalculiX.

Exact solution
--------------
The rollers allow free lateral contraction and remove exactly the rigid-body modes, so the
state is homogeneous uniaxial tension:

    sigma_xx = sigma,  sigma_yy = sigma_zz = sigma_xy = sigma_yz = sigma_zx = 0
    u_x = sigma x / E,  u_y = -nu sigma y / E,  u_z = -nu sigma z / E

The displacement field is linear, so it lies in the span of the C3D10 shape functions on any
conforming mesh (straight or curved, regular or not): the finite element solution must
reproduce it exactly at every node, and the strain and stress exactly at every integration
point. Any departure beyond round-off means a wrong load distribution, a wrong boundary
condition, a wrong constitutive matrix or a broken element.

Source: patch test of Irons; the standard test problem set in R. H. MacNeal & R. L. Harder,
"A proposed standard set of problems to test finite element accuracy", Finite Elements in
Analysis and Design 1 (1985) 3-20 (**unverified**: not looked up; the field above is
elementary and is recomputed in ``test_closed_form_is_consistent``).

Numbers: u_x(L) = sigma L / E = 4.761905e-4 mm, u_y(y = L) = u_z(z = L) = -nu sigma L / E =
-1.428571e-4 mm.

The mesh
--------
The public API offers only ``Study.mesh_size``, so the distortion comes from the unstructured
gmsh tetrahedral mesher itself, with sizes that do not divide L (3.0, 1.3, 0.9). There is no
grid: about 1 % of the element edges are parallel to an axis, edge-length ratios within one
element reach 2.5, and dihedral angles run from about 13 to 169 degrees. The case asserts that
the mesh stays unstructured (the ``assert_unstructured`` fixture in ``benchmarks/conftest.py``),
so a future structured mesher cannot turn it into a trivial test. Faces are flat, so the midside
nodes sit on straight edges; curved element edges cannot be produced through the public API on a
box.

Measured quantities
-------------------
At every node of the mesh: sigma_xx against sigma, the five other stress components against 0,
and the three displacement components against the exact linear field. The public API exposes
nodal stresses only (no integration-point values).

Precision and tolerances
------------------------
The catalog proposes 1e-6 relative. That cannot be met through Tessera's result path, for two
measured reasons:

1. CalculiX writes the ``.frd`` results as ``E12.5`` (six significant digits), so every value
   read back carries a rounding error of up to 5e-6 relative (half a unit in the last digit).
2. CalculiX's nodal stresses are extrapolated from the four integration points of each C3D10
   and averaged. For this constant field every integration point gives sigma_xx = 10.000000
   (checked with ``*EL PRINT``), yet every node of every mesh gives 9.99990 in the ``.frd``: a
   systematic -1e-5 relative bias. It is consistent with extrapolation coefficients rounded to
   five decimals (e.g. the exact 1.927051 and 0.309017 rounded to 1.92705 and 0.30902 give a
   row sum of 0.99999); this explanation is **inferred**, the CalculiX source was
   not read.

Proposed gates (fraction of the exact scale, see ``benchmarks/CATALOG.md`` TS-03):
sigma_xx within 2e-5 relative (bias 1e-5 + rounding 5e-6, rounded up); the other stress
components within 1e-7 sigma absolute (measured below 2e-13 sigma, so four orders of headroom
for other sparse solvers, still five orders below any physical error); displacements within
1e-5 of sigma L / E at every node (measured 1.1e-6; rounding alone allows 5e-6).

Measured through Tessera's public API (gmsh mesher, CalculiX 2.21, SPOOLES), all nodes:

    h [mm]  nodes   max|sxx/s - 1|  max|other|/s  max|u - u_ex|/(sL/E)
    3.0      2053      1.0e-05        4.3e-13         1.06e-06     (~1 s, gate)
    1.3      4462      1.0e-05        7.8e-13         1.06e-06     (~2 s)
    0.9     12565      1.0e-05        1.1e-12         1.08e-06     (~8 s)

Sizes from 4.3 down to 2.3 all give the same 2053-node mesh (the mesher saturates), so the
sweep starts at 3.0. The errors stay at the output-precision floor on every mesh.
"""

from __future__ import annotations

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

REFERENCE = "patch_test"

# Geometry (mm), load (MPa), material (MPa)
SIDE = 10.0  # L: cube 0 <= x, y, z <= L
TRACTION = 10.0  # sigma, uniform in +x on x = L
YOUNGS_MODULUS = 210_000.0  # steel-generic
POISSON_RATIO = 0.30

GATE_MESH_SIZE = 3.0  # mm: coarsest distinct mesh (does not divide L)
CONVERGENCE_MESH_SIZES = (3.0, 1.3, 0.9)  # mm, coarse to fine, none divides L


def exact_displacement(
    points: np.ndarray,
    traction: float = TRACTION,
    youngs_modulus: float = YOUNGS_MODULUS,
    poisson_ratio: float = POISSON_RATIO,
) -> np.ndarray:
    """(N, 3) exact displacements (mm): u = (s x, -nu s y, -nu s z) / E."""
    strain = np.array([1.0, -poisson_ratio, -poisson_ratio]) * traction / youngs_modulus
    return points * strain


def axial_end_displacement(
    side: float = SIDE, traction: float = TRACTION, youngs_modulus: float = YOUNGS_MODULUS
) -> float:
    """u_x (mm) on the loaded face x = L: sigma L / E."""
    return traction * side / youngs_modulus


def lateral_end_displacement(
    side: float = SIDE,
    traction: float = TRACTION,
    youngs_modulus: float = YOUNGS_MODULUS,
    poisson_ratio: float = POISSON_RATIO,
) -> float:
    """u_y on y = L (and u_z on z = L), mm: -nu sigma L / E."""
    return -poisson_ratio * traction * side / youngs_modulus


def _face_with_normal(faces: list[Any], normal: tuple[float, float, float]) -> int:
    matches = [
        f.index for f in faces if f.normal is not None and np.allclose(f.normal, normal, atol=1e-6)
    ]
    assert len(matches) == 1, f"expected one face with normal {normal}, got {matches}"
    return matches[0]


def _solve_patch(
    step_path: Path,
    mesh_size: float,
    workdir: Path,
    assert_unstructured: Callable[[np.ndarray, np.ndarray], None],
) -> tuple[float, float, float, int]:
    """Run the cube through the public API.

    Returns, over all nodes: max |sigma_xx / sigma - 1|, max |other stress| / sigma,
    max |u - u_exact| / (sigma L / E), and the number of nodes.
    """
    from tessera.core import describe_faces, new_project, solve
    from tessera.model import FaceSet, ForceLoad, PartInstance, Support

    p = new_project("patch-test")
    p.assembly.parts.append(
        PartInstance(id="cube", name="cube", geometry=step_path, material_id="steel-generic")
    )
    faces = describe_faces(p, "cube")
    loaded = _face_with_normal(faces, (1.0, 0.0, 0.0))
    loaded_info = next(f for f in faces if f.index == loaded)
    assert loaded_info.centroid[0] == pytest.approx(SIDE)
    assert loaded_info.area == pytest.approx(SIDE * SIDE)

    rollers: tuple[tuple[str, tuple[float, float, float], Literal[1, 2, 3]], ...] = (
        ("roller-x", (-1.0, 0.0, 0.0), 1),
        ("roller-y", (0.0, -1.0, 0.0), 2),
        ("roller-z", (0.0, 0.0, -1.0), 3),
    )
    for name, normal, dof in rollers:
        face = _face_with_normal(faces, normal)
        p.study.supports.append(
            Support(id=name, faces=FaceSet(part_id="cube", faces=(face,)), dofs=(dof,))
        )
    p.study.loads.append(
        ForceLoad(
            id="traction",
            faces=FaceSet(part_id="cube", faces=(loaded,)),
            force=(TRACTION * loaded_info.area, 0.0, 0.0),
        )
    )
    p.study.mesh_size = mesh_size

    workdir.mkdir(parents=True, exist_ok=True)
    result = solve(p, workdir=workdir)

    nodes = result.mesh.nodes
    assert_unstructured(nodes, result.mesh.elements)
    u = result.solution.displacements
    s = result.solution.stresses  # sxx, syy, szz, sxy, syz, szx
    tol = 1e-6 * SIDE
    for axis in range(3):
        on_roller = np.abs(nodes[:, axis]) < tol
        assert on_roller.sum() >= 6, f"no mesh nodes on the {'xyz'[axis]} = 0 face"
        assert np.abs(u[on_roller, axis]).max() < 1e-12, f"{'xyz'[axis]} = 0 face moves normal"

    stress_xx_error = float(np.abs(s[:, 0] / TRACTION - 1.0).max())
    other_stress_error = float(np.abs(s[:, 1:]).max() / TRACTION)
    displacement_error = float(
        np.abs(u - exact_displacement(nodes)).max() / axial_end_displacement()
    )
    return stress_xx_error, other_stress_error, displacement_error, result.mesh.n_nodes


@pytest.fixture(scope="module")
def cube_step(box_step: Callable[..., Path]) -> Path:
    return box_step("patch_cube", (0.0, 0.0, 0.0), (SIDE, SIDE, SIDE))


@pytest.fixture(scope="module")
def patch(
    cube_step: Path,
    tmp_path_factory: pytest.TempPathFactory,
    assert_unstructured: Callable[[np.ndarray, np.ndarray], None],
) -> Callable[[float], tuple[float, float, float, int]]:
    """Cached ``mesh_size -> (sxx error, other-stress error, u error, n_nodes)``."""
    cache: dict[float, tuple[float, float, float, int]] = {}

    def run(mesh_size: float) -> tuple[float, float, float, int]:
        if mesh_size not in cache:
            workdir = tmp_path_factory.mktemp(f"patch-test-h{mesh_size:g}")
            cache[mesh_size] = _solve_patch(cube_step, mesh_size, workdir, assert_unstructured)
        return cache[mesh_size]

    return run


def test_closed_form_is_consistent() -> None:
    """The exact field satisfies Hooke's law for uniaxial stress and the roller conditions."""
    e, nu, s = YOUNGS_MODULUS, POISSON_RATIO, TRACTION
    corner = exact_displacement(np.array([[SIDE, SIDE, SIDE]]))[0]
    assert math.isclose(corner[0], axial_end_displacement(), rel_tol=1e-12)
    assert math.isclose(corner[1], lateral_end_displacement(), rel_tol=1e-12)
    assert math.isclose(corner[2], lateral_end_displacement(), rel_tol=1e-12)
    # Strains from the field back through isotropic Hooke's law give (s, 0, 0).
    eps = corner / SIDE
    lam, mu = e * nu / ((1 + nu) * (1 - 2 * nu)), e / (2 * (1 + nu))
    stress = lam * eps.sum() + 2 * mu * eps
    assert np.allclose(stress, [s, 0.0, 0.0], rtol=0.0, atol=1e-12 * s)
    # Rollers: a point on the plane k = 0 has u_k = 0, whatever its other coordinates.
    on_planes = np.full((3, 3), 0.5 * SIDE)
    np.fill_diagonal(on_planes, 0.0)
    assert np.all(np.diag(exact_displacement(on_planes)) == 0.0)


def test_reference_matches_formula(load_reference: Callable[[str], dict[str, Any]]) -> None:
    """The stored reference can never drift from the closed-form value."""
    ref = load_reference(REFERENCE)
    assert math.isclose(ref["reference_stress_xx_mpa"], TRACTION, rel_tol=1e-9)
    assert math.isclose(ref["reference_displacement_x_mm"], axial_end_displacement(), rel_tol=1e-9)
    assert math.isclose(
        ref["reference_displacement_lateral_mm"], lateral_end_displacement(), rel_tol=1e-9
    )
    # A patch test is exact: tolerances only cover output precision, never discretisation.
    assert 0.0 < ref["rel_tolerance_stress_xx"] <= 1e-4
    assert 0.0 < ref["abs_tolerance_other_stress_mpa"] <= 1e-4 * TRACTION
    assert 0.0 < ref["rel_tolerance_displacement"] <= 1e-4


def _check_against_reference(
    errors: tuple[float, float, float, int], ref: dict[str, Any], mesh_size: float
) -> None:
    stress_xx_error, other_stress_error, displacement_error, _ = errors
    sigma = ref["reference_stress_xx_mpa"]
    assert stress_xx_error < ref["rel_tolerance_stress_xx"], (
        f"h = {mesh_size:g}: max |sigma_xx / sigma - 1| = {stress_xx_error:.3e}"
    )
    assert other_stress_error * sigma < ref["abs_tolerance_other_stress_mpa"], (
        f"h = {mesh_size:g}: max |other stress| = {other_stress_error * sigma:.3e} MPa"
    )
    assert displacement_error < ref["rel_tolerance_displacement"], (
        f"h = {mesh_size:g}: max |u - u_exact| / (sigma L / E) = {displacement_error:.3e}"
    )


def test_patch_test(
    patch: Callable[[float], tuple[float, float, float, int]],
    load_reference: Callable[[str], dict[str, Any]],
) -> None:
    # Solve first: fail (not skip) if the solve breaks, before the reference check can skip.
    errors = patch(GATE_MESH_SIZE)

    ref = load_reference(REFERENCE)
    _check_against_reference(errors, ref, GATE_MESH_SIZE)


def test_patch_test_mesh_convergence(
    patch: Callable[[float], tuple[float, float, float, int]],
    load_reference: Callable[[str], dict[str, Any]],
) -> None:
    # Solve first: fail (not skip) if the solve breaks, before the reference check can skip.
    for size in CONVERGENCE_MESH_SIZES:
        patch(size)

    ref = load_reference(REFERENCE)
    print(f"\nUniform tension patch test, sigma = {TRACTION:g} MPa, all nodes")
    print(f"{'h [mm]':>8} {'nodes':>8} {'|sxx/s - 1|':>12} {'|other|/s':>11} {'|du|/(sL/E)':>12}")
    node_counts = []
    for size in CONVERGENCE_MESH_SIZES:
        e_xx, e_other, e_u, n_nodes = patch(size)
        node_counts.append(n_nodes)
        print(f"{size:>8g} {n_nodes:>8d} {e_xx:>12.3e} {e_other:>11.3e} {e_u:>12.3e}")

    assert node_counts == sorted(node_counts) and len(set(node_counts)) == len(node_counts), (
        "mesh_size is not refining the mesh"
    )
    # An exact patch has no discretisation error to converge: every mesh must sit at the
    # output-precision floor, so the gate applies to each refinement, not just the finest.
    for size in CONVERGENCE_MESH_SIZES:
        _check_against_reference(patch(size), ref, size)
