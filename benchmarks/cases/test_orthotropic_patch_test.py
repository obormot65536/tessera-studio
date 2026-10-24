"""Orthotropic patch tests on an unstructured mesh (TO-01, ADR-009; exact up to output precision).

Problem (mm, N, MPa per ADR-002)
--------------------------------
Cube 0 <= x, y, z <= L with L = 10, made of a project material (``edit.add_material``) with the
catalog's TO constants on the material axes (1 = raster, 2 = across the road, 3 = build):

    E1 = 3000, E2 = 2400, E3 = 1800, G12 = 1000, G13 = 700, G23 = 650 MPa,
    nu12 = 0.35, nu13 = nu23 = 0.30   (major ratios, nu_ij = -eps_j / eps_i under sigma_i)

so nu21 = 0.28, nu31 = 0.18, nu32 = 0.225 and Delta = 0.7427. The material frame comes from the
part's print setup (``edit.set_print_setup``, raster angle 0): build Z gives the identity, build
(1, 0, 0) gives e1 = (0, 0, -1), e2 = (0, 1, 0), e3 = (1, 0, 0) (ADR-009 D2's example, written
out here rather than taken from ``tessera.printsetup``). Every load is a uniform face force of
sigma L^2 = 1000 N (sigma = tau = 10 MPa), spread with consistent nodal forces. Linear static,
C3D10 in CalculiX, problem set up only through ``tessera.core``.

(a) Axial cases (rollers u_x = 0 on x = 0, u_y = 0 on y = 0, u_z = 0 on z = 0, load on the
opposite face): A1, A2, A3 load x, y, z with build Z (material 1, 2, 3); A1-side, A2-side,
A3-side load x, y, z with build (1, 0, 0) (material 3, 2, 1). The state is homogeneous uniaxial
stress, u = eps . x with eps diagonal.

(b) Shear cases, simple shear (ADR-009's replacement for the displacement patch test, which
needs prescribed non-zero displacements Tessera does not have), build Z, other faces free:

    S12: y = 0 fixed; y = L (tau A, 0, 0); x = L (0, tau A, 0); x = 0 (0, -tau A, 0)
    S13: z = 0 fixed; z = L (tau A, 0, 0); x = L (0, 0, tau A); x = 0 (0, 0, -tau A)
    S23: z = 0 fixed; z = L (0, tau A, 0); y = L (0, 0, tau A); y = 0 (0, 0, -tau A)

The exact state is tau_ij = tau with every other stress 0, and u = gamma x_j e_i (gamma =
tau / G_ij, j the fixed face's axis), which is zero on the fixed face, gives exactly the applied
tractions on every loaded face and none on the free faces; forces and moments balance.

Exact solution
--------------
Strains from the orthotropic compliance, eps_mat = S sigma_mat (Voigt 11 22 33 12 23 31,
engineering shear), with sigma_mat = R sigma R^T and eps = R^T eps_mat R (R rows = material axes
in part coordinates). Axial: eps_ii = sigma / E_i, eps_jj = -nu_ij sigma / E_i; shear:
gamma_ij = tau / G_ij. All displacement fields are linear, so the C3D10 solution must reproduce
them at every node and the stresses at every integration point, on any conforming mesh.

Sources: S. G. Lekhnitskii, *Theory of Elasticity of an Anisotropic Body* (Holden-Day 1963 /
Mir 1981), ch. 1 (generalized Hooke's law, engineering constants); R. M. Jones, *Mechanics of
Composite Materials*, 2nd ed. (Taylor & Francis 1999), ch. 2 (orthotropic compliance and the
reciprocity nu_ij / E_i = nu_ji / E_j). Both locations **unverified** (not looked up); every
value is recomputed here from the compliance (``test_closed_form_is_consistent``). Patch test:
Irons; MacNeal & Harder, Finite Elements in Analysis and Design 1 (1985) 3-20 (**unverified**).

Expected values (ADR-009 table, eps in 1e-3, end displacement = eps_load L or gamma L, mm):

    A1       3.333333 -1.166667 -1.000000   0.0333333333
    A2      -1.166667  4.166667 -1.250000   0.0416666667
    A3      -1.000000 -1.250000  5.555556   0.0555555556
    A1-side  5.555556 -1.250000 -1.000000   0.0555555556
    A2-side -1.250000  4.166667 -1.166667   0.0416666667
    A3-side -1.000000 -1.166667  3.333333   0.0333333333
    S12      gamma = 0.01                   0.1
    S13      gamma = 0.0142857143           0.142857142857
    S23      gamma = 0.0153846154           0.153846153846

A2's eps_x = -nu21 sigma / E2 = -nu12 sigma / E1 checks the deck's Poisson convention; the side
cases catch swapped axes, the orientation card order and stresses returned in material axes.

The mesh
--------
As TS-03: the unstructured gmsh tetrahedral mesh at sizes that do not divide L (gate 3.0, sweep
3.0 / 1.3 / 0.9), with ``assert_unstructured`` (``benchmarks/conftest.py``) guarding against a
future structured mesher.

Measured quantities
-------------------
At every node: the loaded stress component against sigma (relative), the five others against 0
(absolute, MPa), and the three displacement components against the exact linear field, scaled
by the case's end displacement. Supported faces are checked to hold exactly.

Precision and tolerances
------------------------
The catalog asks for 1e-6. As measured in TS-03, the ``.frd`` output limits the comparison:
E12.5 rounding (up to 5e-6 relative) and a systematic -1e-5 bias of CalculiX's extrapolated
nodal stresses. ADR-009 therefore proposes TS-03's approved gates: loaded stress 2e-5 relative,
other stresses 1e-6 MPa absolute, displacements 1e-5 of the end displacement. Tightening to 1e-6
needs integration-point output (ADR-008) first.

Measured through Tessera's public API (gmsh, CalculiX 2.21, SPOOLES), all nodes; meshes of
2053 / 4462 / 12565 nodes for h = 3.0 / 1.3 / 0.9 (gate value, then the largest over the sweep):

    case     |s_load/sigma - 1|   max|other| [MPa]      max|du|/u_end        E12.5 bound
    A1       1.0e-05 (all)        3.8e-12 / 9.0e-12     1.51e-06 / 1.52e-06  1.5e-06
    A2       1.0e-05 (all)        4.9e-12 / 1.2e-11     1.21e-06 / 1.24e-06  1.2e-06
    A3       1.0e-05 (all)        4.9e-12 / 1.1e-11     9.13e-07 / 9.32e-07  9.0e-07
    A1-side  1.0e-05 (all)        5.0e-12 / 1.2e-11     9.10e-07 / 9.27e-07  9.0e-07
    A2-side  1.0e-05 (all)        5.4e-12 / 1.1e-11     1.21e-06 / 1.24e-06  1.2e-06
    A3-side  1.0e-05 (all)        4.0e-12 / 1.0e-11     1.51e-06 / 1.52e-06  1.5e-06
    S12      1.0e-05 (all)        1.1e-11 / 3.0e-11     5.21e-07 / 5.28e-07  5.0e-07
    S13      1.0e-05 (all)        1.1e-11 / 3.4e-11     3.45e-06 / 3.51e-06  3.5e-06
    S23      1.0e-05 (all)        1.2e-11 / 3.8e-11     3.26e-06 / 3.28e-06  3.25e-06

The loaded stress reads 9.99990 at every node of every case (the TS-03 extrapolation bias); the
other stresses sit at round-off. The displacement errors are the ``.frd`` rounding itself: an
E12.5 value |u| carries up to 0.5e-5 * 10^floor(log10 |u|), so the largest nodal |u| in each
case gives the bound in the last column (S13 and S23 are higher because their largest u is
above 0.1 mm). The shear fields are therefore reproduced exactly as well: the fixed face and the
face forces are consistent with the linear field. 27 solves (9 cases x 3 meshes) take about
70 s on the development machine; the gate meshes alone about 5 s.
"""

from __future__ import annotations

import math
import shutil
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal, NamedTuple

import numpy as np
import pytest

from tessera.model import OrthotropicDefinition

pytestmark = [
    pytest.mark.benchmark,
    pytest.mark.solver,
    pytest.mark.skipif(shutil.which("ccx") is None, reason="CalculiX (ccx) not on PATH"),
]

REFERENCE = "orthotropic_patch_test"

# Geometry (mm), load (MPa)
SIDE = 10.0  # L: cube 0 <= x, y, z <= L
AREA = SIDE * SIDE
STRESS = 10.0  # sigma = tau

# Catalog TO constants (made up, not user-trusted data), MPa
E1, E2, E3 = 3000.0, 2400.0, 1800.0
G12, G13, G23 = 1000.0, 700.0, 650.0
NU12, NU13, NU23 = 0.35, 0.30, 0.30
MATERIAL = OrthotropicDefinition(
    id="to-catalog", e1=E1, e2=E2, e3=E3, nu12=NU12, nu13=NU13, nu23=NU23,
    g12=G12, g13=G13, g23=G23,
)  # fmt: skip

GATE_MESH_SIZE = 3.0  # mm: coarsest distinct mesh (does not divide L)
CONVERGENCE_MESH_SIZES = (3.0, 1.3, 0.9)  # mm, coarse to fine, none divides L

type Vec3 = tuple[float, float, float]
type Frame = tuple[Vec3, Vec3, Vec3]

# Rows: material axes 1, 2, 3 in part coordinates (ADR-009 D2, raster angle 0).
UPRIGHT: Frame = ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))
ON_SIDE: Frame = ((0.0, 0.0, -1.0), (0.0, 1.0, 0.0), (1.0, 0.0, 0.0))
BUILD_Z: Vec3 = (0.0, 0.0, 1.0)
BUILD_X: Vec3 = (1.0, 0.0, 0.0)

_VOIGT = ((0, 0), (1, 1), (2, 2), (0, 1), (1, 2), (2, 0))  # xx yy zz xy yz zx
_ROLLER_NORMALS: tuple[tuple[Vec3, Literal[1, 2, 3]], ...] = (
    ((-1.0, 0.0, 0.0), 1),
    ((0.0, -1.0, 0.0), 2),
    ((0.0, 0.0, -1.0), 3),
)


def _unit(axis: int, sign: float = 1.0) -> Vec3:
    v = [0.0, 0.0, 0.0]
    v[axis] = sign
    return (v[0], v[1], v[2])


def _force(axis: int, sign: float = 1.0) -> Vec3:
    f = np.array(_unit(axis, sign)) * STRESS * AREA
    return (float(f[0]), float(f[1]), float(f[2]))


@dataclass(frozen=True)
class Case:
    """One TO-01 load case.

    ``component`` (i, j) is the loaded part-axis stress component (sigma_ij = STRESS, the
    others 0). ``fixed_axis`` is None for the axial cases (rollers on the three coordinate
    planes) or j for the shear cases (face x_j = 0 fixed in all dofs). ``loads`` lists
    (outward face normal, total force in N).
    """

    name: str
    build: Vec3
    frame: Frame
    component: tuple[int, int]
    fixed_axis: int | None
    loads: tuple[tuple[Vec3, Vec3], ...]

    @property
    def is_shear(self) -> bool:
        return self.fixed_axis is not None

    @property
    def voigt_index(self) -> int:
        i, j = self.component
        return next(k for k, (a, b) in enumerate(_VOIGT) if {a, b} == {i, j})


def _axial(name: str, build: Vec3, frame: Frame, axis: int) -> Case:
    return Case(name, build, frame, (axis, axis), None, ((_unit(axis), _force(axis)),))


def _shear(name: str, displaced: int, fixed: int) -> Case:
    """Simple shear u = gamma x_fixed e_displaced, face x_fixed = 0 held, build Z."""
    loads = (
        (_unit(fixed), _force(displaced)),
        (_unit(displaced), _force(fixed)),
        (_unit(displaced, -1.0), _force(fixed, -1.0)),
    )
    return Case(name, BUILD_Z, UPRIGHT, (displaced, fixed), fixed, loads)


CASES = (
    _axial("A1", BUILD_Z, UPRIGHT, 0),
    _axial("A2", BUILD_Z, UPRIGHT, 1),
    _axial("A3", BUILD_Z, UPRIGHT, 2),
    _axial("A1-side", BUILD_X, ON_SIDE, 0),
    _axial("A2-side", BUILD_X, ON_SIDE, 1),
    _axial("A3-side", BUILD_X, ON_SIDE, 2),
    _shear("S12", displaced=0, fixed=1),
    _shear("S13", displaced=0, fixed=2),
    _shear("S23", displaced=1, fixed=2),
)
_IDS = [c.name for c in CASES]

# ADR-009's TO-01 tables as printed (eps in 1e-3, 6 decimals; gamma and end u to 10+ digits).
ADR_AXIAL_STRAIN_E3 = {
    "A1": (3.333333, -1.166667, -1.000000),
    "A2": (-1.166667, 4.166667, -1.250000),
    "A3": (-1.000000, -1.250000, 5.555556),
    "A1-side": (5.555556, -1.250000, -1.000000),
    "A2-side": (-1.250000, 4.166667, -1.166667),
    "A3-side": (-1.000000, -1.166667, 3.333333),
}
ADR_SHEAR_STRAIN = {"S12": 0.01, "S13": 0.0142857143, "S23": 0.0153846154}
ADR_END_DISPLACEMENT_MM = {
    "A1": 0.0333333333,
    "A2": 0.0416666667,
    "A3": 0.0555555556,
    "A1-side": 0.0555555556,
    "A2-side": 0.0416666667,
    "A3-side": 0.0333333333,
    "S12": 0.1,
    "S13": 0.142857142857,
    "S23": 0.153846153846,
}


def compliance() -> np.ndarray:
    """6 x 6 compliance (1/MPa), material axes, Voigt 11 22 33 12 23 31, engineering shear."""
    s = np.zeros((6, 6))
    s[0, 0], s[1, 1], s[2, 2] = 1.0 / E1, 1.0 / E2, 1.0 / E3
    s[0, 1] = s[1, 0] = -NU12 / E1
    s[0, 2] = s[2, 0] = -NU13 / E1
    s[1, 2] = s[2, 1] = -NU23 / E2
    s[3, 3], s[4, 4], s[5, 5] = 1.0 / G12, 1.0 / G23, 1.0 / G13
    return s


def _to_voigt(tensor: np.ndarray, shear_factor: float) -> np.ndarray:
    return np.array([tensor[a, b] * (1.0 if a == b else shear_factor) for a, b in _VOIGT])


def _from_voigt(voigt: np.ndarray, shear_factor: float) -> np.ndarray:
    t = np.zeros((3, 3))
    for k, (a, b) in enumerate(_VOIGT):
        t[a, b] = t[b, a] = voigt[k] * (1.0 if a == b else 1.0 / shear_factor)
    return t


def stress_tensor(case: Case) -> np.ndarray:
    """Exact (3, 3) stress in part axes, MPa."""
    i, j = case.component
    sigma = np.zeros((3, 3))
    sigma[i, j] = sigma[j, i] = STRESS
    return sigma


def strain_tensor(case: Case) -> np.ndarray:
    """Exact (3, 3) tensor strain in part axes: eps = R^T (S (R sigma R^T)) R."""
    r = np.array(case.frame)
    stress_mat = _to_voigt(r @ stress_tensor(case) @ r.T, 1.0)
    strain_mat = _from_voigt(compliance() @ stress_mat, 2.0)
    return r.T @ strain_mat @ r


def displacement_gradient(case: Case) -> np.ndarray:
    """H with u = H x: the strain itself (axial) or simple shear (H[i, j] = gamma, else 0)."""
    eps = strain_tensor(case)
    if case.fixed_axis is None:
        return eps
    i, j = case.component
    h = np.zeros((3, 3))
    h[i, j] = 2.0 * eps[i, j]
    return h


def exact_displacement(case: Case, points: np.ndarray) -> np.ndarray:
    """(N, 3) exact displacements (mm)."""
    return points @ displacement_gradient(case).T


def axial_strain(case: Case) -> tuple[float, float, float]:
    d = np.diag(strain_tensor(case))
    return (float(d[0]), float(d[1]), float(d[2]))


def shear_strain(case: Case) -> float:
    """Engineering shear strain gamma = tau / G_ij."""
    i, j = case.component
    return float(2.0 * strain_tensor(case)[i, j])


def end_displacement(case: Case) -> float:
    """Load-direction u (mm) on the loaded face (axial) or the sheared face (shear)."""
    i, j = case.component
    return float(displacement_gradient(case)[i, j] * SIDE)


class Errors(NamedTuple):
    loaded_stress: float  # max |s_load / sigma - 1|
    other_stress_mpa: float  # max |other component|, MPa
    displacement: float  # max |u - u_exact| / end displacement
    n_nodes: int


def _face_with_normal(faces: list[Any], normal: Vec3) -> int:
    matches = [
        f.index for f in faces if f.normal is not None and np.allclose(f.normal, normal, atol=1e-6)
    ]
    assert len(matches) == 1, f"expected one face with normal {normal}, got {matches}"
    return matches[0]


def _solve_case(
    case: Case,
    step_path: Path,
    mesh_size: float,
    workdir: Path,
    assert_unstructured: Callable[[np.ndarray, np.ndarray], None],
) -> Errors:
    from tessera.core import describe_faces, edit, new_project, solve

    project = new_project(f"to-01-{case.name}")
    edit.add_material(project, MATERIAL)
    edit.add_part(project, step_path, part_id="cube", material_id=MATERIAL.id)
    edit.set_print_setup(project, "cube", "fdm", case.build, 0.0)
    faces = describe_faces(project, "cube")

    if case.fixed_axis is None:
        for normal, dof in _ROLLER_NORMALS:
            face = _face_with_normal(faces, normal)
            edit.add_support(project, "cube", [face], dofs=(dof,))
    else:
        face = _face_with_normal(faces, _unit(case.fixed_axis, -1.0))
        edit.add_support(project, "cube", [face], dofs=(1, 2, 3))
    for normal, force in case.loads:
        face = _face_with_normal(faces, normal)
        info = next(f for f in faces if f.index == face)
        assert info.area == pytest.approx(AREA)
        edit.add_force_load(project, "cube", [face], force)
    edit.set_mesh_size(project, mesh_size)

    workdir.mkdir(parents=True, exist_ok=True)
    result = solve(project, workdir=workdir)
    assert result.orthotropic

    nodes = result.mesh.nodes
    assert_unstructured(nodes, result.mesh.elements)
    u = result.solution.displacements
    s = result.solution.stresses  # sxx, syy, szz, sxy, syz, szx
    tol = 1e-6 * SIDE
    if case.fixed_axis is None:
        for axis in range(3):
            on_roller = np.abs(nodes[:, axis]) < tol
            assert on_roller.sum() >= 6, f"no mesh nodes on the {'xyz'[axis]} = 0 face"
            assert np.abs(u[on_roller, axis]).max() < 1e-12, f"{'xyz'[axis]} = 0 face moves"
    else:
        held = np.abs(nodes[:, case.fixed_axis]) < tol
        assert held.sum() >= 6, "no mesh nodes on the fixed face"
        assert np.abs(u[held]).max() < 1e-12, "the fixed face moves"

    k = case.voigt_index
    others = [m for m in range(6) if m != k]
    return Errors(
        loaded_stress=float(np.abs(s[:, k] / STRESS - 1.0).max()),
        other_stress_mpa=float(np.abs(s[:, others]).max()),
        displacement=float(
            np.abs(u - exact_displacement(case, nodes)).max() / abs(end_displacement(case))
        ),
        n_nodes=result.mesh.n_nodes,
    )


@pytest.fixture(scope="module")
def cube_step(box_step: Callable[..., Path]) -> Path:
    return box_step("orthotropic_patch_cube", (0.0, 0.0, 0.0), (SIDE, SIDE, SIDE))


@pytest.fixture(scope="module")
def patch(
    cube_step: Path,
    tmp_path_factory: pytest.TempPathFactory,
    assert_unstructured: Callable[[np.ndarray, np.ndarray], None],
) -> Callable[[Case, float], Errors]:
    """Cached ``(case, mesh_size) -> Errors``."""
    cache: dict[tuple[str, float], Errors] = {}

    def run(case: Case, mesh_size: float) -> Errors:
        key = (case.name, mesh_size)
        if key not in cache:
            workdir = tmp_path_factory.mktemp(f"to-01-{case.name}-h{mesh_size:g}")
            cache[key] = _solve_case(case, cube_step, mesh_size, workdir, assert_unstructured)
        return cache[key]

    return run


def test_closed_form_is_consistent() -> None:
    """The exact fields satisfy the compliance, the supports and the face loads."""
    s = compliance()
    assert np.array_equal(s, s.T)
    nu21, nu31, nu32 = NU12 * E2 / E1, NU13 * E3 / E1, NU23 * E3 / E2
    assert (nu21, nu31, nu32) == pytest.approx((0.28, 0.18, 0.225), rel=1e-12)
    # Reciprocity nu_ij / E_i = nu_ji / E_j, read off the compliance.
    assert math.isclose(-s[1, 0], nu21 / E2, rel_tol=1e-12)
    assert math.isclose(-s[2, 0], nu31 / E3, rel_tol=1e-12)
    assert math.isclose(-s[2, 1], nu32 / E3, rel_tol=1e-12)
    delta = 1 - NU12 * nu21 - NU23 * nu32 - NU13 * nu31 - 2 * nu21 * nu32 * NU13
    assert math.isclose(delta, 0.7427, rel_tol=1e-12)
    assert np.linalg.eigvalsh(s).min() > 0

    for frame, build in ((UPRIGHT, BUILD_Z), (ON_SIDE, BUILD_X)):
        r = np.array(frame)
        assert np.allclose(r @ r.T, np.eye(3), atol=1e-15)
        assert math.isclose(np.linalg.det(r), 1.0, rel_tol=1e-15)
        assert np.array_equal(r[2], build)  # axis 3 is the build direction

    stiffness = np.linalg.inv(s)
    for case in CASES:
        sigma = stress_tensor(case)
        eps = strain_tensor(case)
        h = displacement_gradient(case)
        # The field's strain is eps, and eps goes back through the stiffness to sigma.
        assert np.allclose(0.5 * (h + h.T), eps, rtol=0.0, atol=1e-18)
        r = np.array(case.frame)
        back = _from_voigt(stiffness @ _to_voigt(r @ eps @ r.T, 2.0), 1.0)
        assert np.allclose(r.T @ back @ r, sigma, rtol=0.0, atol=1e-12 * STRESS)
        # Supports: rollers hold u_k = 0 on x_k = 0; the fixed face holds u = 0.
        if case.fixed_axis is None:
            assert np.array_equal(h, np.diag(np.diag(h))), case.name
        else:
            free_columns = [c for c in range(3) if c != case.fixed_axis]
            assert np.all(h[:, free_columns] == 0.0), case.name
        # Face tractions: equal to the applied force on loaded faces, zero on free faces,
        # normal only on roller faces; the fixed face takes the rest.
        applied = {normal: np.array(force) for normal, force in case.loads}
        for axis in range(3):
            for sign in (1.0, -1.0):
                normal = _unit(axis, sign)
                force = sigma @ np.array(normal) * AREA
                if normal in applied:
                    assert np.allclose(force, applied[normal], rtol=0.0, atol=1e-12), case.name
                elif sign < 0 and case.fixed_axis == axis:
                    continue
                elif sign < 0 and case.fixed_axis is None:
                    tangential = np.delete(force, axis)
                    assert np.all(tangential == 0.0), case.name
                else:
                    assert np.all(force == 0.0), (case.name, normal)

    # The ADR's tables, at their printed precision.
    for case in CASES:
        if case.is_shear:
            i, j = case.component
            g = {(0, 1): G12, (0, 2): G13, (1, 2): G23}[(min(i, j), max(i, j))]
            assert math.isclose(shear_strain(case), STRESS / g, rel_tol=1e-12)
            assert shear_strain(case) == pytest.approx(ADR_SHEAR_STRAIN[case.name], abs=1e-10)
        else:
            strain_e3 = np.array(axial_strain(case)) * 1e3
            assert np.allclose(strain_e3, ADR_AXIAL_STRAIN_E3[case.name], rtol=0.0, atol=5e-7)
        assert end_displacement(case) == pytest.approx(
            ADR_END_DISPLACEMENT_MM[case.name], abs=1e-10
        )


def test_reference_matches_formula(load_reference: Callable[[str], dict[str, Any]]) -> None:
    """The stored reference can never drift from the closed-form values."""
    ref = load_reference(REFERENCE)
    assert math.isclose(ref["stress_mpa"], STRESS, rel_tol=1e-9)
    assert set(ref["cases"]) == set(_IDS)
    for case in CASES:
        entry = ref["cases"][case.name]
        if case.is_shear:
            assert math.isclose(entry["shear_strain"], shear_strain(case), rel_tol=1e-9)
        else:
            assert len(entry["strain"]) == 3
            for stored, exact in zip(entry["strain"], axial_strain(case), strict=True):
                assert math.isclose(stored, exact, rel_tol=1e-9), case.name
        assert math.isclose(entry["end_displacement_mm"], end_displacement(case), rel_tol=1e-9)
    # A patch test is exact: tolerances only cover output precision, never discretisation.
    assert 0.0 < ref["rel_tolerance_stress"] <= 1e-4
    assert 0.0 < ref["abs_tolerance_other_stress_mpa"] <= 1e-4 * STRESS
    assert 0.0 < ref["rel_tolerance_displacement"] <= 1e-4


def _check_against_reference(
    case: Case, errors: Errors, ref: dict[str, Any], mesh_size: float
) -> None:
    where = f"{case.name}, h = {mesh_size:g}"
    assert errors.loaded_stress < ref["rel_tolerance_stress"], (
        f"{where}: max |s_load / sigma - 1| = {errors.loaded_stress:.3e}"
    )
    assert errors.other_stress_mpa < ref["abs_tolerance_other_stress_mpa"], (
        f"{where}: max |other stress| = {errors.other_stress_mpa:.3e} MPa"
    )
    assert errors.displacement < ref["rel_tolerance_displacement"], (
        f"{where}: max |u - u_exact| / u_end = {errors.displacement:.3e}"
    )


@pytest.mark.parametrize("case", CASES, ids=_IDS)
def test_orthotropic_patch_test(
    case: Case,
    patch: Callable[[Case, float], Errors],
    load_reference: Callable[[str], dict[str, Any]],
) -> None:
    # Solve first: fail (not skip) if the solve breaks, before the reference check can skip.
    errors = patch(case, GATE_MESH_SIZE)

    ref = load_reference(REFERENCE)
    _check_against_reference(case, errors, ref, GATE_MESH_SIZE)


@pytest.mark.parametrize("case", CASES, ids=_IDS)
def test_orthotropic_patch_test_mesh_convergence(
    case: Case,
    patch: Callable[[Case, float], Errors],
    load_reference: Callable[[str], dict[str, Any]],
) -> None:
    # Solve first: fail (not skip) if the solve breaks, before the reference check can skip.
    for size in CONVERGENCE_MESH_SIZES:
        patch(case, size)

    ref = load_reference(REFERENCE)
    print(f"\nTO-01 {case.name}, sigma = {STRESS:g} MPa, all nodes")
    print(f"{'h [mm]':>8} {'nodes':>8} {'|s/sig - 1|':>12} {'|other| MPa':>12} {'|du|/u_end':>11}")
    node_counts = []
    for size in CONVERGENCE_MESH_SIZES:
        e = patch(case, size)
        node_counts.append(e.n_nodes)
        print(
            f"{size:>8g} {e.n_nodes:>8d} {e.loaded_stress:>12.3e} "
            f"{e.other_stress_mpa:>12.3e} {e.displacement:>11.3e}"
        )

    assert node_counts == sorted(node_counts) and len(set(node_counts)) == len(node_counts), (
        "mesh_size is not refining the mesh"
    )
    # An exact patch has no discretisation error to converge: every mesh must sit at the
    # output-precision floor, so the gate applies to each refinement, not just the finest.
    for size in CONVERGENCE_MESH_SIZES:
        _check_against_reference(case, patch(case, size), ref, size)
