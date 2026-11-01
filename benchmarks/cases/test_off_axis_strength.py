"""Off-axis strength ratio, Tsai-Wu (TO-06, ADR-010 PR-5; exact up to output precision).

Problem (mm, N, MPa per ADR-002)
--------------------------------
Cube 0 <= x, y, z <= L with L = 10 under the uniform uniaxial stress sigma_xx = sigma = 10 MPa,
every other component 0. Elastic constants: the catalog's TO constants (as TO-01; they do not
enter R, since the stress is statically determined). Strengths on the material axes (1 =
raster, 2 = across the road in the layer, 3 = build), **made up for the benchmark, not
user-trusted data**:

    T1 = 40, C1 = 50, T2 = 30, C2 = 45, T3 = 15, C3 = 40, S12 = 20, S13 = S23 = 12 MPa

(a) Build Z, raster angle theta = 0, 30, 45, 60, 90 deg. The material frame comes from
``tessera.printsetup.material_frame((0, 0, 1), theta)``: e1 = (cos, sin, 0), e2 = (-sin, cos,
0), e3 = z (axis 1 turned counter-clockwise about the build direction). The off-axis material
couples sigma_xx to gamma_xy, so rollers on x = 0 and y = 0 would hold the shear and spoil the
uniform state. Instead: equal and opposite face forces sigma L^2 = 1000 N on x = L and x = 0
(consistent nodal forces, ``tessera.physics.uniform_face_load``) and a statically determinate
3-2-1 support on three corner nodes, origin (u = 0), (L, 0, 0) (u_y = u_z = 0) and (0, L, 0)
(u_z = 0). The loads balance, so the supports carry no reaction and the exact field is
u = eps . x plus the rigid rotation the supports fix: uniform stress at every point. Core has
face and edge supports but no point supports, so this case is built at the lower layers
(``tessera.mesh.mesh_step`` -> ``SolveJob`` -> ``CalculixAdapter`` -> ``tessera.post.summarize``),
as ADR-008's pure-bending case does.

(b) Printed on its side, build (1, 0, 0), raster 0: x is material axis 3, so the load is across
the layers and the material frame is aligned with the part axes (no coupling). This case runs
end to end through ``tessera.core.solve`` with rollers on x = 0, y = 0, z = 0 and 1000 N on
x = L (exact, as in TO-01's A3-side).

Exact solution
--------------
Material-axis stresses (Voigt 11 22 33 12 23 31), sigma' = R sigma R^T with R's rows e1 e2 e3:

    (a) (sigma cos^2, sigma sin^2, 0, -sigma sin cos, 0, 0)
    (b) (0, 0, sigma, 0, 0, 0)

Tsai-Wu (ADR-010 D2): F_i s_i + F_ij s_i s_j = 1 at failure, with F_i = 1/T_i - 1/C_i,
F_ii = 1/(T_i C_i), F_44 = 1/S12^2, F_55 = 1/S23^2, F_66 = 1/S13^2 and the interaction terms
F_ij = F*_ij sqrt(F_ii F_jj), F*_ij = -1/2 for the three normal pairs. Scaling every load by R
scales the stress by R, so the strength ratio R is the positive root of

    a R^2 + b R - 1 = 0,   a = F_ij s_i s_j,   b = F_i s_i,   R = (-b + sqrt(b^2 + 4 a)) / (2 a).

(b) reduces to F33 s^2 R^2 + F3 s R = 1, whose root is exactly R = T3 / sigma = 1.5.

Sign of tau_12. A raster at -theta flips only the sign of tau_12. In Tsai-Wu tau_12 appears
only in the quadratic shear term F_66 tau_12^2 (orthotropic: no linear shear term, no
normal-shear interaction), so R is even in tau_12 and the same for +theta and -theta. A
reversed raster-angle convention could therefore not be caught by R; TO-02 (strains) is the
case that checks the sign of the rotation. ``test_closed_form_is_consistent`` asserts both the
printsetup sign used here and the evenness.

Failure mode (ADR-010 D2): the term with the largest maximum-stress ratio |s_k| / strength_k.

    theta   s11/T1   s22/T2   |t12|/S12    mode
    0       0.250    0        0            along the raster, tension
    30      0.1875   0.0833   0.2165       in-layer shear
    45      0.125    0.1667   0.250        in-layer shear
    60      0.0625   0.250    0.2165       across the roads, tension
    90      0        0.333    0            across the roads, tension
    side    s33/T3 = 0.667                 across the layers, tension

Sources: S. W. Tsai and E. M. Wu, "A general theory of strength for anisotropic materials",
J. Composite Materials 5 (1971) 58-80 (**unverified**, not looked up); S. W. Tsai and H. T.
Hahn, *Introduction to Composite Materials* (Technomic 1980) for F*_ij = -1/2
(**unverified**). The stress transformation: R. M. Jones, *Mechanics of Composite Materials*,
2nd ed. (1999), ch. 2 (**unverified**). Every value is recomputed here independently of
``tessera.post`` (``test_closed_form_is_consistent``).

Expected values (ADR-010 table, ``computed``): R = 4.000000, 3.380006, 3.096069, 2.982983,
3.000000 for theta = 0, 30, 45, 60, 90 deg; R = 1.5 on its side.

The mesh
--------
The unstructured gmsh C3D10 mesh of the cube at sizes that do not divide L (gate 3.0, sweep
3.0 / 1.3 / 0.9, as TO-01), guarded by ``assert_unstructured``. The field is exact, so there is
no discretisation error to converge: every refinement must sit at the precision floor.

Measured quantities
-------------------
``StudyResult.strength``: ``lowest.ratio`` (the minimum over every node of the part),
``clear_lowest.ratio`` where it exists ((a) flags the corner supports, so it exists there;
(b) has only rollers, so nothing is flagged and it equals ``lowest``), ``lowest.mode``, and the
lowest ratio of every face. Also every element-nodal ratio (each element's own stress at its
ten nodes, ``tessera.post.element_nodal_stresses`` + ``strength_ratio``), so the maximum
deviation is taken over the whole part, not only at the minimum.

Precision and tolerance
-----------------------
ADR-010 proposes the element-nodal floor approved for ADR-008 (``patch_test_points.toml``,
``rel_tolerance_element_nodal_stress = 1.5e-6``): ``job.dat`` prints E13.6 and the tet4 fit
weights a corner node by at most 2.854, so each element-nodal stress is within 1.43e-6 of the
element's largest stress. R(k s) = R(s) / k, so a common relative error passes into R one to
one; the components here are within a factor 1.3 of sigma, so R stays near the same floor.

Measured (CalculiX 2.21, SPOOLES), max |R / R_exact - 1| over every element-nodal ratio and
over ``lowest``, ``clear_lowest`` and the face minima; meshes of 2053 / 4462 / 12565 nodes for
h = 3.0 / 1.3 / 0.9 (gate value, then the largest over the sweep):

    case        R exact          lowest R (every mesh)   max |R / R_exact - 1|   mode
    raster-0    4.0000000000     4.0000000000            4.5e-13 / 8.2e-13       as table
    raster-30   3.3800055152     3.3800055241            2.6e-09 / 2.6e-09       as table
    raster-45   3.0960693596     3.0960693596            9.3e-14 / 1.3e-12       as table
    raster-60   2.9829826706     2.9829826769            2.1e-09 / 2.1e-09       as table
    raster-90   3.0000000000     3.0000000000            3.6e-13 / 1.7e-12       as table
    side        1.5000000000     1.5000000000            1.1e-13 / 2.1e-13       as table

The 30 and 60 deg cases carry the E13.6 rounding of tau_12 = 4.330127 (exact 4.3301270189),
the same at every node and every mesh; the stresses of the others print exactly, so they sit
at round-off. Every case is about 600 times inside the 1.5e-6 floor. 18 solves take about
40 s on the development machine.

The gate compares against ``references/off_axis_strength.toml`` (proposed in the PR, created
by the maintainer); until it exists the solver cases are skipped after solving.
"""

from __future__ import annotations

import dataclasses
import math
import shutil
from collections.abc import Callable
from pathlib import Path
from typing import Any, Literal, NamedTuple

import numpy as np
import pytest
from numpy.typing import NDArray

from tessera.model import FailureMode, OrthotropicDefinition, OrthotropicStrengthDefinition

pytestmark = [
    pytest.mark.benchmark,
    pytest.mark.solver,
    pytest.mark.skipif(shutil.which("ccx") is None, reason="CalculiX (ccx) not on PATH"),
]

REFERENCE = "off_axis_strength"

# Geometry (mm), load (MPa)
SIDE = 10.0
AREA = SIDE * SIDE
STRESS = 10.0

# Catalog TO constants (made up, not user-trusted data), MPa; R does not depend on them.
E1, E2, E3 = 3000.0, 2400.0, 1800.0
G12, G13, G23 = 1000.0, 700.0, 650.0
NU12, NU13, NU23 = 0.35, 0.30, 0.30

# TO-06 strengths (made up, not user-trusted data), MPa.
T1, C1, T2, C2, T3, C3 = 40.0, 50.0, 30.0, 45.0, 15.0, 40.0
S12, S13, S23 = 20.0, 12.0, 12.0
F_STAR = -0.5
"""F*_ij for the normal-stress pairs (ADR-010 D2), restated here rather than imported."""

STRENGTH = OrthotropicStrengthDefinition(
    t1=T1, c1=C1, t2=T2, c2=C2, t3=T3, c3=C3, s12=S12, s13=S13, s23=S23
)
MATERIAL = OrthotropicDefinition(
    id="to-06", e1=E1, e2=E2, e3=E3, nu12=NU12, nu13=NU13, nu23=NU23,
    g12=G12, g13=G13, g23=G23, strength=STRENGTH,
)  # fmt: skip

GATE_MESH_SIZE = 3.0  # mm: does not divide L
CONVERGENCE_MESH_SIZES = (3.0, 1.3, 0.9)  # mm, coarse to fine, none divides L

type Vec3 = tuple[float, float, float]

BUILD_Z: Vec3 = (0.0, 0.0, 1.0)
BUILD_X: Vec3 = (1.0, 0.0, 0.0)


class Case(NamedTuple):
    name: str
    build: Vec3
    raster_deg: float
    material_stress: tuple[float, float, float, float, float, float]  # 11 22 33 12 23 31
    mode: FailureMode


def _raster(theta: float, mode: FailureMode) -> Case:
    c, s = math.cos(math.radians(theta)), math.sin(math.radians(theta))
    stress = (STRESS * c * c, STRESS * s * s, 0.0, -STRESS * s * c, 0.0, 0.0)
    return Case(f"raster-{theta:g}", BUILD_Z, theta, stress, mode)


RASTER_CASES = (
    _raster(0.0, FailureMode.ALONG_RASTER_TENSION),
    _raster(30.0, FailureMode.IN_LAYER_SHEAR),
    _raster(45.0, FailureMode.IN_LAYER_SHEAR),
    _raster(60.0, FailureMode.ACROSS_ROADS_TENSION),
    _raster(90.0, FailureMode.ACROSS_ROADS_TENSION),
)
SIDE_CASE = Case(
    "side", BUILD_X, 0.0, (0.0, 0.0, STRESS, 0.0, 0.0, 0.0), FailureMode.ACROSS_LAYERS_TENSION
)
CASES = (*RASTER_CASES, SIDE_CASE)
_RASTER_IDS = [c.name for c in RASTER_CASES]
_IDS = [c.name for c in CASES]

ADR_RATIOS = {
    "raster-0": 4.000000,
    "raster-30": 3.380006,
    "raster-45": 3.096069,
    "raster-60": 2.982983,
    "raster-90": 3.000000,
    "side": 1.5,
}
"""ADR-010's TO-06 table as printed (6 decimals)."""


# --- Closed form (independent of tessera.post) --------------------------------------------


def tsai_wu_coefficients() -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """(F_i (6,), F_ij (6, 6)) in Voigt 11 22 33 12 23 31."""
    t = np.array([T1, T2, T3])
    c = np.array([C1, C2, C3])
    linear = np.zeros(6)
    linear[:3] = 1.0 / t - 1.0 / c
    quad = np.zeros((6, 6))
    quad[:3, :3] = np.diag(1.0 / (t * c))
    for i, j in ((0, 1), (0, 2), (1, 2)):
        quad[i, j] = quad[j, i] = F_STAR * math.sqrt(quad[i, i] * quad[j, j])
    quad[3, 3], quad[4, 4], quad[5, 5] = 1.0 / S12**2, 1.0 / S23**2, 1.0 / S13**2
    return linear, quad


def closed_form_ratio(material_stress: tuple[float, ...] | NDArray[np.float64]) -> float:
    """Positive root of a R^2 + b R - 1 = 0 for one material-axis stress (MPa)."""
    s = np.asarray(material_stress, dtype=np.float64)
    linear, quad = tsai_wu_coefficients()
    a = float(s @ quad @ s)
    b = float(linear @ s)
    return (-b + math.sqrt(b * b + 4.0 * a)) / (2.0 * a)


def closed_form_mode(material_stress: tuple[float, ...]) -> FailureMode:
    """The term with the largest maximum-stress ratio (ADR-010 D2)."""
    s = material_stress
    normal = (
        (T1, C1, FailureMode.ALONG_RASTER_TENSION, FailureMode.ALONG_RASTER_COMPRESSION),
        (T2, C2, FailureMode.ACROSS_ROADS_TENSION, FailureMode.ACROSS_ROADS_COMPRESSION),
        (T3, C3, FailureMode.ACROSS_LAYERS_TENSION, FailureMode.ACROSS_LAYERS_COMPRESSION),
    )
    terms = [
        (s[k] / t if s[k] >= 0 else -s[k] / c, tension if s[k] >= 0 else compression)
        for k, (t, c, tension, compression) in enumerate(normal)
    ]
    terms += [
        (abs(s[3]) / S12, FailureMode.IN_LAYER_SHEAR),
        (abs(s[4]) / S23, FailureMode.INTERLAYER_SHEAR),
        (abs(s[5]) / S13, FailureMode.INTERLAYER_SHEAR),
    ]
    return max(terms, key=lambda term: term[0])[1]


def exact_ratio(case: Case) -> float:
    return closed_form_ratio(case.material_stress)


def test_closed_form_is_consistent() -> None:
    """The closed form against the ADR table, the special cases and the frame convention."""
    for case in CASES:
        r = exact_ratio(case)
        # The ADR table at its printed precision.
        assert r == pytest.approx(ADR_RATIOS[case.name], abs=5e-7), case.name
        # R is a root of a R^2 + b R = 1 and scales as 1 / load.
        s = np.array(case.material_stress)
        linear, quad = tsai_wu_coefficients()
        assert math.isclose(float(r * r * (s @ quad @ s) + r * (linear @ s)), 1.0, rel_tol=1e-13)
        assert math.isclose(closed_form_ratio(2.0 * s), r / 2.0, rel_tol=1e-13)
        # Even in tau_12: a raster at -theta gives the same R.
        flipped = s.copy()
        flipped[3] = -flipped[3]
        assert math.isclose(closed_form_ratio(flipped), r, rel_tol=1e-13), case.name
        assert closed_form_mode(case.material_stress) is case.mode, case.name

    # Uniaxial along a material axis: R = T / sigma.
    assert math.isclose(exact_ratio(RASTER_CASES[0]), T1 / STRESS, rel_tol=1e-13)
    assert math.isclose(exact_ratio(RASTER_CASES[-1]), T2 / STRESS, rel_tol=1e-13)
    assert math.isclose(exact_ratio(SIDE_CASE), T3 / STRESS, rel_tol=1e-13)
    # With F* = -1/2 the quadratic form is positive semidefinite, not definite: the normal block
    # vanishes on s_i = sqrt(T_i C_i) (equal scaled stresses, the von Mises analogue of
    # hydrostatic). Every case here has a > 0.
    quad = tsai_wu_coefficients()[1]
    assert np.linalg.eigvalsh(quad).min() > -1e-15 * np.abs(quad).max()
    null = np.zeros(6)
    null[:3] = np.sqrt(np.array([T1, T2, T3]) * np.array([C1, C2, C3]))
    assert abs(null @ quad @ null) < 1e-14
    for case in CASES:
        s = np.array(case.material_stress)
        assert s @ quad @ s > 0, case.name

    # The material-axis stresses are those of the print-setup frame: sigma' = R sigma R^T.
    from tessera.printsetup import material_frame

    sigma = np.zeros((3, 3))
    sigma[0, 0] = STRESS
    for case in CASES:
        rot = material_frame(case.build, case.raster_deg).rotation
        assert np.allclose(rot[2], case.build, atol=1e-15)
        local = rot @ sigma @ rot.T
        voigt = [local[0, 0], local[1, 1], local[2, 2], local[0, 1], local[1, 2], local[2, 0]]
        assert np.allclose(voigt, case.material_stress, rtol=0.0, atol=1e-12), case.name


# --- Solver runs ----------------------------------------------------------------------------


class Outcome(NamedTuple):
    lowest: float
    clear_lowest: float | None
    face_lowest: tuple[float, ...]
    element_nodal: NDArray[np.float64]  # every element-nodal R, (M * 10,)
    mode: FailureMode
    n_nodes: int

    def max_error(self, exact: float) -> float:
        """max |R / R_exact - 1| over every element-nodal ratio and every reported minimum."""
        reported = [self.lowest, *self.face_lowest]
        if self.clear_lowest is not None:
            reported.append(self.clear_lowest)
        values = np.concatenate([self.element_nodal, np.array(reported)])
        return float(np.abs(values / exact - 1.0).max())


def _outcome(job: Any, result: Any) -> Outcome:
    from tessera.post import element_nodal_stresses, strength_ratio

    strength = result.strength
    assert strength is not None, "no strength result: the solve did not give point stresses"
    assert result.solution.point_stresses is not None
    (section,) = job.sections
    (material,) = job.materials
    nodal = element_nodal_stresses(job.mesh, result.solution.point_stresses)
    ratios = strength_ratio(nodal.reshape(-1, 6), material, section.orientation)
    return Outcome(
        lowest=strength.lowest.ratio,
        clear_lowest=None if strength.clear_lowest is None else strength.clear_lowest.ratio,
        face_lowest=tuple(s.ratio for s in strength.face_lowest),
        element_nodal=ratios,
        mode=strength.lowest.mode,
        n_nodes=job.mesh.n_nodes,
    )


def _face_at(mesh: Any, axis: int, value: float) -> str:
    """Name of the mesh surface lying on the plane x_axis = value."""
    names = [
        name
        for name, tris in mesh.surfaces.items()
        if np.abs(mesh.nodes[np.unique(tris), axis] - value).max() < 1e-9 * SIDE
    ]
    assert len(names) == 1, f"expected one face on {'xyz'[axis]} = {value}, got {names}"
    return names[0]


def _corner(nodes: NDArray[np.float64], point: Vec3) -> NDArray[np.int64]:
    distance = np.linalg.norm(nodes - np.array(point), axis=1)
    node = int(np.argmin(distance))
    assert distance[node] < 1e-9 * SIDE, f"no mesh node at {point}"
    return np.array([node], dtype=np.int64)


def _raster_job(case: Case, step_path: Path, mesh_size: float) -> Any:
    """SolveJob for an (a) case: equal and opposite x-face forces, 3-2-1 corner supports."""
    from tessera.materials import get_material
    from tessera.mesh import mesh_step
    from tessera.model import FixedSupport, NodalLoads, Section, SolveJob
    from tessera.physics import uniform_face_load
    from tessera.printsetup import material_frame

    mesh = mesh_step(step_path, part_id="cube", size=mesh_size)
    material = get_material(MATERIAL.id, [MATERIAL])
    assert material.strength is not None
    loads = []
    for value, sign in ((SIDE, 1.0), (0.0, -1.0)):
        tris = mesh.surfaces[_face_at(mesh, 0, value)]
        nodes, forces = uniform_face_load(mesh.nodes, tris, np.array([sign * STRESS * AREA, 0, 0]))
        loads.append(NodalLoads(nodes=nodes, forces=forces))
    node_sets = dict(mesh.node_sets)
    node_sets["support:origin"] = _corner(mesh.nodes, (0.0, 0.0, 0.0))
    node_sets["support:x-axis"] = _corner(mesh.nodes, (SIDE, 0.0, 0.0))
    node_sets["support:y-axis"] = _corner(mesh.nodes, (0.0, SIDE, 0.0))
    supports: tuple[tuple[str, tuple[int, ...]], ...] = (
        ("support:origin", (1, 2, 3)),
        ("support:x-axis", (2, 3)),
        ("support:y-axis", (3,)),
    )
    return SolveJob(
        mesh=dataclasses.replace(mesh, node_sets=node_sets),
        materials=(material,),
        sections=(
            Section(
                element_set="part:cube",
                material=material.name,
                orientation=material_frame(case.build, case.raster_deg),
            ),
        ),
        supports=tuple(FixedSupport(name, dofs) for name, dofs in supports),
        loads=tuple(loads),
        point_stresses=True,
    )


def _solve_raster(
    case: Case,
    step_path: Path,
    mesh_size: float,
    workdir: Path,
    assert_unstructured: Callable[[np.ndarray, np.ndarray], None],
) -> Outcome:
    from tessera.post import summarize
    from tessera.solvers import CalculixAdapter

    job = _raster_job(case, step_path, mesh_size)
    assert_unstructured(job.mesh.nodes, job.mesh.elements)
    workdir.mkdir(parents=True, exist_ok=True)
    solution = CalculixAdapter(timeout=600).run(job, workdir)
    result = summarize(job, solution, workdir=workdir)
    # The 3-2-1 support holds exactly the six rigid-body motions.
    u = solution.displacements
    held = [("support:origin", (0, 1, 2)), ("support:x-axis", (1, 2)), ("support:y-axis", (2,))]
    for name, dofs in held:
        assert np.abs(u[job.mesh.node_sets[name][0], list(dofs)]).max() == 0.0, name
    return _outcome(job, result)


def _face_with_normal(faces: list[Any], normal: Vec3) -> int:
    matches = [
        f.index for f in faces if f.normal is not None and np.allclose(f.normal, normal, atol=1e-6)
    ]
    assert len(matches) == 1, f"expected one face with normal {normal}, got {matches}"
    return matches[0]


class _KeepJob:
    """CalculiX, remembering the last job so its mesh and section can be read back."""

    name = "calculix"

    def __init__(self) -> None:
        from tessera.solvers import CalculixAdapter

        self._inner = CalculixAdapter(timeout=600)
        self.job: Any = None

    def capabilities(self) -> frozenset[Any]:
        return self._inner.capabilities()

    def run(self, job: Any, workdir: Path) -> Any:
        self.job = job
        return self._inner.run(job, workdir)


_ROLLERS: tuple[tuple[Vec3, Literal[1, 2, 3]], ...] = (
    ((-1.0, 0.0, 0.0), 1),
    ((0.0, -1.0, 0.0), 2),
    ((0.0, 0.0, -1.0), 3),
)


def _solve_side(
    step_path: Path,
    mesh_size: float,
    workdir: Path,
    assert_unstructured: Callable[[np.ndarray, np.ndarray], None],
) -> Outcome:
    """(b) through ``tessera.core.solve``: rollers on the coordinate planes, 1000 N on x = L."""
    from tessera.core import describe_faces, edit, new_project, solve
    from tessera.model import FailureCriterion

    project = new_project("to-06-side")
    edit.add_material(project, MATERIAL)
    edit.add_part(project, step_path, part_id="cube", material_id=MATERIAL.id)
    edit.set_print_setup(project, "cube", "fdm", SIDE_CASE.build, SIDE_CASE.raster_deg)
    faces = describe_faces(project, "cube")
    for normal, dof in _ROLLERS:
        edit.add_support(project, "cube", [_face_with_normal(faces, normal)], dofs=(dof,))
    loaded = _face_with_normal(faces, (1.0, 0.0, 0.0))
    edit.add_force_load(project, "cube", [loaded], (STRESS * AREA, 0.0, 0.0))
    edit.set_mesh_size(project, mesh_size)

    adapter = _KeepJob()
    workdir.mkdir(parents=True, exist_ok=True)
    result = solve(project, workdir=workdir, adapter=adapter)
    assert adapter.job.point_stresses, "core did not ask for point stresses"
    assert result.strength is not None
    assert result.strength.criteria == {MATERIAL.id: FailureCriterion.TSAI_WU}
    assert result.strength.assumed == {}
    assert_unstructured(result.mesh.nodes, result.mesh.elements)
    return _outcome(adapter.job, result)


@pytest.fixture(scope="module")
def cube_step(box_step: Callable[..., Path]) -> Path:
    return box_step("off_axis_strength_cube", (0.0, 0.0, 0.0), (SIDE, SIDE, SIDE))


@pytest.fixture(scope="module")
def strength_run(
    cube_step: Path,
    tmp_path_factory: pytest.TempPathFactory,
    assert_unstructured: Callable[[np.ndarray, np.ndarray], None],
) -> Callable[[Case, float], Outcome]:
    """Cached ``(case, mesh_size) -> Outcome``."""
    cache: dict[tuple[str, float], Outcome] = {}

    def run(case: Case, mesh_size: float) -> Outcome:
        key = (case.name, mesh_size)
        if key not in cache:
            workdir = tmp_path_factory.mktemp(f"to-06-{case.name}-h{mesh_size:g}")
            if case is SIDE_CASE:
                cache[key] = _solve_side(cube_step, mesh_size, workdir, assert_unstructured)
            else:
                cache[key] = _solve_raster(case, cube_step, mesh_size, workdir, assert_unstructured)
        return cache[key]

    return run


def test_reference_matches_formula(load_reference: Callable[[str], dict[str, Any]]) -> None:
    """The stored reference can never drift from the closed-form values."""
    ref = load_reference(REFERENCE)
    assert math.isclose(ref["stress_mpa"], STRESS, rel_tol=1e-9)
    assert set(ref["cases"]) == set(_IDS)
    for case in CASES:
        entry = ref["cases"][case.name]
        assert math.isclose(entry["strength_ratio"], exact_ratio(case), rel_tol=1e-9), case.name
        assert entry["mode"] == case.mode.value, case.name
    # The field is exact: the tolerance covers output precision only, never discretisation.
    assert 0.0 < ref["rel_tolerance_strength_ratio"] <= 1e-5


def _check_against_reference(
    case: Case, outcome: Outcome, ref: dict[str, Any], mesh_size: float
) -> None:
    where = f"{case.name}, h = {mesh_size:g}"
    entry = ref["cases"][case.name]
    expected = entry["strength_ratio"]
    tol = ref["rel_tolerance_strength_ratio"]
    assert abs(outcome.lowest / expected - 1.0) <= tol, f"{where}: lowest R {outcome.lowest!r}"
    if outcome.clear_lowest is not None:
        assert abs(outcome.clear_lowest / expected - 1.0) <= tol, (
            f"{where}: clear lowest R {outcome.clear_lowest!r}"
        )
    error = outcome.max_error(expected)
    assert error <= tol, f"{where}: max |R / R_exact - 1| = {error:.3e}"
    assert outcome.mode.value == entry["mode"], f"{where}: mode {outcome.mode}"


@pytest.mark.parametrize("case", CASES, ids=_IDS)
def test_off_axis_strength(
    case: Case,
    strength_run: Callable[[Case, float], Outcome],
    load_reference: Callable[[str], dict[str, Any]],
) -> None:
    # Solve first: fail (not skip) if the solve breaks, before the reference check can skip.
    outcome = strength_run(case, GATE_MESH_SIZE)
    assert outcome.mode is case.mode
    # Every face of the cube reports its lowest ratio. The raster cases flag the nodes near
    # their corner supports, so clear_lowest is a separate check there; the side case has only
    # rollers, nothing is flagged and clear_lowest is lowest itself.
    assert len(outcome.face_lowest) == 6
    assert outcome.clear_lowest is not None
    if case is SIDE_CASE:
        assert outcome.clear_lowest == outcome.lowest

    ref = load_reference(REFERENCE)
    _check_against_reference(case, outcome, ref, GATE_MESH_SIZE)


@pytest.mark.parametrize("case", CASES, ids=_IDS)
def test_off_axis_strength_mesh_convergence(
    case: Case,
    strength_run: Callable[[Case, float], Outcome],
    load_reference: Callable[[str], dict[str, Any]],
) -> None:
    # Solve first: fail (not skip) if the solve breaks, before the reference check can skip.
    for size in CONVERGENCE_MESH_SIZES:
        strength_run(case, size)

    ref = load_reference(REFERENCE)
    exact = exact_ratio(case)
    print(f"\nTO-06 {case.name}, sigma = {STRESS:g} MPa, R exact = {exact:.10f}")
    print(f"{'h [mm]':>8} {'nodes':>8} {'lowest R':>14} {'clear R':>14} {'max |R/R0-1|':>13}")
    node_counts = []
    for size in CONVERGENCE_MESH_SIZES:
        o = strength_run(case, size)
        node_counts.append(o.n_nodes)
        clear = float("nan") if o.clear_lowest is None else o.clear_lowest
        print(
            f"{size:>8g} {o.n_nodes:>8d} {o.lowest:>14.10f} {clear:>14.10f} "
            f"{o.max_error(exact):>13.3e}"
        )

    assert node_counts == sorted(node_counts) and len(set(node_counts)) == len(node_counts), (
        "mesh_size is not refining the mesh"
    )
    # An exact field has no discretisation error: every mesh must sit at the precision floor.
    for size in CONVERGENCE_MESH_SIZES:
        _check_against_reference(case, strength_run(case, size), ref, size)
