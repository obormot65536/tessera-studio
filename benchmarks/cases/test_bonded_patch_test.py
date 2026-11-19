"""Bonded interface patch test, non-matching meshes (TC-01, ADR-011; exact up to output precision).

Written test-first: it fails until ADR-011 PR-5 adds ``edit.place_part`` and ``edit.add_bond``
and lets ``tessera.core.solve`` take several parts, then skips until the maintainer approves
``benchmarks/references/bonded_patch_test.toml``.

Problem (mm, N, MPa per ADR-002)
--------------------------------
Two blocks of side L = 10, bonded on their common face z = L, each from its own STEP file and
meshed on its own (ADR-011 D2: parts are never merged, only joined by the bond):

- ``base``: 0 <= x, y, z <= L, steel, at its place (identity placement).
- ``cap``: a cube 0 <= x, y, z <= L in its own axes, placed with ``edit.place_part``: turned
  90 deg about global z, then moved by (L, 0, L), so it fills 0 <= x, y <= L, L <= z <= 2L. The
  turn exercises ADR-011 D1 (nodes and the material frame move with the part) and makes the
  interface meshes differ (see "The mesh").

The bond joins the base's face z = L and the cap's face z = L (``edit.add_bond``). Loads are
uniform face forces of sigma A, sigma = 10 MPa, spread with consistent nodal forces; supports
are rollers (``dofs``), given in global axes (ADR-011 D1). Linear static, C3D10 in CalculiX.

Cases (cap material on a steel base):

    z-steel   z tension, free sides     rollers u_x = 0 on x = 0 and u_y = 0 on y = 0 (both
                                        parts), u_z = 0 on z = 0 (base); +z force on z = 2L
    x-steel   x tension (in-plane)      same rollers; +x force on x = L of both parts
    z-pla     z tension, confined       rollers on x = 0, x = L, y = 0, y = L of both parts,
                                        u_z = 0 on z = 0; +z force on z = 2L
    z-ortho   z tension, confined       as z-pla, cap orthotropic (TO constants), printed on its
                                        side (build direction part x, raster angle 0)

Materials are project definitions (``edit.add_material``), made-up test constants, not
material data a user would trust (rule 6 does not apply to them): steel E = 210 000 MPa,
nu = 0.30 with a made-up 250 MPa yield strength, only so that the solve asks for
integration-point stresses (ADR-011 D2, ADR-010); a PLA-like isotropic E = 3500 MPa,
nu = 0.36; the catalog's TO constants (E1 3000, E2 2400, E3 1800, G12 1000, G13 700, G23 650 MPa,
nu12 0.35, nu13 = nu23 0.30).

Why the PLA and orthotropic cases are laterally confined
--------------------------------------------------------
ADR-011 asks for "z tension, PLA on steel" and "an orthotropic block on steel". With free sides
(as z-steel) these have **no** uniform-stress solution: under sigma_zz alone the lateral strain
is -nu sigma / E, -1.43e-5 for steel and -1.03e-3 for the PLA, 72 times larger. A tied
interface forces equal in-plane strains on both sides (continuity of u along the face), so the
blocks shear each other near the bond and the stress there is far from uniform (an error of
order sigma, not 1e-4). Rollers on all four sides of both blocks make the in-plane strains zero
in both, which is compatible, and the state is again exactly uniform in each block (uniaxial
strain). With the same material on both sides (z-steel, x-steel) the free-sided field is
compatible and is kept, as measured in ADR-011 D3.

Exact solution
--------------
In each block the stress is uniform, so the displacement is linear, u = eps . x (base) and
u = eps_cap . (x - x0) + eps_base . x0 with x0 = (0, 0, L) (cap). The fields are continuous on
z = L because both blocks have the same in-plane strains there (eps_xx, eps_yy, eps_xy, and
eps_xz = eps_yz = 0), and the traction sigma . e_z is the same on both sides.

- Free sides, load along axis a: sigma = sigma e_a e_a in both blocks; eps = S sigma.
  z-steel: u = (-nu x, -nu y, z) sigma / E; x-steel: u = (x, -nu y, -nu z) sigma / E.
- Confined, z load (eps_xx = eps_yy = 0, no shear): sigma_zz = sigma and sigma_xx, sigma_yy from
  [S_xx S_xy; S_xy S_yy] [sigma_xx; sigma_yy] = -[S_xz; S_yz] sigma, eps_zz = S_zx sigma_xx +
  S_zy sigma_yy + S_zz sigma (S the compliance in global axes). Isotropic: sigma_xx = sigma_yy =
  nu / (1 - nu) sigma and eps_zz = sigma / M with the constrained modulus
  M = E (1 - nu) / ((1 + nu) (1 - 2 nu)): steel 4.285714 MPa, eps_zz = 3.537415e-5; PLA-like
  5.625 MPa, eps_zz = 1.7e-3.
- Orthotropic cap: the print setup gives the part-axis frame (ADR-009 D2; build part x: rows
  e1 = (0, 0, -1), e2 = (0, 1, 0), e3 = (1, 0, 0)). The placement turns it (ADR-011 D1: axes_global
  = axes_part @ R^T): e1 = (0, 0, -1), e2 = (-1, 0, 0), e3 = (0, 1, 0). So the load runs along
  material axis 1 (raster), sigma_xx = sigma_2 = 3.581769 MPa, sigma_yy = sigma_3 = 2.605898 MPa
  and eps_zz = 2.654870e-3. Were the frame not turned with the part, sigma_xx and sigma_yy would
  swap: the case checks D1's frame rotation.

``test_closed_form_is_consistent`` recomputes all of this two ways (global compliance; and
material axes / the constrained-modulus formula) and checks compatibility, traction continuity,
the supports and the face loads.

Sources: patch test of Irons; R. H. MacNeal & R. L. Harder, "A proposed standard set of problems
to test finite element accuracy", Finite Elements in Analysis and Design 1 (1985) 3-20
(**unverified**). Hooke's law and uniaxial strain: S. P. Timoshenko & J. N. Goodier, *Theory of
Elasticity*, 3rd ed., McGraw-Hill 1970, ch. 1 (**unverified** location). Orthotropic
compliance: S. G. Lekhnitskii, *Theory of Elasticity of an Anisotropic Body* (Mir 1981), ch. 1;
R. M. Jones, *Mechanics of Composite Materials*, 2nd ed. (1999), ch. 2 (**unverified**, as
TO-01). Interface conditions of a perfectly bonded bimaterial (continuity of displacement and
traction): any elasticity text, e.g. Timoshenko & Goodier (**unverified** location).

The mesh
--------
The public API has one mesh size for the study (``Study.mesh_size``), so the two blocks cannot
be given different element sizes as in ADR-011 D3's hand-built measurement (3.0 and 2.1 mm).
Measured with ``tessera.mesh.mesh_step``: two 10 mm cubes meshed separately at the same size
give **identical** interface nodes when the cap is only moved (100 % coincident at 3.0, 2.1,
1.3; 85 % at 0.9). Turned 90 deg about z, the cap's interface triangulation differs: 31 / 23 /
16 % of its interface nodes coincide with the base's at 3.0 / 1.3 / 0.9, all on the face's
perimeter, whose 1D division is the same. The case asserts at most 50 % coincide on every mesh,
so the interface stays non-matching. Sizes as TS-03 (gate 3.0, sweep 3.0 / 1.3 / 0.9, none
divides L), each block checked with ``assert_unstructured``.

Measured quantities
-------------------
On every mesh: the six stress components at every integration point of each block (the
solver's own element stresses, ADR-008, so no extrapolation bias) and at every node, against
the block's exact stress, relative to sigma; the three displacement components at every node
against the exact field, relative to the largest exact displacement over the two blocks. Roller
faces are checked to hold, the placement to have moved the cap, and the parts' node sets to be
disjoint. A jump in u across the bond shows as a displacement error on one side, because both
sides are checked against the same continuous field.

Proposed tolerances (ADR-011 test plan; maintainer-owned reference file)
------------------------------------------------------------------------
z-steel, z-pla, z-ortho: 1e-4 for stress (points and nodes) and displacement. x-steel: 1e-5 for
integration-point stress and displacement, 2e-5 for nodal stress, because CalculiX's nodal
stresses carry TS-03's measured -1e-5 extrapolation bias and a 1e-5 nodal gate cannot pass.
ADR-011 D3 measured, at the proposed penalty stiffness on two steel cubes: 1.0e-7 stress and
1.6e-6 displacement (z tension), about 1e-6 (x tension). The PLA-like and orthotropic cases are
not measured yet (ADR-011 PR-3).

Measured through Tessera's public API: not yet; the case fails at ``edit.place_part`` until
ADR-011 PR-5.
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

from tessera.model import (
    IsotropicDefinition,
    IsotropicStrengthDefinition,
    OrthotropicDefinition,
    Project,
)

pytestmark = [
    pytest.mark.benchmark,
    pytest.mark.solver,
    pytest.mark.skipif(shutil.which("ccx") is None, reason="CalculiX (ccx) not on PATH"),
]

REFERENCE = "bonded_patch_test"

type Vec3 = tuple[float, float, float]
type Frame = tuple[Vec3, Vec3, Vec3]
type Turn = tuple[Literal["x", "y", "z"], float]

# Geometry (mm), load (MPa)
SIDE = 10.0  # L: both blocks are L-cubes
AREA = SIDE * SIDE
STRESS = 10.0  # sigma

# Placement of the cap (ADR-011 D6): turned about global z, then moved.
CAP_TURNS: tuple[Turn, ...] = (("z", 90.0),)
CAP_ROTATION = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])  # R for z 90 deg
CAP_TRANSLATION: Vec3 = (SIDE, 0.0, SIDE)
CAP_ORIGIN = np.array([0.0, 0.0, SIDE])  # x0: the base point of the cap's field
BOND_ID = "glue"

GATE_MESH_SIZE = 3.0  # mm: coarsest distinct mesh (does not divide L)
CONVERGENCE_MESH_SIZES = (3.0, 1.3, 0.9)  # mm, coarse to fine, none divides L
MAX_COINCIDENT_FRACTION = 0.5  # of the cap's interface nodes on a base interface node

# Made-up test materials (not user-trusted data).
STEEL = IsotropicDefinition(
    id="tc-steel",
    youngs_modulus=210_000.0,
    poisson_ratio=0.30,
    strength=IsotropicStrengthDefinition(tensile=250.0),  # only to get point stresses
)
PLA = IsotropicDefinition(id="tc-pla", youngs_modulus=3500.0, poisson_ratio=0.36)
E1, E2, E3 = 3000.0, 2400.0, 1800.0
G12, G13, G23 = 1000.0, 700.0, 650.0
NU12, NU13, NU23 = 0.35, 0.30, 0.30
ORTHO = OrthotropicDefinition(
    id="tc-ortho", e1=E1, e2=E2, e3=E3, nu12=NU12, nu13=NU13, nu23=NU23,
    g12=G12, g13=G13, g23=G23,
)  # fmt: skip

IDENTITY: Frame = ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))
# Rows: material axes in part coordinates for build direction part x, raster 0 (ADR-009 D2,
# written out as in TO-01 rather than taken from tessera.printsetup).
ON_SIDE: Frame = ((0.0, 0.0, -1.0), (0.0, 1.0, 0.0), (1.0, 0.0, 0.0))
BUILD_X: Vec3 = (1.0, 0.0, 0.0)

_VOIGT = ((0, 0), (1, 1), (2, 2), (0, 1), (1, 2), (2, 0))  # xx yy zz xy yz zx


@dataclass(frozen=True)
class Block:
    """A block's material and, for an orthotropic one, its print frame in part axes."""

    material: IsotropicDefinition | OrthotropicDefinition
    build: Vec3 | None = None
    part_frame: Frame = IDENTITY


@dataclass(frozen=True)
class Case:
    name: str
    load_axis: int  # 0 = x, 2 = z
    confined: bool  # rollers on all four sides of both blocks
    cap: Block


BASE = Block(STEEL)
CASES = (
    Case("z-steel", 2, False, Block(STEEL)),
    Case("x-steel", 0, False, Block(STEEL)),
    Case("z-pla", 2, True, Block(PLA)),
    Case("z-ortho", 2, True, Block(ORTHO, BUILD_X, ON_SIDE)),
)
_IDS = [c.name for c in CASES]


# --- closed form -----------------------------------------------------------------------------


def material_compliance(material: IsotropicDefinition | OrthotropicDefinition) -> np.ndarray:
    """6 x 6 compliance (1/MPa), material axes, Voigt xx yy zz xy yz zx, engineering shear."""
    s = np.zeros((6, 6))
    if isinstance(material, IsotropicDefinition):
        e, nu = material.youngs_modulus, material.poisson_ratio
        s[:3, :3] = -nu / e
        np.fill_diagonal(s[:3, :3], 1.0 / e)
        s[3, 3] = s[4, 4] = s[5, 5] = 2.0 * (1.0 + nu) / e
        return s
    m = material
    s[0, 0], s[1, 1], s[2, 2] = 1.0 / m.e1, 1.0 / m.e2, 1.0 / m.e3
    s[0, 1] = s[1, 0] = -m.nu12 / m.e1
    s[0, 2] = s[2, 0] = -m.nu13 / m.e1
    s[1, 2] = s[2, 1] = -m.nu23 / m.e2
    s[3, 3], s[4, 4], s[5, 5] = 1.0 / m.g12, 1.0 / m.g23, 1.0 / m.g13
    return s


def _to_voigt(tensor: np.ndarray, shear_factor: float) -> np.ndarray:
    return np.array([tensor[a, b] * (1.0 if a == b else shear_factor) for a, b in _VOIGT])


def _from_voigt(voigt: np.ndarray, shear_factor: float) -> np.ndarray:
    t = np.zeros((3, 3))
    for k, (a, b) in enumerate(_VOIGT):
        t[a, b] = t[b, a] = voigt[k] * (1.0 if a == b else 1.0 / shear_factor)
    return t


def global_frame(block: Block, rotation: np.ndarray) -> np.ndarray:
    """Rows: material axes in global coordinates, axes_part @ R^T (ADR-011 D1)."""
    return np.array(block.part_frame) @ rotation.T


def global_compliance(block: Block, rotation: np.ndarray) -> np.ndarray:
    """6 x 6 compliance in global axes (Voigt, engineering shear), column by column."""
    q = global_frame(block, rotation)
    s_mat = material_compliance(block.material)
    columns = []
    for k in range(6):
        unit = np.zeros(6)
        unit[k] = 1.0
        stress_mat = _to_voigt(q @ _from_voigt(unit, 1.0) @ q.T, 1.0)
        strain_mat = _from_voigt(s_mat @ stress_mat, 2.0)
        columns.append(_to_voigt(q.T @ strain_mat @ q, 2.0))
    return np.array(columns).T


def _rotation(which: str) -> np.ndarray:
    return CAP_ROTATION if which == "cap" else np.eye(3)


def _block(case: Case, which: str) -> Block:
    return case.cap if which == "cap" else BASE


def exact_stress(case: Case, which: str) -> np.ndarray:
    """(6,) exact Voigt stress of the block (``"base"`` or ``"cap"``), MPa."""
    s = global_compliance(_block(case, which), _rotation(which))
    sigma = np.zeros(6)
    sigma[case.load_axis] = STRESS
    if case.confined:
        # eps_xx = eps_yy = 0, shears 0: solve for sigma_xx, sigma_yy.
        sigma[:2] = np.linalg.solve(s[:2, :2], -s[:2, 2] * STRESS)
    return sigma


def exact_strain(case: Case, which: str) -> np.ndarray:
    """(3, 3) exact tensor strain of the block."""
    s = global_compliance(_block(case, which), _rotation(which))
    return _from_voigt(s @ exact_stress(case, which), 2.0)


def exact_displacement(case: Case, which: str, points: np.ndarray) -> np.ndarray:
    """(N, 3) exact displacements (mm) of points of the block, in global axes."""
    if which == "base":
        return points @ exact_strain(case, "base").T
    base_at_origin = exact_strain(case, "base") @ CAP_ORIGIN
    return (points - CAP_ORIGIN) @ exact_strain(case, "cap").T + base_at_origin


def _box_corners(low: Vec3, high: Vec3) -> np.ndarray:
    return np.array([[x, y, z] for x in (low[0], high[0]) for y in (low[1], high[1])
                     for z in (low[2], high[2])])  # fmt: skip


BASE_BOX: tuple[Vec3, Vec3] = ((0.0, 0.0, 0.0), (SIDE, SIDE, SIDE))
CAP_BOX: tuple[Vec3, Vec3] = ((0.0, 0.0, SIDE), (SIDE, SIDE, 2.0 * SIDE))


def max_displacement(case: Case) -> float:
    """Largest |u_exact| over both blocks (mm); fields are linear, so it is at a corner."""
    base = exact_displacement(case, "base", _box_corners(*BASE_BOX))
    cap = exact_displacement(case, "cap", _box_corners(*CAP_BOX))
    return float(np.linalg.norm(np.vstack([base, cap]), axis=1).max())


def _supports(case: Case) -> list[tuple[str, Vec3, Literal[1, 2, 3]]]:
    """(part, outward global normal, dof): rollers, in global axes."""
    rollers: list[tuple[str, Vec3, Literal[1, 2, 3]]] = []
    for part in ("base", "cap"):
        rollers.append((part, (-1.0, 0.0, 0.0), 1))
        rollers.append((part, (0.0, -1.0, 0.0), 2))
        if case.confined:
            rollers.append((part, (1.0, 0.0, 0.0), 1))
            rollers.append((part, (0.0, 1.0, 0.0), 2))
    rollers.append(("base", (0.0, 0.0, -1.0), 3))
    return rollers


def _loads(case: Case) -> list[tuple[str, Vec3, Vec3]]:
    """(part, outward global normal, total force in N)."""
    if case.load_axis == 2:
        return [("cap", (0.0, 0.0, 1.0), (0.0, 0.0, STRESS * AREA))]
    return [(part, (1.0, 0.0, 0.0), (STRESS * AREA, 0.0, 0.0)) for part in ("base", "cap")]


# --- solve through the public API ------------------------------------------------------------


def _place_part(project: Project, part_id: str, translation: Vec3, turns: tuple[Turn, ...]) -> None:
    """ADR-011 D6 ``edit.place_part``; an AttributeError until ADR-011 PR-5 adds it."""
    from tessera.core import edit

    edit.place_part(project, part_id, translation=translation, rotations=turns)  # pyright: ignore[reportAttributeAccessIssue]


def _add_bond(
    project: Project, part_a: str, faces_a: list[int], part_b: str, faces_b: list[int]
) -> None:
    """ADR-011 D6 ``edit.add_bond``; an AttributeError until ADR-011 PR-5 adds it."""
    from tessera.core import edit

    edit.add_bond(project, part_a, faces_a, part_b, faces_b, connection_id=BOND_ID)  # pyright: ignore[reportAttributeAccessIssue]


def _face_with_normal(faces: list[Any], normal: np.ndarray) -> int:
    matches = [
        f.index for f in faces if f.normal is not None and np.allclose(f.normal, normal, atol=1e-6)
    ]
    assert len(matches) == 1, f"expected one face with normal {normal}, got {matches}"
    return matches[0]


class Errors(NamedTuple):
    point_stress: float  # max |sigma_ip - sigma_exact| / sigma, all components
    nodal_stress: float  # max |sigma_node - sigma_exact| / sigma, all components
    displacement: float  # max |u - u_exact| / max |u_exact|
    n_nodes: int
    coincident: float  # fraction of the cap's interface nodes on a base interface node


def _part_nodes(result: Any, part: str) -> tuple[np.ndarray, np.ndarray]:
    """(element indices, node indices) of a part in the merged mesh (``part:<id>``, D2)."""
    elements = np.asarray(result.mesh.element_sets[f"part:{part}"])
    nodes = np.unique(result.mesh.elements[elements])
    return elements, nodes


def _solve_case(
    case: Case,
    steps: dict[str, Path],
    mesh_size: float,
    workdir: Path,
    assert_unstructured: Callable[[np.ndarray, np.ndarray], None],
) -> Errors:
    from tessera.core import describe_faces, edit, new_project, solve

    project = new_project(f"tc-01-{case.name}")
    for material in {
        BASE.material.id: BASE.material,
        case.cap.material.id: case.cap.material,
    }.values():
        edit.add_material(project, material)
    edit.add_part(project, steps["base"], part_id="base", material_id=BASE.material.id)
    edit.add_part(project, steps["cap"], part_id="cap", material_id=case.cap.material.id)
    if case.cap.build is not None:
        edit.set_print_setup(project, "cap", "fdm", case.cap.build, 0.0)
    # Faces in part axes, read before placing; a global normal n is R^T n in part axes.
    faces = {part: describe_faces(project, part) for part in ("base", "cap")}

    def face(part: str, global_normal: Vec3) -> int:
        return _face_with_normal(faces[part], _rotation(part).T @ np.array(global_normal))

    _place_part(project, "cap", CAP_TRANSLATION, CAP_TURNS)
    _add_bond(project, "base", [face("base", (0.0, 0.0, 1.0))], "cap", [face("cap", (0, 0, -1.0))])
    for part, normal, dof in _supports(case):
        edit.add_support(project, part, [face(part, normal)], dofs=(dof,))
    for part, normal, force in _loads(case):
        edit.add_force_load(project, part, [face(part, normal)], force)
    edit.set_mesh_size(project, mesh_size)

    workdir.mkdir(parents=True, exist_ok=True)
    return _measure(case, solve(project, workdir=workdir), mesh_size, assert_unstructured)


def _measure(
    case: Case,
    result: Any,
    mesh_size: float,
    assert_unstructured: Callable[[np.ndarray, np.ndarray], None],
) -> Errors:
    """Errors of a solved case against the exact fields, after checking the setup held."""
    nodes = result.mesh.nodes
    u = result.solution.displacements
    points = result.solution.point_stresses
    assert points is not None, "a part with a strength must give integration-point stresses"
    tol = 1e-6 * SIDE
    point_error = nodal_error = displacement_error = 0.0
    interface: dict[str, np.ndarray] = {}
    part_nodes: dict[str, np.ndarray] = {}
    for part, (low, high) in (("base", BASE_BOX), ("cap", CAP_BOX)):
        elements, ids = _part_nodes(result, part)
        part_nodes[part] = ids
        assert_unstructured(nodes, result.mesh.elements[elements])
        x = nodes[ids]
        assert np.all(x >= np.array(low) - tol) and np.all(x <= np.array(high) + tol), (
            f"{part}: nodes outside its placed box"
        )
        interface[part] = x[np.abs(x[:, 2] - SIDE) < tol]
        assert len(interface[part]) >= 6, f"{part}: no nodes on the bonded face"

        sigma_ip = points.stresses[elements]
        assert not np.isnan(sigma_ip).any(), f"{part}: missing integration-point stresses"
        sigma_ex = exact_stress(case, part)
        point_error = max(point_error, float(np.abs(sigma_ip - sigma_ex).max() / STRESS))
        sigma_node = result.solution.stresses[ids]
        nodal_error = max(nodal_error, float(np.abs(sigma_node - sigma_ex).max() / STRESS))
        du = u[ids] - exact_displacement(case, part, x)
        displacement_error = max(displacement_error, float(np.abs(du).max()))

        for roller_part, normal, dof in _supports(case):
            if roller_part != part:
                continue
            axis = dof - 1
            plane = low[axis] if normal[axis] < 0 else high[axis]
            held = ids[np.abs(x[:, axis] - plane) < tol]
            assert len(held) >= 6, f"{part}: no nodes on the roller face {normal}"
            assert np.abs(u[held, axis]).max() < 1e-12, f"{part}: roller face {normal} moves"

    assert not np.intersect1d(part_nodes["base"], part_nodes["cap"]).size, "parts share nodes"
    gaps = np.linalg.norm(interface["cap"][:, None, :] - interface["base"][None, :, :], axis=2)
    coincident = float((gaps.min(axis=1) < 1e-6).mean())
    assert coincident <= MAX_COINCIDENT_FRACTION, (
        f"h = {mesh_size:g}: {coincident:.0%} of the interface nodes match: not a non-matching test"
    )
    return Errors(
        point_stress=point_error,
        nodal_stress=nodal_error,
        displacement=displacement_error / max_displacement(case),
        n_nodes=result.mesh.n_nodes,
        coincident=coincident,
    )


@pytest.fixture(scope="module")
def block_steps(box_step: Callable[..., Path]) -> dict[str, Path]:
    """Both blocks as L-cubes at the origin; the cap is placed by ``place_part``."""
    return {
        "base": box_step("bonded_patch_base", (0.0, 0.0, 0.0), (SIDE, SIDE, SIDE)),
        "cap": box_step("bonded_patch_cap", (0.0, 0.0, 0.0), (SIDE, SIDE, SIDE)),
    }


@pytest.fixture(scope="module")
def patch(
    block_steps: dict[str, Path],
    tmp_path_factory: pytest.TempPathFactory,
    assert_unstructured: Callable[[np.ndarray, np.ndarray], None],
) -> Callable[[Case, float], Errors]:
    """Cached ``(case, mesh_size) -> Errors``."""
    cache: dict[tuple[str, float], Errors] = {}

    def run(case: Case, mesh_size: float) -> Errors:
        key = (case.name, mesh_size)
        if key not in cache:
            workdir = tmp_path_factory.mktemp(f"tc-01-{case.name}-h{mesh_size:g}")
            cache[key] = _solve_case(case, block_steps, mesh_size, workdir, assert_unstructured)
        return cache[key]

    return run


# --- tests -----------------------------------------------------------------------------------


def _isotropic(block: Block) -> tuple[float, float]:
    assert isinstance(block.material, IsotropicDefinition)
    return block.material.youngs_modulus, block.material.poisson_ratio


def test_closed_form_is_consistent() -> None:
    """The exact fields are compatible across the bond and satisfy every boundary condition."""
    # The placement maps the cap's own cube onto L <= z <= 2L.
    r = CAP_ROTATION
    assert np.allclose(r @ r.T, np.eye(3)) and math.isclose(np.linalg.det(r), 1.0)
    angle = math.radians(CAP_TURNS[0][1])
    assert np.allclose(r[:2, :2], [[math.cos(angle), -math.sin(angle)],
                                   [math.sin(angle), math.cos(angle)]])  # fmt: skip
    placed = _box_corners(*BASE_BOX) @ r.T + np.array(CAP_TRANSLATION)
    assert np.allclose(np.sort(placed, axis=0), np.sort(_box_corners(*CAP_BOX), axis=0))

    for material in (STEEL, PLA, ORTHO):
        s = material_compliance(material)
        assert np.array_equal(s, s.T) and np.linalg.eigvalsh(s).min() > 0.0

    for case in CASES:
        strains = {w: exact_strain(case, w) for w in ("base", "cap")}
        stresses = {w: _from_voigt(exact_stress(case, w), 1.0) for w in ("base", "cap")}
        # Interface z = L: equal in-plane strains (u continuous) and equal traction sigma . e_z.
        assert np.allclose(strains["base"][:, :2], strains["cap"][:, :2], rtol=0, atol=1e-15)
        assert np.allclose(stresses["base"][:, 2], stresses["cap"][:, 2], rtol=0, atol=1e-12)
        ends = np.array([[3.0, 7.0, SIDE], [9.0, 1.0, SIDE]])
        assert np.allclose(exact_displacement(case, "base", ends),
                           exact_displacement(case, "cap", ends), rtol=0, atol=1e-15)  # fmt: skip
        for which, (low, high) in (("base", BASE_BOX), ("cap", CAP_BOX)):
            sigma = stresses[which]
            # Back through the global compliance: the strain gives the stress.
            back = np.linalg.solve(global_compliance(_block(case, which), _rotation(which)),
                                   _to_voigt(strains[which], 2.0))  # fmt: skip
            assert np.allclose(back, exact_stress(case, which), rtol=0, atol=1e-9 * STRESS)
            assert np.allclose(sigma[[0, 1, 0], [1, 2, 2]], 0.0, atol=1e-12 * STRESS)
            rollers = {(n, d) for p, n, d in _supports(case) if p == which}
            loads = {n: np.array(f) for p, n, f in _loads(case) if p == which}
            for axis in range(3):
                for sign, plane in ((-1.0, low[axis]), (1.0, high[axis])):
                    normal_t = tuple(sign * float(axis == k) for k in range(3))
                    normal: Vec3 = (normal_t[0], normal_t[1], normal_t[2])
                    traction = sigma @ np.array(normal) * AREA
                    if axis == 2 and plane == SIDE:
                        continue  # the bonded face
                    if normal in loads:
                        assert np.allclose(traction, loads[normal], rtol=0, atol=1e-9)
                    elif (normal, axis + 1) in rollers:
                        assert np.allclose(np.delete(traction, axis), 0.0, atol=1e-9)
                        point = np.array([[1.0, 2.0, 3.0]]) + np.array(low)
                        point[0, axis] = plane
                        assert abs(exact_displacement(case, which, point)[0, axis]) < 1e-15
                    else:
                        assert np.allclose(traction, 0.0, atol=1e-9), (case.name, which, normal)

    by_name = {c.name: c for c in CASES}
    # Free sides: the TS-03 fields.
    e, nu = _isotropic(BASE)
    z_steel = exact_displacement(by_name["z-steel"], "cap", np.array([[SIDE, SIDE, 2 * SIDE]]))
    assert np.allclose(z_steel, np.array([[-nu * SIDE, -nu * SIDE, 2 * SIDE]]) * STRESS / e)
    x_steel = exact_displacement(by_name["x-steel"], "cap", np.array([[SIDE, SIDE, 2 * SIDE]]))
    assert np.allclose(x_steel, np.array([[SIDE, -nu * SIDE, -nu * 2 * SIDE]]) * STRESS / e)
    # Confined isotropic: constrained modulus and nu / (1 - nu).
    for case, which in ((by_name["z-pla"], "base"), (by_name["z-pla"], "cap")):
        e, nu = _isotropic(_block(case, which))
        modulus = e * (1 - nu) / ((1 + nu) * (1 - 2 * nu))
        assert math.isclose(exact_strain(case, which)[2, 2], STRESS / modulus, rel_tol=1e-12)
        lateral = nu / (1 - nu) * STRESS
        assert np.allclose(exact_stress(case, which), [lateral, lateral, STRESS, 0, 0, 0])
    # Without confinement, PLA on steel would not be compatible: why z-pla is confined.
    free_pla = Case("free-pla", 2, False, Block(PLA))
    assert not np.allclose(
        exact_strain(free_pla, "base")[0, 0], exact_strain(free_pla, "cap")[0, 0]
    )

    # Confined orthotropic, worked in material axes. The frame turns with the part.
    ortho = by_name["z-ortho"]
    q = global_frame(ortho.cap, CAP_ROTATION)
    assert np.array_equal(q, [[0.0, 0.0, -1.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    s = material_compliance(ORTHO)
    # Axis 1 is along -z: sigma_1 = sigma, eps_2 = eps_3 = 0 (axes 2, 3 along -x, +y).
    s23 = np.linalg.solve(s[1:3, 1:3], -s[1:3, 0] * STRESS)
    eps_1 = s[0, 0] * STRESS + s[0, 1:3] @ s23
    assert math.isclose(exact_strain(ortho, "cap")[2, 2], eps_1, rel_tol=1e-12)
    assert np.allclose(exact_stress(ortho, "cap"), [s23[0], s23[1], STRESS, 0, 0, 0], atol=1e-12)
    assert exact_stress(ortho, "cap")[:2] == pytest.approx((3.581769437, 2.605898123), abs=1e-9)
    assert eps_1 == pytest.approx(2.654870420e-3, rel=1e-9)
    # A frame that did not turn with the part would swap sigma_xx and sigma_yy, visibly.
    unturned = global_compliance(ortho.cap, np.eye(3))
    swapped = np.linalg.solve(unturned[:2, :2], -unturned[:2, 2] * STRESS)
    assert np.allclose(swapped, s23[::-1]) and abs(s23[0] - s23[1]) > 0.05 * STRESS


def test_reference_matches_formula(load_reference: Callable[[str], dict[str, Any]]) -> None:
    """The stored reference can never drift from the closed-form values."""
    ref = load_reference(REFERENCE)
    assert math.isclose(ref["stress_mpa"], STRESS, rel_tol=1e-9)
    assert set(ref["cases"]) == set(_IDS)
    for case in CASES:
        entry = ref["cases"][case.name]
        for which in ("base", "cap"):
            stored = entry[f"stress_{which}_mpa"]
            assert len(stored) == 6
            for got, exact in zip(stored, exact_stress(case, which), strict=True):
                assert math.isclose(got, exact, rel_tol=1e-9, abs_tol=1e-9 * STRESS), case.name
        assert math.isclose(entry["max_displacement_mm"], max_displacement(case), rel_tol=1e-9)
        # A patch test is exact: tolerances only cover output precision and the bond's
        # penalty (sigma / K), never discretisation.
        for key in ("rel_tolerance_stress", "rel_tolerance_nodal_stress"):
            assert 0.0 < entry[key] <= 1e-4, (case.name, key)
        assert 0.0 < entry["rel_tolerance_displacement"] <= 1e-4, case.name


def _check_against_reference(case: Case, errors: Errors, ref: dict[str, Any], h: float) -> None:
    entry = ref["cases"][case.name]
    where = f"{case.name}, h = {h:g}"
    assert errors.point_stress < entry["rel_tolerance_stress"], (
        f"{where}: integration points max |sigma - exact| / sigma = {errors.point_stress:.3e}"
    )
    assert errors.nodal_stress < entry["rel_tolerance_nodal_stress"], (
        f"{where}: nodes max |sigma - exact| / sigma = {errors.nodal_stress:.3e}"
    )
    assert errors.displacement < entry["rel_tolerance_displacement"], (
        f"{where}: max |u - u_exact| / max |u_exact| = {errors.displacement:.3e}"
    )


@pytest.mark.parametrize("case", CASES, ids=_IDS)
def test_bonded_patch_test(
    case: Case,
    patch: Callable[[Case, float], Errors],
    load_reference: Callable[[str], dict[str, Any]],
) -> None:
    # Solve first: fail (not skip) while the assembly API is missing, before the reference
    # check can skip.
    errors = patch(case, GATE_MESH_SIZE)

    ref = load_reference(REFERENCE)
    _check_against_reference(case, errors, ref, GATE_MESH_SIZE)


@pytest.mark.parametrize("case", CASES, ids=_IDS)
def test_bonded_patch_test_mesh_convergence(
    case: Case,
    patch: Callable[[Case, float], Errors],
    load_reference: Callable[[str], dict[str, Any]],
) -> None:
    # Solve first: fail (not skip) while the assembly API is missing.
    for size in CONVERGENCE_MESH_SIZES:
        patch(case, size)

    ref = load_reference(REFERENCE)
    print(f"\nTC-01 {case.name}, sigma = {STRESS:g} MPa, all points and nodes of both blocks")
    print(
        f"{'h [mm]':>8} {'nodes':>8} {'matching':>9} {'|ds| ip/s':>11} {'|ds| node/s':>12} "
        f"{'|du|/u_max':>11}"
    )
    node_counts = []
    for size in CONVERGENCE_MESH_SIZES:
        e = patch(case, size)
        node_counts.append(e.n_nodes)
        print(
            f"{size:>8g} {e.n_nodes:>8d} {e.coincident:>9.0%} {e.point_stress:>11.3e} "
            f"{e.nodal_stress:>12.3e} {e.displacement:>11.3e}"
        )

    assert node_counts == sorted(node_counts) and len(set(node_counts)) == len(node_counts), (
        "mesh_size is not refining the mesh"
    )
    # An exact patch has no discretisation error to converge: every mesh must pass the gate.
    for size in CONVERGENCE_MESH_SIZES:
        _check_against_reference(case, patch(case, size), ref, size)
