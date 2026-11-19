"""Bonded block on part of a plate face (TC-01b, ADR-011; equilibrium and flags, no stress value).

Written test-first: it fails until ADR-011 PR-5 adds ``edit.place_part`` and ``edit.add_bond``
and lets ``tessera.core.solve`` take several parts. The flag test needs no reference; the
force-balance tests then skip until the maintainer approves
``benchmarks/references/bonded_partial_cover.toml``.

Problem (mm, N, MPa per ADR-002)
--------------------------------
A steel plate 0 <= x, y <= W = 30, 0 <= z <= T = 5, and a steel block of side B = 10, each from
its own STEP file at the origin. The block is moved with ``edit.place_part`` by (10, 10, T), so
its bottom face covers the middle of the plate's top face: 10 <= x, y <= 20, z = T. The bond
joins the plate's top face and the block's bottom face (``edit.add_bond``, the plate named
first, so the slave choice cannot come from the argument order). The plate's bottom face is held
in all directions; the block's top face carries a uniform force F = 1000 N in +z (10 MPa).
Steel is a project material, E = 210 000 MPa, nu = 0.30, with a made-up 250 MPa yield strength
only so that the solve asks for integration-point stresses (ADR-011 D2, ADR-010).

What ADR-011 expects (D3, D5)
-----------------------------
- Only the block's face is covered by the other face, so it is the slave, whatever the mesh
  sizes; the plate's nodes outside the footprint stay unbonded. The public result does not name
  the slave (``StudyResult`` carries no job), so the case checks what that choice implies: the
  solve is accepted (were the plate's whole face the bonded side, the gap check would refuse it:
  its nodes outside the footprint are up to 10 mm from the block), and every node of the block's
  face moves with the plate face under it (no jump).
- The bond's boundary is the footprint's perimeter. ``FlagKind.BOND_EDGE`` (bond ``glue``) flags
  nodes within the flag distance d of it, on both parts (the inside corner of the block on the
  plate, which neither part shows on its own, lies there), and nothing in the middle of the bond.

Equilibrium check (mean-stress theorem)
---------------------------------------
For a body without body forces, with sigma_kj,k = 0 and the divergence theorem,

    int_V sigma_ij dV = int_V sigma_kj d(x_i)/dx_k dV = oint_S x_i t_j dS.

Block: t_z = F / B^2 on its top z = T + B and the bond's traction on its bottom z = T, which must
total -F for the block to be in equilibrium, so int sigma_zz dV = F (T + B) - F T = F B
(10 000 N mm), whatever the bond's traction distribution. Plate: +F through the bond at z = T,
the reactions at z = 0 add nothing: int sigma_zz dV = F T (5 000 N mm). Together they say that the
whole load passes through the bond into the plate and the support. The finite element solution
satisfies the same identities exactly (take the virtual displacement z e_z, which the C3D10
shape functions contain): the integral is the four-point sum of the element's own
integration-point stresses (``SolveResult.point_stresses``, ADR-008), exact because the stress
of a straight-edged C3D10 is linear.

Source: the mean-stress theorem, e.g. M. E. Gurtin, "The linear theory of elasticity",
*Handbuch der Physik* VIa/2, Springer 1972 (**unverified** location); derived above in one line.

Displacement jump: the displacement of each node of the block's face against the plate's
quadratic (TRI6) displacement at the same point, interpolated in the plate's top-face triangles,
relative to the largest nodal displacement. The penalty bond lets the faces separate by
sigma / K (ADR-011 D3: about 1e-5 of the part's displacement or less).

The mesh
--------
One study mesh size (the public API has no per-part size). The interface does not match: the
plate's face is meshed with the plate's outline and the block's face on its own. Gate 1.5 mm
(flag distance d = 2h: measured on the two meshes merged by hand, h = 1.82 and d = 3.6 mm, so 15
nodes of the bonded face lie more than 1.1 d from its edge and must not be flagged; at 2.0 mm
there are none), sweep 3.0 / 2.0 / 1.5.

Proposed tolerances (maintainer-owned reference file): force balance 1e-3 relative (output
precision, the contact step's iteration residual, and the bond forces acting at slave nodes that
have moved by u, an error of about u_z / B ~ 5e-5); jump 1e-4 of the largest displacement. Not
measured yet (ADR-011 PR-3 and PR-5).
"""

from __future__ import annotations

import shutil
from collections.abc import Callable
from pathlib import Path
from typing import Any, NamedTuple

import numpy as np
import pytest

from tessera.model import IsotropicDefinition, IsotropicStrengthDefinition, Project

pytestmark = [
    pytest.mark.benchmark,
    pytest.mark.solver,
    pytest.mark.skipif(shutil.which("ccx") is None, reason="CalculiX (ccx) not on PATH"),
]

REFERENCE = "bonded_partial_cover"

type Vec3 = tuple[float, float, float]

# Geometry (mm), load (N)
PLATE_WIDTH = 30.0  # W
PLATE_THICKNESS = 5.0  # T
BLOCK_SIDE = 10.0  # B
BLOCK_TRANSLATION: Vec3 = (10.0, 10.0, PLATE_THICKNESS)
FOOTPRINT_LOW = np.array([10.0, 10.0])
FOOTPRINT_HIGH = np.array([20.0, 20.0])
FORCE = 1000.0  # F, +z on the block's top face
BOND_ID = "glue"

GATE_MESH_SIZE = 1.5  # mm
CONVERGENCE_MESH_SIZES = (3.0, 2.0, 1.5)  # mm, coarse to fine
INSIDE, OUTSIDE = 0.9, 1.1  # margins on the flag distance d for the edge-flag checks

STEEL = IsotropicDefinition(
    id="tc-steel",
    youngs_modulus=210_000.0,
    poisson_ratio=0.30,
    strength=IsotropicStrengthDefinition(tensile=250.0),  # only to get point stresses
)

# Local faces of a C3D10 (corners, then the midside nodes of edges (0,1), (1,2), (2,0) of the
# face triangle), from the node order corners 0-3, edges (0,1) (1,2) (2,0) (0,3) (1,3) (2,3).
_TET10_FACES = (
    (0, 1, 2, 4, 5, 6),
    (0, 1, 3, 4, 8, 7),
    (1, 2, 3, 5, 9, 8),
    (0, 2, 3, 6, 9, 7),
)
_TET10_EDGES = ((0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3))


def block_integral() -> float:
    """int sigma_zz dV over the block, N mm: F B."""
    return FORCE * BLOCK_SIDE


def plate_integral() -> float:
    """int sigma_zz dV over the plate, N mm: F T."""
    return FORCE * PLATE_THICKNESS


def footprint_distance(points: np.ndarray) -> np.ndarray:
    """(N,) distance (mm) from each point to the footprint's perimeter at z = T."""
    corners = [
        (FOOTPRINT_LOW[0], FOOTPRINT_LOW[1]),
        (FOOTPRINT_HIGH[0], FOOTPRINT_LOW[1]),
        (FOOTPRINT_HIGH[0], FOOTPRINT_HIGH[1]),
        (FOOTPRINT_LOW[0], FOOTPRINT_HIGH[1]),
    ]
    ends = [np.array([*c, PLATE_THICKNESS]) for c in corners]
    best = np.full(len(points), np.inf)
    for a, b in zip(ends, ends[1:] + ends[:1], strict=True):
        s = np.clip((points - a) @ (b - a) / ((b - a) @ (b - a)), 0.0, 1.0)
        best = np.minimum(best, np.linalg.norm(points - (a + s[:, None] * (b - a)), axis=1))
    return best


def _place_part(project: Project, part_id: str, translation: Vec3) -> None:
    """ADR-011 D6 ``edit.place_part``; an AttributeError until ADR-011 PR-5 adds it."""
    from tessera.core import edit

    edit.place_part(project, part_id, translation=translation, rotations=())  # pyright: ignore[reportAttributeAccessIssue]


def _add_bond(
    project: Project, part_a: str, faces_a: list[int], part_b: str, faces_b: list[int]
) -> None:
    """ADR-011 D6 ``edit.add_bond``; an AttributeError until ADR-011 PR-5 adds it."""
    from tessera.core import edit

    edit.add_bond(project, part_a, faces_a, part_b, faces_b, connection_id=BOND_ID)  # pyright: ignore[reportAttributeAccessIssue]


def _face_with_normal(faces: list[Any], normal: Vec3) -> int:
    matches = [
        f.index for f in faces if f.normal is not None and np.allclose(f.normal, normal, atol=1e-6)
    ]
    assert len(matches) == 1, f"expected one face with normal {normal}, got {matches}"
    return matches[0]


def _solve(steps: dict[str, Path], mesh_size: float, workdir: Path) -> Any:
    from tessera.core import describe_faces, edit, new_project, solve

    project = new_project("tc-01b")
    edit.add_material(project, STEEL)
    edit.add_part(project, steps["plate"], part_id="plate", material_id=STEEL.id)
    edit.add_part(project, steps["block"], part_id="block", material_id=STEEL.id)
    plate = describe_faces(project, "plate")
    block = describe_faces(project, "block")
    top = _face_with_normal(block, (0.0, 0.0, 1.0))
    assert next(f for f in block if f.index == top).area == pytest.approx(BLOCK_SIDE**2)

    _place_part(project, "block", BLOCK_TRANSLATION)
    _add_bond(
        project,
        "plate",
        [_face_with_normal(plate, (0.0, 0.0, 1.0))],
        "block",
        [_face_with_normal(block, (0.0, 0.0, -1.0))],
    )
    edit.add_support(project, "plate", [_face_with_normal(plate, (0.0, 0.0, -1.0))])
    edit.add_force_load(project, "block", [top], (0.0, 0.0, FORCE))
    edit.set_mesh_size(project, mesh_size)

    workdir.mkdir(parents=True, exist_ok=True)
    return solve(project, workdir=workdir)


def _part(result: Any, part: str) -> tuple[np.ndarray, np.ndarray]:
    """(element indices, node indices) of a part in the merged mesh (``part:<id>``, D2)."""
    elements = np.asarray(result.mesh.element_sets[f"part:{part}"])
    return elements, np.unique(result.mesh.elements[elements])


def _sigma_zz_integral(result: Any, elements: np.ndarray) -> float:
    """int sigma_zz dV (N mm): element volume times the mean of its four point stresses."""
    points = result.solution.point_stresses
    assert points is not None, "a part with a strength must give integration-point stresses"
    assert points.stresses.shape[1] == 4, "expected the four-point C3D10 rule (ADR-008)"
    connectivity = result.mesh.elements[elements]
    x = result.mesh.nodes[connectivity]  # (M, 10, 3)
    # Straight edges (flat faces): the midside nodes are midpoints, so the stress is linear in
    # each element and the four-point rule with equal weights integrates it exactly.
    for k, (i, j) in enumerate(_TET10_EDGES):
        assert np.allclose(x[:, 4 + k], 0.5 * (x[:, i] + x[:, j]), rtol=0.0, atol=1e-9)
    edges = x[:, 1:4] - x[:, :1]
    volumes = np.abs(np.linalg.det(edges)) / 6.0
    sigma_zz = points.stresses[elements][:, :, 2]
    assert not np.isnan(sigma_zz).any(), "missing integration-point stresses"
    return float(volumes @ sigma_zz.mean(axis=1))


def _plate_face_displacement(result: Any, plate_elements: np.ndarray, at: np.ndarray) -> np.ndarray:
    """(P, 3) plate displacement at points ``at`` of its top face z = T (TRI6 interpolation)."""
    nodes, u = result.mesh.nodes, result.solution.displacements
    tol = 1e-6 * PLATE_THICKNESS
    triangles = []
    for element in result.mesh.elements[plate_elements]:
        for face in _TET10_FACES:
            tri = element[list(face)]
            if np.all(np.abs(nodes[tri[:3], 2] - PLATE_THICKNESS) < tol):
                triangles.append(tri)
    tris = np.array(triangles)  # (K, 6)
    a, b, c = (nodes[tris[:, k], :2] for k in range(3))  # (K, 2) each
    det = (b[:, 0] - a[:, 0]) * (c[:, 1] - a[:, 1]) - (c[:, 0] - a[:, 0]) * (b[:, 1] - a[:, 1])
    p = at[:, None, :2]  # (P, 1, 2)
    l1 = (p[..., 0] - a[:, 0]) * (c[:, 1] - a[:, 1]) - (c[:, 0] - a[:, 0]) * (p[..., 1] - a[:, 1])
    l2 = (b[:, 0] - a[:, 0]) * (p[..., 1] - a[:, 1]) - (p[..., 0] - a[:, 0]) * (b[:, 1] - a[:, 1])
    l1, l2 = l1 / det, l2 / det
    l0 = 1.0 - l1 - l2
    inside = np.minimum(np.minimum(l0, l1), l2)  # (P, K)
    k = inside.argmax(axis=1)
    rows = np.arange(len(at))
    assert np.all(inside[rows, k] > -1e-9), "a block face node is not over the plate's face"
    w0, w1, w2 = l0[rows, k], l1[rows, k], l2[rows, k]
    shape = np.stack(
        [w0 * (2 * w0 - 1), w1 * (2 * w1 - 1), w2 * (2 * w2 - 1), 4 * w0 * w1, 4 * w1 * w2,
         4 * w2 * w0],
        axis=1,
    )  # fmt: skip
    return np.einsum("pk,pkj->pj", shape, u[tris[k]])


class Balance(NamedTuple):
    block: float  # int sigma_zz dV / (F B) - 1
    plate: float  # int sigma_zz dV / (F T) - 1
    jump: float  # max |u_block - u_plate| on the block's face / max |u|
    n_nodes: int
    peak: float  # raw von Mises peak, MPa (singular at the bond edge: grows)


def _balance(result: Any) -> Balance:
    plate_elements, plate_nodes = _part(result, "plate")
    block_elements, block_nodes = _part(result, "block")
    assert not np.intersect1d(plate_nodes, block_nodes).size, "parts share nodes"
    x = result.mesh.nodes[block_nodes]
    low = np.array([*FOOTPRINT_LOW, PLATE_THICKNESS])
    high = low + BLOCK_SIDE
    assert np.all(x >= low - 1e-6) and np.all(x <= high + 1e-6), "the block was not placed"

    on_face = block_nodes[np.abs(x[:, 2] - PLATE_THICKNESS) < 1e-6 * PLATE_THICKNESS]
    assert len(on_face) >= 6, "no block nodes on the bonded face"
    u = result.solution.displacements
    plate_u = _plate_face_displacement(result, plate_elements, result.mesh.nodes[on_face])
    jump = float(np.abs(u[on_face] - plate_u).max() / np.linalg.norm(u, axis=1).max())
    return Balance(
        block=_sigma_zz_integral(result, block_elements) / block_integral() - 1.0,
        plate=_sigma_zz_integral(result, plate_elements) / plate_integral() - 1.0,
        jump=jump,
        n_nodes=result.mesh.n_nodes,
        peak=result.peak_von_mises.value,
    )


@pytest.fixture(scope="module")
def part_steps(box_step: Callable[..., Path]) -> dict[str, Path]:
    """Plate and block, each at the origin; the block is placed by ``place_part``."""
    return {
        "plate": box_step(
            "bonded_cover_plate", (0.0, 0.0, 0.0), (PLATE_WIDTH, PLATE_WIDTH, PLATE_THICKNESS)
        ),
        "block": box_step(
            "bonded_cover_block", (0.0, 0.0, 0.0), (BLOCK_SIDE, BLOCK_SIDE, BLOCK_SIDE)
        ),
    }


@pytest.fixture(scope="module")
def cover(
    part_steps: dict[str, Path], tmp_path_factory: pytest.TempPathFactory
) -> Callable[[float], Any]:
    """Cached ``mesh_size -> StudyResult``."""
    cache: dict[float, Any] = {}

    def run(mesh_size: float) -> Any:
        if mesh_size not in cache:
            workdir = tmp_path_factory.mktemp(f"tc-01b-h{mesh_size:g}")
            cache[mesh_size] = _solve(part_steps, mesh_size, workdir)
        return cache[mesh_size]

    return run


def test_closed_form_is_consistent() -> None:
    """The placed block sits inside the plate's top face, and the integrals follow from F."""
    low = np.array(BLOCK_TRANSLATION)
    assert np.allclose(low[:2], FOOTPRINT_LOW) and np.allclose(low[:2] + BLOCK_SIDE, FOOTPRINT_HIGH)
    assert low[2] == PLATE_THICKNESS
    assert np.all(FOOTPRINT_LOW > 0.0) and np.all(FOOTPRINT_HIGH < PLATE_WIDTH)
    assert block_integral() == pytest.approx(10_000.0) and plate_integral() == pytest.approx(
        5_000.0
    )
    # The mean-stress theorem on the block's uniform far field: sigma_zz = F / B^2 over B^3.
    assert block_integral() == pytest.approx(FORCE / BLOCK_SIDE**2 * BLOCK_SIDE**3)
    probe = np.array(
        [[15.0, 15.0, PLATE_THICKNESS], [10.0, 12.0, PLATE_THICKNESS], [12.0, 10.0, 0.0]]
    )
    assert np.allclose(footprint_distance(probe), [5.0, 0.0, PLATE_THICKNESS])


def test_reference_matches_formula(load_reference: Callable[[str], dict[str, Any]]) -> None:
    ref = load_reference(REFERENCE)
    assert ref["force_n"] == pytest.approx(FORCE, rel=1e-9)
    assert ref["block_sigma_zz_integral_n_mm"] == pytest.approx(block_integral(), rel=1e-9)
    assert ref["plate_sigma_zz_integral_n_mm"] == pytest.approx(plate_integral(), rel=1e-9)
    assert 0.0 < ref["rel_tolerance_force_balance"] <= 1e-3
    assert 0.0 < ref["rel_tolerance_jump"] <= 1e-4


def test_partial_cover_bond_is_accepted_and_its_edge_flagged(
    cover: Callable[[float], Any],
) -> None:
    """ADR-011 D3 and D5, no reference value: the solve is accepted, nodes near the footprint's
    edge carry ``BOND_EDGE`` on both parts and the middle of the bond does not, and the raw peak
    (the inside corner of the block on the plate) is flagged."""
    from tessera.model import FlagKind, StressFlag

    result = cover(GATE_MESH_SIZE)
    flags = result.node_flags
    d = flags.distance
    edge = StressFlag(FlagKind.BOND_EDGE, bond=BOND_ID)
    assert edge in flags.masks, f"no bond-edge flag; flags are {list(flags.masks)}"
    mask = flags.masks[edge]
    distance = footprint_distance(result.mesh.nodes)

    near = distance <= INSIDE * d
    for part in ("plate", "block"):
        _, ids = _part(result, part)
        assert near[ids].any(), f"no {part} nodes near the bond's edge"
        assert mask[ids][near[ids]].all(), f"{part} nodes near the bond's edge are not flagged"
    far = distance >= OUTSIDE * d
    assert not mask[far].any(), "nodes away from the bond's edge are flagged as at it"
    on_bond = np.abs(result.mesh.nodes[:, 2] - PLATE_THICKNESS) < 1e-6
    middle = (
        on_bond
        & np.all(result.mesh.nodes[:, :2] > FOOTPRINT_LOW, axis=1)
        & np.all(result.mesh.nodes[:, :2] < FOOTPRINT_HIGH, axis=1)
        & far
    )
    assert middle.any(), "the mesh has no nodes in the middle of the bond to check"

    peak = result.peak_von_mises
    assert peak.part in ("plate", "block")
    assert footprint_distance(np.array([peak.location]))[0] <= d, "raw peak not at the bond edge"
    assert edge in peak.flags
    for spot in result.face_peaks:
        assert spot.part in ("plate", "block")
        spot_distance = footprint_distance(np.array([spot.location]))[0]
        if spot_distance <= INSIDE * d:
            assert edge in spot.flags, f"face {spot.face} of {spot.part!r}: {spot.flags}"
        elif spot_distance >= OUTSIDE * d:
            assert edge not in spot.flags, f"face {spot.face} of {spot.part!r}: {spot.flags}"
    clear = result.clear_peak_von_mises
    assert clear is not None and not clear.flagged


def _check_against_reference(balance: Balance, ref: dict[str, Any], mesh_size: float) -> None:
    tol = ref["rel_tolerance_force_balance"]
    assert abs(balance.block) < tol, (
        f"h = {mesh_size:g}: block int sigma_zz off by {balance.block:.3e}"
    )
    assert abs(balance.plate) < tol, (
        f"h = {mesh_size:g}: plate int sigma_zz off by {balance.plate:.3e}"
    )
    assert balance.jump < ref["rel_tolerance_jump"], (
        f"h = {mesh_size:g}: displacement jump across the bond {balance.jump:.3e} of max |u|"
    )


def test_partial_cover_force_balance(
    cover: Callable[[float], Any], load_reference: Callable[[str], dict[str, Any]]
) -> None:
    # Solve first: fail (not skip) while the assembly API is missing.
    balance = _balance(cover(GATE_MESH_SIZE))

    ref = load_reference(REFERENCE)
    _check_against_reference(balance, ref, GATE_MESH_SIZE)


def test_partial_cover_mesh_convergence(
    cover: Callable[[float], Any], load_reference: Callable[[str], dict[str, Any]]
) -> None:
    # Solve first: fail (not skip) while the assembly API is missing.
    rows = [(size, _balance(cover(size))) for size in CONVERGENCE_MESH_SIZES]

    ref = load_reference(REFERENCE)
    print(f"\nTC-01b block on plate, F = {FORCE:g} N; integrals against F B and F T")
    print(f"{'h [mm]':>8} {'nodes':>8} {'block':>11} {'plate':>11} {'jump':>10} {'peak MPa':>9}")
    for size, b in rows:
        print(
            f"{size:>8g} {b.n_nodes:>8d} {b.block:>11.3e} {b.plate:>11.3e} {b.jump:>10.2e} "
            f"{b.peak:>9.2f}"
        )
    node_counts = [b.n_nodes for _, b in rows]
    assert node_counts == sorted(node_counts) and len(set(node_counts)) == len(node_counts), (
        "mesh_size is not refining the mesh"
    )
    # Equilibrium holds on every mesh; the raw peak is singular at the bond edge and only
    # printed.
    for size, b in rows:
        _check_against_reference(b, ref, size)
