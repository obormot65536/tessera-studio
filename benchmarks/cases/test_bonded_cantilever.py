"""Cantilever cut at mid-span into two bonded parts, tip deflection (TC-02, ADR-011; gate 2 %).

Written test-first: it fails until ADR-011 PR-5 adds ``edit.place_part`` and ``edit.add_bond``
and lets ``tessera.core.solve`` take several parts, then skips until the maintainer approves
``benchmarks/references/bonded_cantilever.toml``.

Problem (mm, N, MPa per ADR-002)
--------------------------------
The TS-01 beam (``test_cantilever.py``): steel (``steel-generic``, E = 210 000 MPa,
nu = 0.30), length L = 200 along +x, section b = h = 10 centred on the x axis, root face x = 0
clamped, a total force P = 100 N in -z spread over the tip face x = L. Here it is two parts,
each from its own STEP file, meshed on its own and joined by a bond on the cut x = L / 2:

- ``root``: 0 <= x <= L / 2, at its place.
- ``tip``: a 100 x 10 x 10 box at the origin in its own axes, placed with ``edit.place_part``:
  turned 90 deg about global x, then moved by (L / 2, b / 2, -h / 2), so it fills
  L / 2 <= x <= L, |y|, |z| <= 5. Turning the square section about the beam axis leaves the
  beam unchanged, but the tip half's nodes, faces and supports now go through ADR-011 D1.

The bond joins the root's face x = L / 2 and the tip's face at its own x = 0 (``edit.add_bond``).

Reference
---------
Timoshenko beam, as TS-01: delta = P L^3 / (3 E I) + P L / (kappa G A), I = b h^3 / 12,
A = b h, G = E / (2 (1 + nu)), Cowper's kappa = 10 (1 + nu) / (12 + 11 nu): 1.526724 mm.
Sources: S. P. Timoshenko and J. M. Gere, *Mechanics of Materials* (end-loaded cantilever);
Timoshenko and Goodier, *Theory of Elasticity*, 3rd ed., sec. 21; G. R. Cowper, *J. Appl.
Mech.* 33 (1966) 335-340. Full discussion in ``test_cantilever.py``: the clamped root face makes
the 3D model about 0.35 % stiffer than the beam.

A perfectly bonded cut changes nothing in the exact solution, so the bonded beam must also agree
with the single-part TS-01 beam solved here at the same mesh size, which takes the clamp and
discretisation effects out of the comparison: what is left is the bond's own error.

The mesh
--------
One study mesh size. The cut faces do not match: measured with ``tessera.mesh.mesh_step``,
38 % / 31 % of the tip half's cut-face nodes coincide with the root half's at 2.5 / 1.8 mm
(43 % at 5.0 mm), the rest are elsewhere; the case asserts at most 50 %. Gate 2.5 mm (as TS-01),
sweep 5.0 / 2.5 / 1.8.

Measured quantity: mean u_z of the nodes on the tip face x = L.

Proposed tolerances (maintainer-owned reference file; unchanged from the catalog and ADR-011):
within 2 % of 1.526724 mm (the v0.1 gate), and within 0.5 % of the single-part result at the
same mesh size, at the gate and finer. ADR-011 D3 measured -0.40 % against the reference and
within 0.03 % of the single part at 2.5 / 2.0 mm with every join tried, because the bond's error
stays local. At 5.0 mm (two elements through the depth) the agreement is printed, not gated.
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

from tessera.model import Project

pytestmark = [
    pytest.mark.benchmark,
    pytest.mark.solver,
    pytest.mark.skipif(shutil.which("ccx") is None, reason="CalculiX (ccx) not on PATH"),
]

REFERENCE = "bonded_cantilever"
SINGLE_PART_REFERENCE = "cantilever_tip_deflection"  # TS-01

type Vec3 = tuple[float, float, float]
type Turn = tuple[Literal["x", "y", "z"], float]

# Geometry (mm), load (N), material (MPa): TS-01
LENGTH = 200.0
WIDTH = 10.0  # b, along y
HEIGHT = 10.0  # h, along z (load direction)
CUT = LENGTH / 2.0
TIP_FORCE = 100.0  # total, applied in -z
YOUNGS_MODULUS = 210_000.0  # steel-generic
POISSON_RATIO = 0.30

TIP_TURNS: tuple[Turn, ...] = (("x", 90.0),)
TIP_ROTATION = np.array([[1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, 1.0, 0.0]])  # R for x 90 deg
TIP_TRANSLATION: Vec3 = (CUT, WIDTH / 2.0, -HEIGHT / 2.0)
BOND_ID = "cut"

GATE_MESH_SIZE = 2.5  # mm: 4 elements through the depth, as TS-01
CONVERGENCE_MESH_SIZES = (5.0, 2.5, 1.8)  # mm, coarse to fine
CONVERGENCE_SLACK = 0.002  # allowed relative-error increase between refinements (noise)
MAX_COINCIDENT_FRACTION = 0.5  # of the tip half's cut-face nodes on a root half's node


def timoshenko_tip_deflection(
    length: float = LENGTH,
    width: float = WIDTH,
    height: float = HEIGHT,
    force: float = TIP_FORCE,
    youngs_modulus: float = YOUNGS_MODULUS,
    poisson_ratio: float = POISSON_RATIO,
) -> float:
    """Tip deflection magnitude (mm) of an end-loaded rectangular Timoshenko cantilever."""
    inertia = width * height**3 / 12.0
    area = width * height
    shear_modulus = youngs_modulus / (2.0 * (1.0 + poisson_ratio))
    kappa = 10.0 * (1.0 + poisson_ratio) / (12.0 + 11.0 * poisson_ratio)  # Cowper (1966)
    bending = force * length**3 / (3.0 * youngs_modulus * inertia)
    shear = force * length / (kappa * shear_modulus * area)
    return bending + shear


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


def _face_with_normal(faces: list[Any], normal: Vec3) -> int:
    matches = [
        f.index for f in faces if f.normal is not None and np.allclose(f.normal, normal, atol=1e-6)
    ]
    assert len(matches) == 1, f"expected one face with normal {normal}, got {matches}"
    return matches[0]


def _tip_deflection(result: Any) -> tuple[float, int]:
    """(mean u_z of the tip-face nodes, number of nodes)."""
    nodes = result.mesh.nodes
    on_tip = np.isclose(nodes[:, 0], LENGTH, rtol=0.0, atol=1e-6 * LENGTH)
    assert on_tip.sum() >= 6, "no mesh nodes found on the tip face"
    return float(result.solution.displacements[on_tip, 2].mean()), result.mesh.n_nodes


def _check_root_clamped(result: Any) -> None:
    nodes = result.mesh.nodes
    on_root = np.isclose(nodes[:, 0], 0.0, rtol=0.0, atol=1e-6 * LENGTH)
    assert on_root.sum() >= 6, "no mesh nodes found on the root face"
    assert np.abs(result.solution.displacements[on_root]).max() < 1e-12, "root is not clamped"


def _solve_single(step: Path, mesh_size: float, workdir: Path) -> Any:
    """The TS-01 beam as one part (works today; the comparison baseline)."""
    from tessera.core import describe_faces, edit, new_project, solve

    project = new_project("tc-02-single")
    edit.add_part(project, step, part_id="beam", material_id="steel-generic")
    faces = describe_faces(project, "beam")
    edit.add_support(project, "beam", [_face_with_normal(faces, (-1.0, 0.0, 0.0))])
    edit.add_force_load(
        project, "beam", [_face_with_normal(faces, (1.0, 0.0, 0.0))], (0.0, 0.0, -TIP_FORCE)
    )
    edit.set_mesh_size(project, mesh_size)
    workdir.mkdir(parents=True, exist_ok=True)
    result = solve(project, workdir=workdir)
    _check_root_clamped(result)
    return result


def _solve_bonded(steps: dict[str, Path], mesh_size: float, workdir: Path) -> Any:
    """The beam as two bonded halves, the tip half placed by ``place_part``."""
    from tessera.core import describe_faces, edit, new_project, solve

    project = new_project("tc-02-bonded")
    edit.add_part(project, steps["root"], part_id="root", material_id="steel-generic")
    edit.add_part(project, steps["tip"], part_id="tip", material_id="steel-generic")
    root = describe_faces(project, "root")
    tip = describe_faces(project, "tip")  # part axes; the turn about x keeps +-x normals
    root_end = _face_with_normal(root, (1.0, 0.0, 0.0))
    tip_start = _face_with_normal(tip, (-1.0, 0.0, 0.0))
    tip_end = _face_with_normal(tip, (1.0, 0.0, 0.0))
    assert next(f for f in root if f.index == root_end).centroid[0] == pytest.approx(CUT)
    assert next(f for f in tip if f.index == tip_end).area == pytest.approx(WIDTH * HEIGHT)

    _place_part(project, "tip", TIP_TRANSLATION, TIP_TURNS)
    _add_bond(project, "root", [root_end], "tip", [tip_start])
    edit.add_support(project, "root", [_face_with_normal(root, (-1.0, 0.0, 0.0))])
    edit.add_force_load(project, "tip", [tip_end], (0.0, 0.0, -TIP_FORCE))
    edit.set_mesh_size(project, mesh_size)
    workdir.mkdir(parents=True, exist_ok=True)
    result = solve(project, workdir=workdir)

    _check_root_clamped(result)
    parts = {}
    for part in ("root", "tip"):
        elements = np.asarray(result.mesh.element_sets[f"part:{part}"])
        parts[part] = np.unique(result.mesh.elements[elements])
    assert not np.intersect1d(parts["root"], parts["tip"]).size, "parts share nodes"
    nodes = result.mesh.nodes
    tip_x = nodes[parts["tip"]]
    low, high = np.array([CUT, -WIDTH / 2, -HEIGHT / 2]), np.array([LENGTH, WIDTH / 2, HEIGHT / 2])
    assert np.all(tip_x >= low - 1e-6) and np.all(tip_x <= high + 1e-6), "tip half not placed"
    faces = {}
    for part, ids in parts.items():
        x = nodes[ids]
        faces[part] = x[np.isclose(x[:, 0], CUT, rtol=0.0, atol=1e-6 * LENGTH)][:, 1:]
        assert len(faces[part]) >= 6, f"no {part} nodes on the cut"
    gaps = np.linalg.norm(faces["tip"][:, None, :] - faces["root"][None, :, :], axis=2)
    coincident = float((gaps.min(axis=1) < 1e-6).mean())
    assert coincident <= MAX_COINCIDENT_FRACTION, (
        f"h = {mesh_size:g}: {coincident:.0%} of the cut's nodes match: not a non-matching test"
    )
    return result


@pytest.fixture(scope="module")
def beam_steps(box_step: Callable[..., Path]) -> dict[str, Path]:
    return {
        "single": box_step("cantilever", (0.0, -WIDTH / 2, -HEIGHT / 2), (LENGTH, WIDTH, HEIGHT)),
        "root": box_step("bonded_cantilever_root", (0.0, -WIDTH / 2, -HEIGHT / 2),
                         (CUT, WIDTH, HEIGHT)),
        "tip": box_step("bonded_cantilever_tip", (0.0, 0.0, 0.0), (LENGTH - CUT, WIDTH, HEIGHT)),
    }  # fmt: skip


@pytest.fixture(scope="module")
def deflection(
    beam_steps: dict[str, Path], tmp_path_factory: pytest.TempPathFactory
) -> Callable[[str, float], tuple[float, int]]:
    """Cached ``(kind, mesh_size) -> (mean tip u_z, n_nodes)``, kind ``single`` or ``bonded``."""
    cache: dict[tuple[str, float], tuple[float, int]] = {}

    def run(kind: str, mesh_size: float) -> tuple[float, int]:
        key = (kind, mesh_size)
        if key not in cache:
            workdir = tmp_path_factory.mktemp(f"tc-02-{kind}-h{mesh_size:g}")
            if kind == "single":
                result = _solve_single(beam_steps["single"], mesh_size, workdir)
            else:
                result = _solve_bonded(beam_steps, mesh_size, workdir)
            cache[key] = _tip_deflection(result)
        return cache[key]

    return run


def test_closed_form_is_consistent() -> None:
    """The placed tip half fills the second half of the TS-01 beam."""
    r = TIP_ROTATION
    assert np.allclose(r @ r.T, np.eye(3)) and math.isclose(np.linalg.det(r), 1.0)
    angle = math.radians(TIP_TURNS[0][1])
    assert np.allclose(r[1:, 1:], [[math.cos(angle), -math.sin(angle)],
                                   [math.sin(angle), math.cos(angle)]])  # fmt: skip
    corners = np.array([[x, y, z] for x in (0.0, LENGTH - CUT) for y in (0.0, WIDTH)
                        for z in (0.0, HEIGHT)])  # fmt: skip
    placed = corners @ r.T + np.array(TIP_TRANSLATION)
    assert np.allclose(placed.min(axis=0), [CUT, -WIDTH / 2, -HEIGHT / 2])
    assert np.allclose(placed.max(axis=0), [LENGTH, WIDTH / 2, HEIGHT / 2])
    assert np.allclose(r @ [1.0, 0.0, 0.0], [1.0, 0.0, 0.0])  # the tip half's +x stays +x
    assert timoshenko_tip_deflection() == pytest.approx(1.526724, abs=5e-7)


def test_reference_matches_formula(load_reference: Callable[[str], dict[str, Any]]) -> None:
    """The stored reference can never drift from the closed form or from TS-01's."""
    ref = load_reference(REFERENCE)
    assert math.isclose(
        ref["reference_tip_deflection_mm"], timoshenko_tip_deflection(), rel_tol=1e-9
    )
    single = load_reference(SINGLE_PART_REFERENCE)
    assert math.isclose(
        ref["reference_tip_deflection_mm"], single["reference_tip_deflection_mm"], rel_tol=1e-9
    )
    assert 0.0 < ref["rel_tolerance"] <= 0.02  # the v0.1 gate
    assert 0.0 < ref["rel_tolerance_to_single_part"] <= 0.005


def test_bonded_cantilever_tip_deflection(
    deflection: Callable[[str, float], tuple[float, int]],
    load_reference: Callable[[str], dict[str, Any]],
) -> None:
    # Solve first: fail (not skip) while the assembly API is missing.
    bonded, _ = deflection("bonded", GATE_MESH_SIZE)
    single, _ = deflection("single", GATE_MESH_SIZE)

    ref = load_reference(REFERENCE)
    assert bonded < 0.0, "tip must deflect in the load direction (-z)"
    assert -bonded == pytest.approx(ref["reference_tip_deflection_mm"], rel=ref["rel_tolerance"])
    assert bonded == pytest.approx(single, rel=ref["rel_tolerance_to_single_part"])


def test_bonded_cantilever_mesh_convergence(
    deflection: Callable[[str, float], tuple[float, int]],
    load_reference: Callable[[str], dict[str, Any]],
) -> None:
    # Solve first: fail (not skip) while the assembly API is missing.
    rows = []
    for size in CONVERGENCE_MESH_SIZES:
        bonded, n_nodes = deflection("bonded", size)
        single, _ = deflection("single", size)
        rows.append((size, n_nodes, -bonded, -single))

    ref = load_reference(REFERENCE)
    expected = ref["reference_tip_deflection_mm"]
    tolerance = ref["rel_tolerance"]
    print(f"\nBonded cantilever tip deflection, reference {expected:.6f} mm (Timoshenko)")
    print(
        f"{'h [mm]':>8} {'nodes':>8} {'bonded':>10} {'vs ref':>9} {'single':>10} {'vs single':>10}"
    )
    for size, n_nodes, bonded, single in rows:
        print(
            f"{size:>8g} {n_nodes:>8d} {bonded:>10.6f} {bonded / expected - 1:>9.3%} "
            f"{single:>10.6f} {bonded / single - 1:>10.4%}"
        )

    node_counts = [r[1] for r in rows]
    assert node_counts == sorted(node_counts) and len(set(node_counts)) == len(node_counts), (
        "mesh_size is not refining the mesh"
    )
    errors = [abs(r[2] - expected) / expected for r in rows]
    for coarse, fine in itertools.pairwise(errors):
        # As TS-01: the error must not grow under refinement unless it is already well inside.
        assert fine <= coarse + CONVERGENCE_SLACK or fine < tolerance / 4
    assert errors[-1] < tolerance
    for size, _, bonded, single in rows:
        if size <= GATE_MESH_SIZE:
            assert bonded == pytest.approx(single, rel=ref["rel_tolerance_to_single_part"]), (
                f"h = {size:g}: bonded {bonded:.6f} mm vs single part {single:.6f} mm"
            )
