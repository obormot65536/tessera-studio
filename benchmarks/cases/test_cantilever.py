"""Cantilever beam, tip deflection under a tip shear load (v0.1 validation case, gate: 2 %).

Problem (mm, N, MPa per ADR-002)
--------------------------------
Prismatic steel beam along +x: length L = 200, width b = 10 (y), height h = 10 (z),
E = 210 000 MPa, nu = 0.30 (material ``steel-generic``). The root face x = 0 is fully clamped
(all three displacement components fixed); the tip face x = L carries a total force
P = 100 N in -z, spread uniformly over the face. Linear static analysis, quadratic tetrahedra
(C3D10) in CalculiX. The measured quantity is the mean u_z of all mesh nodes on the tip face.

Reference: Timoshenko beam theory
---------------------------------
    delta = P L^3 / (3 E I) + P L / (kappa G A)

with I = b h^3 / 12, A = b h, G = E / (2 (1 + nu)) and Cowper's shear coefficient for a
rectangular section kappa = 10 (1 + nu) / (12 + 11 nu).

- Bending term: S. P. Timoshenko and J. M. Gere, *Mechanics of Materials*, cantilever with an
  end load (any edition; e.g. 4th ed., PWS 1997, Appendix G, case 1).
- Shear term and Timoshenko beam: S. P. Timoshenko, "On the correction for shear of the
  differential equation for transverse vibrations of prismatic bars", *Phil. Mag.* 41 (1921)
  744-746; static end-load deflection in Timoshenko and Goodier, *Theory of Elasticity*,
  3rd ed., McGraw-Hill 1970, section 21 (cantilever loaded at the end).
- kappa: G. R. Cowper, "The shear coefficient in Timoshenko's beam theory", *J. Appl. Mech.*
  33 (1966) 335-340, Table 1 (rectangle).

Numbers: I = 833.33 mm^4, bending 1.523810 mm, shear 0.002914 mm (0.19 %),
delta = 1.526724 mm.

Why these dimensions
--------------------
L/h = 20 makes the beam slender: the shear term is only 0.19 % of the total, so the choice of
shear coefficient (Cowper vs. 5/6) changes the reference by < 0.01 %, and the 3D end effects
that beam theory ignores are confined to about one section depth at each end (Saint-Venant),
i.e. a small fraction of the span. The deflection (L/131) is small, so linear theory holds.

Effect of the fully clamped root face
-------------------------------------
Beam theory clamps only the centroidal line and leaves the root section free to contract
(Poisson) and warp (shear). Fixing every node of the root face suppresses both, so the 3D model
is slightly stiffer than the beam reference. The effect scales with the shear term, i.e. with
(h/L)^2. An independent gmsh + CalculiX C3D10 run of this exact problem converged to
about -0.35 % of the Timoshenko value (-0.64 % at h_mesh = 10, -0.39 % at h_mesh = 2.5),
well inside the 2 % gate. A uniform tip traction instead of the parabolic shear distribution of
the elasticity solution only changes the deflection locally at the tip and averages out in the
mean tip-face u_z.
"""

from __future__ import annotations

import itertools
import math
import shutil
from collections.abc import Callable
from pathlib import Path
from typing import Any

import numpy as np
import pytest

pytestmark = [
    pytest.mark.benchmark,
    pytest.mark.solver,
    pytest.mark.skipif(shutil.which("ccx") is None, reason="CalculiX (ccx) not on PATH"),
]

REFERENCE = "cantilever_tip_deflection"

# Geometry (mm), load (N), material (MPa)
LENGTH = 200.0
WIDTH = 10.0  # b, along y
HEIGHT = 10.0  # h, along z (load direction)
TIP_FORCE = 100.0  # total, applied in -z
YOUNGS_MODULUS = 210_000.0  # steel-generic
POISSON_RATIO = 0.30

GATE_MESH_SIZE = 2.5  # mm: 4 elements through the depth
CONVERGENCE_MESH_SIZES = (10.0, 5.0, 2.5)  # mm, coarse to fine
CONVERGENCE_SLACK = 0.002  # allowed relative-error increase between refinements (noise)


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


def _face_with_normal(faces: list[Any], normal: tuple[float, float, float]) -> int:
    matches = [
        f.index for f in faces if f.normal is not None and np.allclose(f.normal, normal, atol=1e-6)
    ]
    assert len(matches) == 1, f"expected one face with normal {normal}, got {matches}"
    return matches[0]


def _solve_tip_deflection(step_path: Path, mesh_size: float, workdir: Path) -> tuple[float, int]:
    """Run the cantilever through the public API; return (mean tip u_z, number of nodes)."""
    from tessera.core import describe_faces, new_project, solve
    from tessera.model import FaceSet, ForceLoad, PartInstance, Support

    p = new_project("cantilever")
    p.assembly.parts.append(
        PartInstance(id="beam", name="beam", geometry=step_path, material_id="steel-generic")
    )
    faces = describe_faces(p, "beam")
    root = _face_with_normal(faces, (-1.0, 0.0, 0.0))
    tip = _face_with_normal(faces, (1.0, 0.0, 0.0))
    tip_info = next(f for f in faces if f.index == tip)
    assert tip_info.centroid[0] == pytest.approx(LENGTH)
    assert tip_info.area == pytest.approx(WIDTH * HEIGHT)

    p.study.supports.append(Support(id="root", faces=FaceSet(part_id="beam", faces=(root,))))
    p.study.loads.append(
        ForceLoad(
            id="tip", faces=FaceSet(part_id="beam", faces=(tip,)), force=(0.0, 0.0, -TIP_FORCE)
        )
    )
    p.study.mesh_size = mesh_size

    workdir.mkdir(parents=True, exist_ok=True)
    result = solve(p, workdir=workdir)

    nodes = result.mesh.nodes
    on_tip = np.isclose(nodes[:, 0], LENGTH, rtol=0.0, atol=1e-6 * LENGTH)
    assert on_tip.sum() >= 6, "no mesh nodes found on the tip face"
    on_root = np.isclose(nodes[:, 0], 0.0, rtol=0.0, atol=1e-6 * LENGTH)
    assert np.abs(result.solution.displacements[on_root]).max() < 1e-12, "root is not clamped"
    return float(result.solution.displacements[on_tip, 2].mean()), result.mesh.n_nodes


@pytest.fixture(scope="module")
def beam_step(box_step: Callable[..., Path]) -> Path:
    return box_step("cantilever", (0.0, -WIDTH / 2, -HEIGHT / 2), (LENGTH, WIDTH, HEIGHT))


@pytest.fixture(scope="module")
def tip_deflection(
    beam_step: Path, tmp_path_factory: pytest.TempPathFactory
) -> Callable[[float], tuple[float, int]]:
    """Cached ``mesh_size -> (mean tip u_z, n_nodes)`` so the gate and the sweep share runs."""
    cache: dict[float, tuple[float, int]] = {}

    def run(mesh_size: float) -> tuple[float, int]:
        if mesh_size not in cache:
            workdir = tmp_path_factory.mktemp(f"cantilever-h{mesh_size:g}")
            cache[mesh_size] = _solve_tip_deflection(beam_step, mesh_size, workdir)
        return cache[mesh_size]

    return run


def test_reference_matches_formula(load_reference: Callable[[str], dict[str, Any]]) -> None:
    """The stored reference can never drift from the closed-form value."""
    ref = load_reference(REFERENCE)
    assert math.isclose(
        ref["reference_tip_deflection_mm"], timoshenko_tip_deflection(), rel_tol=1e-9
    )
    assert 0.0 < ref["rel_tolerance"] <= 0.02  # spec gate: within 2 %


def test_cantilever_tip_deflection(
    tip_deflection: Callable[[float], tuple[float, int]],
    load_reference: Callable[[str], dict[str, Any]],
) -> None:
    # Fail (not skip) while the public API is missing, before the reference check can skip.
    from tessera.core import describe_faces, solve  # noqa: F401

    ref = load_reference(REFERENCE)
    expected = ref["reference_tip_deflection_mm"]

    u_z, _ = tip_deflection(GATE_MESH_SIZE)

    assert u_z < 0.0, "tip must deflect in the load direction (-z)"
    assert -u_z == pytest.approx(expected, rel=ref["rel_tolerance"])


def test_cantilever_mesh_convergence(
    tip_deflection: Callable[[float], tuple[float, int]],
    load_reference: Callable[[str], dict[str, Any]],
) -> None:
    # Fail (not skip) while the public API is missing, before the reference check can skip.
    from tessera.core import describe_faces, solve  # noqa: F401

    ref = load_reference(REFERENCE)
    expected = ref["reference_tip_deflection_mm"]
    tolerance = ref["rel_tolerance"]

    rows = []
    for size in CONVERGENCE_MESH_SIZES:
        u_z, n_nodes = tip_deflection(size)
        rows.append((size, n_nodes, -u_z, abs(-u_z - expected) / expected))

    print(f"\nCantilever tip deflection, reference {expected:.6f} mm (Timoshenko)")
    print(f"{'h [mm]':>8} {'nodes':>8} {'delta [mm]':>12} {'rel. error':>11}")
    for size, n_nodes, delta, error in rows:
        print(f"{size:>8g} {n_nodes:>8d} {delta:>12.6f} {error:>11.4%}")

    node_counts = [r[1] for r in rows]
    assert node_counts == sorted(node_counts) and len(set(node_counts)) == len(node_counts), (
        "mesh_size is not refining the mesh"
    )
    errors = [r[3] for r in rows]
    for coarse, fine in itertools.pairwise(errors):
        # Error must not grow under refinement, unless it already sits inside the gate (the
        # converged value differs from beam theory by ~0.35 % because of the clamped root).
        assert fine <= coarse + CONVERGENCE_SLACK or fine < tolerance / 4
    assert errors[-1] < tolerance
