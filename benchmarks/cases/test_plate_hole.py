"""Plate with a central hole in tension: stress concentration vs Howland (v0.1 gate: 2 %).

Problem (mm, N, MPa per ADR-002)
--------------------------------
Thin steel strip along x, width W = 50 (y), thickness t = 1 (z), with a central through hole of
diameter d = 10 (d/W = 0.2, t/r = 0.2), E = 210 000 MPa, nu = 0.30 (``steel-generic``). A
uniform gross tension sigma = 100 MPa acts along x at both ends, at x = +-L_half = +-100
(= 2 W from the hole centre). Linear static, quadratic tetrahedra (C3D10) in CalculiX.

Modelled as one eighth (x >= 0, y >= 0, z >= 0) with three symmetry planes, each a roller
support that fixes only the normal displacement: u_x = 0 on x = 0, u_y = 0 on y = 0, u_z = 0 on
z = 0 (the mid-plane of the plate). The end face x = L_half carries the force
sigma * (W/2) * (t/2) = 1250 N in +x, spread uniformly.

Measured quantity
-----------------
sigma_xx (load direction) at the single mesh node at (0, r, 0): the hole edge on the net
section (the plane through the hole centre perpendicular to the load), at mid-thickness. That
point is a vertex of the eighth-model geometry, so a node always sits exactly there, whatever
the mesh size. There sigma_xx = sigma_theta_theta is the peak hoop stress.

    Kt_n = sigma_xx(0, r, 0) / sigma_net,   sigma_net = sigma * W / (W - d) = 125 MPa

i.e. Kt is formed on the NET section, consistent with the reference below.

Reference: finite-width Kirsch problem, Howland's strip solution (PROVISIONAL)
-----------------------------------------------------------------------------
Infinite plate: Kt = 3 (G. Kirsch, "Die Theorie der Elastizitaet und die Beduerfnisse der
Festigkeitslehre", Z. VDI 42 (1898) 797-807). For a strip of finite width the reference is
Howland's solution as reported in S. P. Timoshenko and J. N. Goodier, *Theory of Elasticity*,
3rd ed., McGraw-Hill 1970 (sec. 35), from R. C. J. Howland, "On the stresses in the
neighbourhood of a circular hole in a strip under tension", Phil. Trans. R. Soc. A 229 (1930)
49-86. At d/W = 0.2 the gross-section factor is Kt_g = 3.14, so

    Kt_n = Kt_g * (1 - d/W) = 3.14 * 0.8 = 2.512

PROVISIONAL: the maintainer approved this source; the tabulated 3.14 is still to be confirmed
in the book and may be refined.

Accuracy of the reference at d/W = 0.2: published values span 2.5065-2.5190 (0.5 %).
Heywood, *Designing by Photoelasticity*, Chapman & Hall 1952: Kt_n = 2 + (1 - d/W)^3 = 2.512;
older Peterson cubic 3 - 3.14 a + 3.667 a^2 - 1.527 a^3 = 2.5065; Pilkey & Pilkey, *Peterson's
Stress Concentration Factors*, 3rd ed., Wiley 2008, Chart 4.1 fit = 2.5190. An independent
fine-mesh gmsh + CalculiX plane-stress check (t/r = 0.04, local refinement to r/100)
converged to Kt_n = 2.5217, +0.4 % from Howland.

Why these dimensions
--------------------
- d/W = 0.2: the finite-width correction is small (Kt_g 3.14 vs 3), so the uncertainty of
  the correction is well under 1 % of Kt. At larger d/W the published values diverge more.
- L_half = 2 W: the hole disturbance has decayed long before the loaded end. The same
  independent check gave the same Kt to 0.01 % at L_half = 4 W, and +0.1 % at L_half = W.
- t/r = 0.2 (thin plate, close to plane stress). In 3D, Kt at mid-thickness rises above the
  plane-stress value and Kt at the free surfaces drops (E. Sternberg and M. A. Sadowsky,
  J. Appl. Mech. 16 (1949) 27-38; E. S. Folias and J.-J. Wang, Comput. Mech. 6 (1990)
  379-391). The deviation grows with t/r. Independent converged C3D10 runs of this eighth
  model (Kt_n at mid-thickness / surface):
  t/r = 0.04: 2.5217 / 2.5211; t/r = 0.2: 2.527 / 2.512; t/r = 1: 2.592 / 2.397;
  t/r = 2: 2.631 / 2.270. At t/r = 0.2 the mid-thickness value is +0.2 % from plane stress.

Why symmetry (roller) supports, not fully fixed faces
-----------------------------------------------------
A fully fixed x = 0 face clamps the net section, the very place being measured: it suppresses
the Poisson contraction there and creates a singular corner at the hole edge. The independent
check gave Kt_n = 2.24 at mid-thickness (-11 %) and 3.35 at the surface (+33 %). The other
fully-fixed option is a full plate clamped at one end. It is accurate (2.5215 with local
refinement) but costs 8x the nodes. With the single global mesh size it gave -3.5 % at h = 2
and -0.8 % at h = 1 (145 k nodes, ~1 min), and h = 0.5 (~600 k nodes) is out of budget. Also
no geometry vertex guarantees a node at the peak. Hence this case needs
``Support(..., dofs=...)``.

Independent check of this exact setup (global mesh size, as Tessera meshes it):
h = 2: 2.5043 (-0.58 %), h = 1: 2.5093 (-0.39 %), h = 0.5: 2.5227 (+0.14 %).
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

REFERENCE = "plate_hole_kt"

# Geometry (mm), load (MPa), material (MPa)
WIDTH = 50.0  # W, along y (full plate)
HOLE_DIAMETER = 10.0  # d
THICKNESS = 1.0  # t, along z (full plate)
HALF_LENGTH = 100.0  # hole centre to loaded end, along x (load direction)
GROSS_STRESS = 100.0  # remote tension sigma on the gross section W * t
YOUNGS_MODULUS = 210_000.0  # steel-generic
POISSON_RATIO = 0.30

GATE_MESH_SIZE = 1.0  # mm: about 8 elements on the quarter hole arc
CONVERGENCE_MESH_SIZES = (2.0, 1.0, 0.5)  # mm, coarse to fine
CONVERGENCE_SLACK = 0.002  # allowed relative-error increase between refinements (noise)


def net_stress(
    gross_stress: float = GROSS_STRESS, width: float = WIDTH, diameter: float = HOLE_DIAMETER
) -> float:
    """Nominal net-section stress (MPa) at the hole."""
    return gross_stress * width / (width - diameter)


HOWLAND_KT_GROSS = {0.2: 3.14}
"""Howland's gross-section Kt by d/W, as tabulated in Timoshenko & Goodier (provisional)."""


def finite_width_kt_net(diameter: float = HOLE_DIAMETER, width: float = WIDTH) -> float:
    """Net-section Kt from Howland's tabulated gross-section value (Timoshenko & Goodier)."""
    ratio = round(diameter / width, 6)
    return HOWLAND_KT_GROSS[ratio] * (1.0 - diameter / width)


def write_plate_hole_step(path: Path) -> Path:
    """Eighth of the plate (x, y, z >= 0) minus the hole, as a STEP file (gmsh OpenCASCADE)."""
    import gmsh

    owns_session = not gmsh.isInitialized()
    if owns_session:
        gmsh.initialize(interruptible=False)
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add(f"benchmark-{path.stem}")
    try:
        occ = gmsh.model.occ
        plate = occ.addBox(0.0, 0.0, 0.0, HALF_LENGTH, WIDTH / 2, THICKNESS / 2)
        hole = occ.addCylinder(0.0, 0.0, -1.0, 0.0, 0.0, THICKNESS + 2.0, HOLE_DIAMETER / 2)
        occ.cut([(3, plate)], [(3, hole)])
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


def _require_roller_support() -> None:
    """Fail (not skip) until the public API has symmetry/roller supports (``Support.dofs``)."""
    from tessera.core import describe_faces, solve  # noqa: F401
    from tessera.model import FaceSet, Support

    support = Support(id="probe", faces=FaceSet(part_id="p", faces=(1,)), dofs=(1,))
    assert tuple(support.dofs) == (1,)


def _solve_kt(
    step_path: Path, mesh_size: float, workdir: Path, *, check_mesh: bool = False
) -> tuple[float, float, int, Any]:
    """Run the eighth model through the public API.

    Returns (Kt_n at mid-thickness, Kt_n at the surface, number of nodes, the StudyResult).
    """
    from tessera.core import describe_faces, new_project, solve
    from tessera.model import FaceSet, ForceLoad, PartInstance, Support

    p = new_project("plate-hole")
    p.assembly.parts.append(
        PartInstance(id="plate", name="plate", geometry=step_path, material_id="steel-generic")
    )
    faces = describe_faces(p, "plate")
    sym_x = _face_with_normal(faces, (-1.0, 0.0, 0.0))  # net section, x = 0
    sym_y = _face_with_normal(faces, (0.0, -1.0, 0.0))  # plate centre line, y = 0
    sym_z = _face_with_normal(faces, (0.0, 0.0, -1.0))  # mid-plane, z = 0
    end = _face_with_normal(faces, (1.0, 0.0, 0.0))
    end_info = next(f for f in faces if f.index == end)
    assert end_info.centroid[0] == pytest.approx(HALF_LENGTH)
    assert end_info.area == pytest.approx(WIDTH / 2 * THICKNESS / 2)

    for name, face, dof in (("sym-x", sym_x, 1), ("sym-y", sym_y, 2), ("sym-z", sym_z, 3)):
        p.study.supports.append(
            Support(id=name, faces=FaceSet(part_id="plate", faces=(face,)), dofs=(dof,))
        )
    end_force = GROSS_STRESS * end_info.area
    p.study.loads.append(
        ForceLoad(
            id="tension",
            faces=FaceSet(part_id="plate", faces=(end,)),
            force=(end_force, 0.0, 0.0),
        )
    )
    p.study.mesh_size = mesh_size

    workdir.mkdir(parents=True, exist_ok=True)
    result = solve(p, workdir=workdir, check_mesh=check_mesh)

    nodes = result.mesh.nodes
    u = result.solution.displacements
    tol = 1e-6 * HALF_LENGTH
    on_x0 = np.abs(nodes[:, 0]) < tol
    # Rollers, not clamps: normal displacement zero, in-plane (Poisson) contraction free.
    assert np.abs(u[on_x0, 0]).max() < 1e-12, "x = 0 symmetry plane moves along x"
    assert np.abs(u[on_x0, 1]).max() > 1e-6, "x = 0 face is clamped in y, expected a roller"
    on_y0 = np.abs(nodes[:, 1]) < tol
    assert np.abs(u[on_y0, 1]).max() < 1e-12, "y = 0 symmetry plane moves along y"
    on_z0 = np.abs(nodes[:, 2]) < tol
    assert np.abs(u[on_z0, 2]).max() < 1e-12, "z = 0 symmetry plane moves along z"

    radius = HOLE_DIAMETER / 2
    at_edge = on_x0 & (np.abs(nodes[:, 1] - radius) < tol)
    mid = np.flatnonzero(at_edge & on_z0)
    surface = np.flatnonzero(at_edge & (np.abs(nodes[:, 2] - THICKNESS / 2) < tol))
    assert mid.size == 1, f"expected one node at (0, r, 0), found {mid.size}"
    assert surface.size == 1, f"expected one node at (0, r, t/2), found {surface.size}"

    sigma_xx = result.solution.stresses[:, 0]
    s_net = net_stress()
    return (
        float(sigma_xx[mid[0]] / s_net),
        float(sigma_xx[surface[0]] / s_net),
        result.mesh.n_nodes,
        result,
    )


def _raw_von_mises_peak(stresses: np.ndarray) -> tuple[float, int]:
    """(largest von Mises stress, its node), computed here independently of tessera.post."""
    sxx, syy, szz, sxy, syz, szx = stresses.T
    vm = np.sqrt(
        0.5 * ((sxx - syy) ** 2 + (syy - szz) ** 2 + (szz - sxx) ** 2)
        + 3.0 * (sxy**2 + syz**2 + szx**2)
    )
    node = int(np.nanargmax(vm))
    return float(vm[node]), node


@pytest.fixture(scope="module")
def plate_step(tmp_path_factory: pytest.TempPathFactory) -> Path:
    directory = tmp_path_factory.mktemp("plate-hole-geometry")
    return write_plate_hole_step(directory / "plate_hole_eighth.step")


@pytest.fixture(scope="module")
def plate_runs(
    plate_step: Path, tmp_path_factory: pytest.TempPathFactory
) -> Callable[[float], tuple[float, float, int, Any]]:
    """Cached ``mesh_size -> (Kt_n mid, Kt_n surface, n_nodes, StudyResult)``."""
    cache: dict[float, tuple[float, float, int, Any]] = {}

    def run(mesh_size: float) -> tuple[float, float, int, Any]:
        if mesh_size not in cache:
            workdir = tmp_path_factory.mktemp(f"plate-hole-h{mesh_size:g}")
            cache[mesh_size] = _solve_kt(plate_step, mesh_size, workdir)
        return cache[mesh_size]

    return run


@pytest.fixture(scope="module")
def plate_kt(
    plate_runs: Callable[[float], tuple[float, float, int, Any]],
) -> Callable[[float], tuple[float, float, int]]:
    """``mesh_size -> (Kt_n mid, Kt_n surface, n_nodes)``, shared by gate and sweep."""
    return lambda mesh_size: plate_runs(mesh_size)[:3]


def test_reference_matches_formula(load_reference: Callable[[str], dict[str, Any]]) -> None:
    """The stored reference can never drift from the closed-form value."""
    ref = load_reference(REFERENCE)
    assert math.isclose(ref["reference_kt_net"], finite_width_kt_net(), rel_tol=1e-9)
    assert 0.0 < ref["rel_tolerance"] <= 0.02  # spec gate: within 2 %


def test_plate_hole_kt(
    plate_kt: Callable[[float], tuple[float, float, int]],
    load_reference: Callable[[str], dict[str, Any]],
) -> None:
    # Fail (not skip) while the roller support is missing, before the reference check can skip.
    _require_roller_support()

    ref = load_reference(REFERENCE)
    kt_mid, _, _ = plate_kt(GATE_MESH_SIZE)

    assert kt_mid == pytest.approx(ref["reference_kt_net"], rel=ref["rel_tolerance"])


def test_plate_hole_mesh_convergence(
    plate_kt: Callable[[float], tuple[float, float, int]],
    load_reference: Callable[[str], dict[str, Any]],
) -> None:
    # Fail (not skip) while the roller support is missing, before the reference check can skip.
    _require_roller_support()

    ref = load_reference(REFERENCE)
    expected = ref["reference_kt_net"]
    tolerance = ref["rel_tolerance"]

    rows = []
    for size in CONVERGENCE_MESH_SIZES:
        kt_mid, kt_surface, n_nodes = plate_kt(size)
        rows.append((size, n_nodes, kt_mid, kt_surface, abs(kt_mid - expected) / expected))

    print(
        f"\nPlate with hole, Kt_n reference {expected:.5f}"
        " (Howland via Timoshenko & Goodier, d/W = 0.2)"
    )
    print(f"{'h [mm]':>8} {'nodes':>8} {'Kt_n mid':>10} {'Kt_n surf':>10} {'rel. error':>11}")
    for size, n_nodes, kt_mid, kt_surface, error in rows:
        print(f"{size:>8g} {n_nodes:>8d} {kt_mid:>10.5f} {kt_surface:>10.5f} {error:>11.4%}")

    node_counts = [r[1] for r in rows]
    assert node_counts == sorted(node_counts) and len(set(node_counts)) == len(node_counts), (
        "mesh_size is not refining the mesh"
    )
    errors = [r[4] for r in rows]
    for coarse, fine in itertools.pairwise(errors):
        # Error must not grow under refinement, unless it already sits well inside the gate
        # (the converged 3D mid-thickness value is ~0.3 % above the plane-stress reference).
        assert fine <= coarse + CONVERGENCE_SLACK or fine < tolerance / 4
    assert errors[-1] < tolerance


def test_plate_hole_peak_stress_is_clear_at_the_hole_edge(
    plate_runs: Callable[[float], tuple[float, float, int, Any]],
) -> None:
    """ADR-007: symmetry rollers hold one direction and the plate has no inside edges, so no
    node is flagged; the raw peak is the hole-edge stress and the clear peak is the same spot.
    No reference value: this checks where and why, not how much."""
    _require_roller_support()
    result = plate_runs(GATE_MESH_SIZE)[3]

    assert not result.node_flags.flagged.any()
    peak = result.peak_von_mises
    assert not peak.flagged
    value, node = _raw_von_mises_peak(result.solution.stresses)
    assert (peak.value, peak.node) == (pytest.approx(value), node)
    assert result.max_von_mises == pytest.approx(value)
    assert result.clear_peak_von_mises == peak
    x, y, _ = peak.location
    tol = 1e-6 * HALF_LENGTH
    assert math.hypot(x, y) == pytest.approx(HOLE_DIAMETER / 2, abs=tol), "peak not on the hole"
    assert abs(x) < 0.1 * HOLE_DIAMETER, "peak not on the net section of the hole edge"


def test_plate_hole_check_mesh_settles_the_hole_peak(plate_step: Path, tmp_path: Path) -> None:
    """ADR-007 two-mesh check: the hole edge is a real, finite stress raiser, so its peak
    settles (measured -1.3 % from 2 mm to 2/sqrt(2) mm, +0.22 % from 1 mm), unlike a singular
    peak. Run at the coarsest sweep mesh to keep it cheap.

    The hole edge has under one element across the modelled 0.5 mm there, far below the
    thin-wall minimum of 5, but the peak is not flagged, and only flagged spots are guarded
    (ADR-007 Amendment 2): a smooth stress raiser is resolved by quadratic elements even then.
    """
    from tessera.model import Verdict

    _require_roller_support()
    *_, result = _solve_kt(plate_step, CONVERGENCE_MESH_SIZES[0], tmp_path, check_mesh=True)
    check = result.mesh_check
    assert check is not None
    assert not result.peak_von_mises.flagged
    assert check.peak_von_mises.verdict is Verdict.SETTLED
    # The wall is the modelled half thickness: the z = 0 symmetry plane is a one-direction
    # support, which is never taken as a mirror (it could be a table the plate rests on).
    wall = result.peak_von_mises.wall
    assert wall is not None and wall.thickness == pytest.approx(THICKNESS / 2, rel=0.01)
    across = check.peak_von_mises.elements_across
    assert across is not None and across < check.wall_elements_min
    assert check.clear_peak_von_mises is not None
    assert check.clear_peak_von_mises.verdict is Verdict.SETTLED
    assert check.max_displacement.verdict is Verdict.SETTLED
