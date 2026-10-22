"""core.solve wiring, tested with a fake solver so it runs without CalculiX."""

import math
from pathlib import Path

import gmsh
import numpy as np
import pytest

from tessera.core import MESH_CHECK_REFINEMENT, describe_faces, solve
from tessera.core.project_io import new_project
from tessera.mesh import mesh_step
from tessera.model import (
    AnalysisKind,
    FaceSet,
    ForceLoad,
    InputError,
    PartInstance,
    SolveJob,
    SolveResult,
    Support,
    Verdict,
)


@pytest.fixture(scope="module")
def bracket_step(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """An L-shaped bracket with a hole: no two faces share area and centroid."""
    path = tmp_path_factory.mktemp("geom") / "bracket.step"
    gmsh.initialize(interruptible=False)
    gmsh.option.setNumber("General.Terminal", 0)
    a = gmsh.model.occ.addBox(0, 0, 0, 60, 20, 5)
    b = gmsh.model.occ.addBox(0, 0, 0, 5, 20, 40)
    (fused, _) = gmsh.model.occ.fuse([(3, a)], [(3, b)])
    hole = gmsh.model.occ.addCylinder(40, 10, -1, 0, 0, 7, 4)
    gmsh.model.occ.cut(fused, [(3, hole)])
    gmsh.model.occ.synchronize()
    gmsh.write(str(path))
    gmsh.finalize()
    return path


class FakeSolver:
    name = "fake"

    def __init__(self) -> None:
        self.job: SolveJob | None = None

    def capabilities(self) -> frozenset[AnalysisKind]:
        return frozenset({AnalysisKind.STATIC})

    def run(self, job: SolveJob, workdir: Path) -> SolveResult:
        self.job = job
        n = job.mesh.n_nodes
        disp = np.zeros((n, 3))
        disp[:, 2] = -job.mesh.nodes[:, 0] * 1e-3
        return SolveResult(displacements=disp, stresses=np.zeros((n, 6)), solver=self.name)


def _project(step: Path, root: int, tip: int):
    p = new_project("bracket")
    p.assembly.parts.append(
        PartInstance(id="b", name="b", geometry=step, material_id="steel-generic")
    )
    p.study.supports.append(Support(id="wall", faces=FaceSet(part_id="b", faces=(root,))))
    p.study.loads.append(
        ForceLoad(id="push", faces=FaceSet(part_id="b", faces=(tip,)), force=(0.0, 0.0, -50.0))
    )
    p.study.mesh_size = 4.0
    return p


def _root_and_tip(step: Path) -> tuple[int, int]:
    faces = describe_faces(_project(step, 1, 1), "b")
    root = next(f for f in faces if f.normal == pytest.approx((-1.0, 0.0, 0.0)))
    plus_x = [f for f in faces if f.normal == pytest.approx((1.0, 0.0, 0.0))]
    tip = max(plus_x, key=lambda f: f.centroid[0])  # end of the arm, not the leg's inner face
    return root.index, tip.index


def test_describe_faces_matches_mesh_surfaces(bracket_step: Path) -> None:
    # Loads and supports reach the right faces only if both imports number faces alike.
    faces = describe_faces(_project(bracket_step, 1, 1), "b")
    mesh = mesh_step(bracket_step, part_id="b", size=3.0)
    assert {f"face:{f.index}" for f in faces} == set(mesh.surfaces)
    x = mesh.nodes
    for f in faces:
        tris = mesh.surfaces[f"face:{f.index}"]
        cross = np.cross(x[tris[:, 1]] - x[tris[:, 0]], x[tris[:, 2]] - x[tris[:, 0]])
        areas = 0.5 * np.linalg.norm(cross, axis=1)
        centroid = (areas[:, None] * x[tris[:, :3]].mean(axis=1)).sum(axis=0) / areas.sum()
        assert areas.sum() == pytest.approx(f.area, rel=2e-2)
        assert centroid == pytest.approx(f.centroid, abs=0.5)


def test_solve_builds_the_expected_job(bracket_step: Path, tmp_path: Path) -> None:
    root, tip = _root_and_tip(bracket_step)
    fake = FakeSolver()
    result = solve(_project(bracket_step, root, tip), workdir=tmp_path, adapter=fake)
    job = fake.job
    assert job is not None
    assert [s.element_set for s in job.sections] == ["part:b"]
    assert job.materials[0].youngs_modulus == 210_000.0
    support_nodes = job.mesh.node_sets["support:wall"]
    assert np.allclose(job.mesh.nodes[support_nodes, 0], 0.0)
    total = sum(load.forces.sum(axis=0) for load in job.loads)
    assert total == pytest.approx([0.0, 0.0, -50.0])
    loaded = np.concatenate([load.nodes for load in job.loads])
    assert np.allclose(job.mesh.nodes[loaded, 0], 60.0)
    assert result.max_displacement == pytest.approx(0.06)
    assert result.workdir == tmp_path


def test_solve_rejects_unknown_face(bracket_step: Path, tmp_path: Path) -> None:
    with pytest.raises(InputError, match="no face"):
        solve(_project(bracket_step, 999, 1), workdir=tmp_path, adapter=FakeSolver())


def test_solve_rejects_adapter_without_capability(bracket_step: Path, tmp_path: Path) -> None:
    class NoStatic(FakeSolver):
        def capabilities(self) -> frozenset[AnalysisKind]:
            return frozenset()

    root, tip = _root_and_tip(bracket_step)
    with pytest.raises(InputError, match="cannot run"):
        solve(_project(bracket_step, root, tip), workdir=tmp_path, adapter=NoStatic())


def test_unknown_part_and_material_are_input_errors(bracket_step: Path) -> None:
    p = _project(bracket_step, 1, 1)
    with pytest.raises(InputError, match="no part 'ghost'"):
        describe_faces(p, "ghost")
    p.assembly.parts[0].material_id = "unobtainium"
    with pytest.raises(InputError, match="unknown material"):
        solve(p, adapter=FakeSolver())


def test_roller_support_passes_only_its_dofs(bracket_step: Path, tmp_path: Path) -> None:
    root, tip = _root_and_tip(bracket_step)
    project = _project(bracket_step, root, tip)
    project.study.supports[0] = Support(
        id="wall", faces=FaceSet(part_id="b", faces=(root,)), dofs=(1,)
    )
    fake = FakeSolver()
    solve(project, workdir=tmp_path, adapter=fake)
    assert fake.job is not None
    assert [s.dofs for s in fake.job.supports] == [(1,)]


def _job(project) -> SolveJob:
    fake = FakeSolver()
    solve(project, adapter=fake)
    assert fake.job is not None
    return fake.job


def test_pressure_pushes_into_plane_and_curved_faces(bracket_step: Path) -> None:
    from tessera.model import PressureLoad

    root, tip = _root_and_tip(bracket_step)
    faces = describe_faces(_project(bracket_step, root, tip), "b")
    bore = [f.index for f in faces if f.kind == "Cylinder"]
    project = _project(bracket_step, root, tip)
    project.study.loads = [
        PressureLoad(id="tip", faces=FaceSet(part_id="b", faces=(tip,)), pressure=3.0),
        PressureLoad(id="bore", faces=FaceSet(part_id="b", faces=tuple(bore)), pressure=7.0),
    ]
    job = _job(project)
    tip_load, bore_load = job.loads
    # Tip face: x = 60, 20 x 5 mm, outward +x, so 3 MPa pushes in -x with 300 N.
    np.testing.assert_allclose(tip_load.forces.sum(axis=0), [-300.0, 0.0, 0.0], atol=1e-9)
    # A full bore balances itself, while pushing outward from its axis (40, 10) by p * 2 pi r t.
    np.testing.assert_allclose(bore_load.forces.sum(axis=0), 0.0, atol=1e-6)
    away = job.mesh.nodes[bore_load.nodes, :2] - [40.0, 10.0]
    away /= np.linalg.norm(away, axis=1)[:, None]
    outward = np.einsum("ij,ij->", bore_load.forces[:, :2], away)
    assert outward == pytest.approx(7.0 * 2 * np.pi * 4 * 5, rel=1e-3)


def test_gravity_weighs_the_whole_part(bracket_step: Path) -> None:
    root, tip = _root_and_tip(bracket_step)
    project = _project(bracket_step, root, tip)
    project.study.loads = []
    project.study.gravity = (0.0, 0.0, -9810.0)
    (weight,) = _job(project).loads
    volume = 60 * 20 * 5 + 5 * 20 * 40 - 5 * 20 * 5 - np.pi * 4**2 * 5
    expected = -7.85e-9 * volume * 9810.0
    total = weight.forces.sum(axis=0)
    assert total[:2] == pytest.approx([0.0, 0.0], abs=1e-12)
    assert total[2] == pytest.approx(expected, rel=1e-4)


def test_gravity_needs_a_material_density(
    bracket_step: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    import tessera.materials as materials
    from tessera.model import IsotropicMaterial

    weightless = IsotropicMaterial(name="weightless", youngs_modulus=1000.0, poisson_ratio=0.3)
    monkeypatch.setattr(materials, "get_isotropic", lambda _id: weightless)
    root, tip = _root_and_tip(bracket_step)
    project = _project(bracket_step, root, tip)
    project.study.gravity = (0.0, 0.0, -9810.0)
    with pytest.raises(InputError, match="no density"):
        solve(project, adapter=FakeSolver())


def test_moment_load_twists_the_tip_face_without_a_net_force(bracket_step: Path) -> None:
    from tessera.model import MomentLoad

    root, tip = _root_and_tip(bracket_step)
    project = _project(bracket_step, root, tip)
    project.study.loads = [
        MomentLoad(id="twist", faces=FaceSet(part_id="b", faces=(tip,)), moment=(2e4, 0.0, 0.0))
    ]
    job = _job(project)
    (twist,) = job.loads
    np.testing.assert_allclose(twist.forces.sum(axis=0), 0.0, atol=1e-8)
    # Tip face x = 60, 0 <= y <= 20, 0 <= z <= 5: centroid (60, 10, 2.5).
    lever = job.mesh.nodes[twist.nodes] - [60.0, 10.0, 2.5]
    np.testing.assert_allclose(np.cross(lever, twist.forces).sum(axis=0), [2e4, 0, 0], atol=1e-6)


def test_loads_and_supports_on_an_unknown_part_are_refused(bracket_step: Path) -> None:
    # Appending to a built project bypasses pydantic validation; solve must still check.
    root, tip = _root_and_tip(bracket_step)
    project = _project(bracket_step, root, tip)
    project.study.loads.append(
        ForceLoad(id="typo", faces=FaceSet(part_id="ghost", faces=(tip,)), force=(1.0, 0, 0))
    )
    with pytest.raises(InputError, match="unknown part 'ghost'"):
        solve(project, adapter=FakeSolver())
    project.study.loads.pop()
    project.study.supports.append(Support(id="s2", faces=FaceSet(part_id="ghost", faces=(1,))))
    with pytest.raises(InputError, match="unknown part 'ghost'"):
        solve(project, adapter=FakeSolver())


def test_duplicate_part_ids_added_later_are_refused(bracket_step: Path) -> None:
    project = _project(bracket_step, 1, 1)
    project.assembly.parts.append(project.assembly.parts[0].model_copy())
    with pytest.raises(InputError, match="unique"):
        describe_faces(project, "b")


def test_missing_geometry_names_the_part(tmp_path: Path) -> None:
    project = _project(tmp_path / "gone.step", 1, 1)
    with pytest.raises(FileNotFoundError, match="geometry of part 'b' not found"):
        describe_faces(project, "b")
    with pytest.raises(FileNotFoundError, match="geometry of part 'b' not found"):
        solve(project, adapter=FakeSolver())


class StressAlongX(FakeSolver):
    """Von Mises stress equal to x (uniaxial), so the peak sits at the tip of the arm."""

    def run(self, job: SolveJob, workdir: Path) -> SolveResult:
        result = super().run(job, workdir)
        stresses = np.zeros((job.mesh.n_nodes, 6))
        stresses[:, 0] = job.mesh.nodes[:, 0]
        return SolveResult(displacements=result.displacements, stresses=stresses, solver="fake")


def test_solve_flags_the_clamp_and_the_inside_corner_but_not_the_hole(
    bracket_step: Path, tmp_path: Path
) -> None:
    from tessera.model import FlagKind, StressFlag

    root, tip = _root_and_tip(bracket_step)
    result = solve(_project(bracket_step, root, tip), workdir=tmp_path, adapter=StressAlongX())
    flags = result.node_flags
    x = result.mesh.nodes
    wall = StressFlag(FlagKind.SUPPORT, support="wall")
    corner = StressFlag(FlagKind.SHARP_INSIDE_EDGE)
    assert list(flags.masks) == [wall, corner]
    # Distances go to the nodes of the clamped face x = 0, so the band ends near x = d.
    assert flags.masks[wall][x[:, 0] < 0.9 * flags.distance].all()
    assert not flags.masks[wall][x[:, 0] > flags.distance + 1e-9].any()
    # The inside edge where the arm (z <= 5) meets the leg (x <= 5) runs along y at x = z = 5.
    to_corner = np.hypot(x[:, 0] - 5.0, x[:, 2] - 5.0)
    assert flags.masks[corner][to_corner < 0.5 * flags.distance].all()
    assert not flags.masks[corner][to_corner > 1.5 * flags.distance].any()

    faces = describe_faces(_project(bracket_step, root, tip), "b")
    hole = next(f.index for f in faces if f.kind == "Cylinder")
    assert not flags.flagged[result.mesh.node_sets[f"face:{hole}"]].any()

    assert result.max_von_mises == pytest.approx(60.0)
    peak = result.peak_von_mises
    assert peak.value == pytest.approx(60.0)
    assert result.mesh.nodes[peak.node, 0] == pytest.approx(60.0)
    assert peak.location[0] == pytest.approx(60.0) and not peak.flagged
    assert result.clear_peak_von_mises == peak
    assert result.face_peaks[0].value == pytest.approx(60.0)
    assert {p.face for p in result.face_peaks} == {f.index for f in faces}


class RecordingSolver(FakeSolver):
    """Keeps every job and working directory, in the order solved."""

    def __init__(self) -> None:
        super().__init__()
        self.runs: list[tuple[SolveJob, Path]] = []

    def run(self, job: SolveJob, workdir: Path) -> SolveResult:
        self.runs.append((job, workdir))
        return super().run(job, workdir)


def test_solve_checks_nothing_unless_asked(bracket_step: Path, tmp_path: Path) -> None:
    root, tip = _root_and_tip(bracket_step)
    fake = RecordingSolver()
    result = solve(_project(bracket_step, root, tip), workdir=tmp_path, adapter=fake)
    assert len(fake.runs) == 1
    assert result.mesh_check is None


def test_check_mesh_solves_again_on_a_finer_mesh(bracket_step: Path, tmp_path: Path) -> None:
    root, tip = _root_and_tip(bracket_step)
    fake = RecordingSolver()
    project = _project(bracket_step, root, tip)
    result = solve(project, workdir=tmp_path, adapter=fake, check_mesh=True)
    (coarse_job, coarse_dir), (fine_job, fine_dir) = fake.runs
    assert (coarse_dir, fine_dir) == (tmp_path / "coarse", tmp_path)
    assert result.workdir == tmp_path
    # Everything reported comes from the fine solve.
    assert result.mesh is fine_job.mesh
    check = result.mesh_check
    assert check is not None
    assert check.coarse_nodes == coarse_job.mesh.n_nodes < fine_job.mesh.n_nodes
    # Elements 1/sqrt(2) as long: about twice the elements in 3D would be 2.8x; allow for gmsh.
    assert check.fine_edge_length / check.coarse_edge_length == pytest.approx(
        MESH_CHECK_REFINEMENT, rel=0.15
    )
    assert project.study.mesh_size == 4.0, "the project is not changed"
    # The fake bends the bracket linearly in x: the same tip deflection on both meshes.
    assert check.max_displacement.verdict is Verdict.SETTLED


def test_check_mesh_scales_every_element_size(
    bracket_step: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    import tessera.core.study as study

    calls: list[tuple[float | None, float]] = []
    mesh_part = study._mesh_part

    def recording(part: PartInstance, size: float | None, size_factor: float = 1.0):
        calls.append((size, size_factor))
        return mesh_part(part, size, size_factor)

    monkeypatch.setattr(study, "_mesh_part", recording)
    root, tip = _root_and_tip(bracket_step)
    project = _project(bracket_step, root, tip)
    for size in (4.0, None):  # the same settings, scaled, also when gmsh picks the size
        project.study.mesh_size = size
        calls.clear()
        solve(project, workdir=tmp_path, adapter=RecordingSolver(), check_mesh=True)
        assert calls == [(size, 1.0), (size, MESH_CHECK_REFINEMENT)]


class RecordingCalculix(FakeSolver):
    """Stands in for ``tessera.solvers.CalculixAdapter``; records the timeout of each run."""

    def __init__(self, timeout: float | None = None) -> None:
        super().__init__()
        self.timeout = timeout
        self.run_timeouts: list[float | None] = []

    def run(self, job: SolveJob, workdir: Path) -> SolveResult:
        self.run_timeouts.append(self.timeout)
        return super().run(job, workdir)


@pytest.fixture
def recording_calculix(monkeypatch: pytest.MonkeyPatch) -> list[RecordingCalculix]:
    """Every adapter ``solve`` builds, in order."""
    import tessera.solvers

    made: list[RecordingCalculix] = []

    def make(timeout: float | None = None) -> RecordingCalculix:
        made.append(RecordingCalculix(timeout))
        return made[-1]

    monkeypatch.setattr(tessera.solvers, "CalculixAdapter", make)
    return made


def test_the_built_in_solver_gets_a_finite_default_timeout(
    bracket_step: Path, tmp_path: Path, recording_calculix: list[RecordingCalculix]
) -> None:
    # ADR-004: a hung solve must not block the GUI or MCP server forever.
    from tessera.core import SOLVE_TIMEOUT

    root, tip = _root_and_tip(bracket_step)
    solve(_project(bracket_step, root, tip), workdir=tmp_path)
    (adapter,) = recording_calculix
    assert adapter.run_timeouts == [SOLVE_TIMEOUT]
    assert math.isfinite(SOLVE_TIMEOUT) and SOLVE_TIMEOUT > 0


def test_the_timeout_can_be_set_per_solve_and_applies_to_each_mesh_of_the_check(
    bracket_step: Path, tmp_path: Path, recording_calculix: list[RecordingCalculix]
) -> None:
    root, tip = _root_and_tip(bracket_step)
    project = _project(bracket_step, root, tip)
    project.study.mesh_size = 8.0
    solve(project, workdir=tmp_path, check_mesh=True, timeout=42.5)
    (adapter,) = recording_calculix
    assert adapter.run_timeouts == [42.5, 42.5]


@pytest.mark.parametrize("timeout", [0.0, -5.0, math.nan, math.inf, True, False, "60", [60]])
def test_the_timeout_must_be_a_positive_number_of_seconds(
    bracket_step: Path, recording_calculix: list[RecordingCalculix], timeout: object
) -> None:
    with pytest.raises(InputError, match="timeout must be a positive number of seconds"):
        solve(_project(bracket_step, 1, 1), timeout=timeout)  # type: ignore[arg-type]
    assert recording_calculix == []


def test_a_timeout_for_an_adapter_you_pass_belongs_on_that_adapter(bracket_step: Path) -> None:
    with pytest.raises(ValueError, match="set it on the adapter"):
        solve(_project(bracket_step, 1, 1), adapter=FakeSolver(), timeout=10.0)


@pytest.mark.parametrize(("hung_run", "which"), [(1, "coarse"), (2, "fine")])
def test_a_timeout_in_the_mesh_check_names_the_solve(
    bracket_step: Path, tmp_path: Path, hung_run: int, which: str
) -> None:
    from tessera.model import SolverTimeoutError

    class Hangs(FakeSolver):
        runs = 0

        def run(self, job: SolveJob, workdir: Path) -> SolveResult:
            self.runs += 1
            if self.runs == hung_run:
                raise SolverTimeoutError("CalculiX did not finish within 9 s", 9.0, "log")
            return super().run(job, workdir)

    root, tip = _root_and_tip(bracket_step)
    project = _project(bracket_step, root, tip)
    project.study.mesh_size = 8.0
    with pytest.raises(SolverTimeoutError) as caught:
        solve(project, workdir=tmp_path, adapter=Hangs(), check_mesh=True)
    assert str(caught.value) == (
        f"CalculiX did not finish within 9 s (in the {which} solve of the mesh check)"
    )
    assert (caught.value.timeout, caught.value.log) == (9.0, "log")


def test_a_timeout_without_the_mesh_check_keeps_its_message(
    bracket_step: Path, tmp_path: Path
) -> None:
    from tessera.model import SolverTimeoutError

    class Hangs(FakeSolver):
        def run(self, job: SolveJob, workdir: Path) -> SolveResult:
            raise SolverTimeoutError("CalculiX did not finish within 9 s", 9.0)

    root, tip = _root_and_tip(bracket_step)
    with pytest.raises(SolverTimeoutError, match=r"^CalculiX did not finish within 9 s$"):
        solve(_project(bracket_step, root, tip), workdir=tmp_path, adapter=Hangs())
