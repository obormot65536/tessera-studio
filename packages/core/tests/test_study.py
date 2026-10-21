"""core.solve wiring, tested with a fake solver so it runs without CalculiX."""

from pathlib import Path

import gmsh
import numpy as np
import pytest

from tessera.core import describe_faces, solve
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
