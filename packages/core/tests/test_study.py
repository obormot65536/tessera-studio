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
