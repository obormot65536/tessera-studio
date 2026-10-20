"""core.describe_faces / core.solve on STL parts, with a fake solver (no CalculiX needed)."""

from pathlib import Path

import gmsh
import numpy as np
import pytest

from tessera.core import describe_faces, solve
from tessera.core.project_io import new_project
from tessera.model import (
    AnalysisKind,
    FaceSet,
    ForceLoad,
    PartInstance,
    SolveJob,
    SolveResult,
    Support,
)


@pytest.fixture(scope="module")
def box_files(tmp_path_factory: pytest.TempPathFactory) -> tuple[Path, Path]:
    d = tmp_path_factory.mktemp("geom")
    step, stl = d / "box.step", d / "box.stl"
    gmsh.initialize(interruptible=False)
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.occ.addBox(0, 0, 0, 40, 10, 10)
        gmsh.model.occ.synchronize()
        gmsh.write(str(step))
        gmsh.model.mesh.generate(2)
        gmsh.write(str(stl))
    finally:
        gmsh.finalize()
    return step, stl


class FakeSolver:
    name = "fake"

    def __init__(self) -> None:
        self.job: SolveJob | None = None

    def capabilities(self) -> frozenset[AnalysisKind]:
        return frozenset({AnalysisKind.STATIC})

    def run(self, job: SolveJob, workdir: Path) -> SolveResult:
        self.job = job
        n = job.mesh.n_nodes
        return SolveResult(
            displacements=np.zeros((n, 3)), stresses=np.zeros((n, 6)), solver=self.name
        )


def _project(geometry: Path, root: int = 1, tip: int = 1):
    p = new_project("box")
    p.assembly.parts.append(
        PartInstance(id="b", name="b", geometry=geometry, material_id="steel-generic")
    )
    p.study.supports.append(Support(id="wall", faces=FaceSet(part_id="b", faces=(root,))))
    p.study.loads.append(
        ForceLoad(id="push", faces=FaceSet(part_id="b", faces=(tip,)), force=(0.0, 0.0, -50.0))
    )
    p.study.mesh_size = 4.0
    return p


def _root_and_tip(geometry: Path) -> tuple[int, int]:
    faces = describe_faces(_project(geometry), "b")
    root = next(f for f in faces if f.normal == pytest.approx((-1.0, 0.0, 0.0)))
    tip = next(f for f in faces if f.normal == pytest.approx((1.0, 0.0, 0.0)))
    return root.index, tip.index


def test_describe_faces_on_an_stl_part(box_files: tuple[Path, Path]) -> None:
    _, stl = box_files
    faces = describe_faces(_project(stl), "b")
    assert len(faces) == 6
    assert all(f.kind == "Plane" for f in faces)
    assert sorted(round(f.area, 6) for f in faces) == pytest.approx([100, 100, 400, 400, 400, 400])


def test_stl_and_step_box_give_the_same_solve_job(
    box_files: tuple[Path, Path], tmp_path: Path
) -> None:
    totals = {}
    for name, geometry in zip(("step", "stl"), box_files, strict=True):
        root, tip = _root_and_tip(geometry)
        fake = FakeSolver()
        solve(_project(geometry, root, tip), workdir=tmp_path / name, adapter=fake)
        job = fake.job
        assert job is not None
        assert [s.element_set for s in job.sections] == ["part:b"]
        support_nodes = job.mesh.node_sets["support:wall"]
        assert support_nodes.size > 0
        assert np.allclose(job.mesh.nodes[support_nodes, 0], 0.0)
        loaded = np.concatenate([load.nodes for load in job.loads])
        assert np.allclose(job.mesh.nodes[loaded, 0], 40.0)
        totals[name] = sum(load.forces.sum(axis=0) for load in job.loads)
    assert totals["stl"] == pytest.approx([0.0, 0.0, -50.0])
    assert totals["stl"] == pytest.approx(totals["step"])
