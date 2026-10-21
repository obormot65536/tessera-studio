"""Editing a project through tessera.core: every change is checked before it is applied."""

from pathlib import Path

import gmsh
import pytest

from tessera.core import edit, new_project
from tessera.model import ForceLoad, InputError, MomentLoad, PressureLoad, Project


@pytest.fixture(scope="module")
def box_step(tmp_path_factory: pytest.TempPathFactory) -> Path:
    path = tmp_path_factory.mktemp("geom") / "box.step"
    gmsh.initialize(interruptible=False)
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.occ.addBox(0, 0, 0, 40, 10, 10)
        gmsh.model.occ.synchronize()
        gmsh.write(str(path))
    finally:
        gmsh.finalize()
    return path


@pytest.fixture
def project(box_step: Path) -> Project:
    p = new_project("box")
    edit.add_part(p, box_step, material_id="steel-generic")
    return p


def test_add_part_names_it_after_the_file_and_checks_the_material(box_step: Path) -> None:
    p = new_project("p")
    part = edit.add_part(p, box_step)
    assert (part.id, part.name, part.material_id) == ("box", "box", None)
    second = edit.add_part(p, box_step, material_id="aluminium-6061-t6")
    assert second.id == "box-2"
    with pytest.raises(InputError, match="unknown material 'unobtainium'"):
        edit.add_part(p, box_step, material_id="unobtainium")
    assert len(p.assembly.parts) == 2, "a refused change leaves the project as it was"


def test_add_part_refuses_missing_or_unsupported_geometry(tmp_path: Path) -> None:
    p = new_project("p")
    with pytest.raises(FileNotFoundError, match="not found"):
        edit.add_part(p, tmp_path / "gone.step")
    (tmp_path / "part.obj").write_text("v 0 0 0\n", encoding="utf-8")
    with pytest.raises(InputError, match="STEP or STL"):
        edit.add_part(p, tmp_path / "part.obj")


def test_set_material(project: Project) -> None:
    edit.set_material(project, "box", "aluminium-6061-t6")
    assert project.assembly.parts[0].material_id == "aluminium-6061-t6"
    with pytest.raises(InputError, match="no part 'ghost'"):
        edit.set_material(project, "ghost", "steel-generic")


def test_supports_check_faces_and_get_unique_ids(project: Project) -> None:
    fixed = edit.add_support(project, "box", [1])
    roller = edit.add_support(project, "box", [2, 3], dofs=(3,))
    assert (fixed.id, fixed.dofs) == ("support-1", (1, 2, 3))
    assert (roller.id, roller.faces.faces, roller.dofs) == ("support-2", (2, 3), (3,))
    with pytest.raises(InputError, match="face 7"):
        edit.add_support(project, "box", [7])  # a box has faces 1..6
    with pytest.raises(InputError, match="already used"):
        edit.add_support(project, "box", [4], support_id="support-1")
    assert len(project.study.supports) == 2


def test_each_load_kind(project: Project) -> None:
    force = edit.add_force_load(project, "box", [2], (0.0, 0.0, -100.0))
    pressure = edit.add_pressure_load(project, "box", [3], 0.5)
    moment = edit.add_moment_load(project, "box", [2], (1e3, 0.0, 0.0), load_id="twist")
    assert isinstance(force, ForceLoad) and force.id == "load-1"
    assert isinstance(pressure, PressureLoad) and pressure.id == "load-2"
    assert isinstance(moment, MomentLoad) and moment.id == "twist"
    assert [load.id for load in project.study.loads] == ["load-1", "load-2", "twist"]


def test_ids_are_unique_across_supports_and_loads(project: Project) -> None:
    edit.add_support(project, "box", [1], support_id="end")
    with pytest.raises(InputError, match="already used"):
        edit.add_force_load(project, "box", [2], (1.0, 0, 0), load_id="end")


def test_study_settings(project: Project) -> None:
    edit.set_gravity(project, (0.0, 0.0, -9810.0))
    edit.set_mesh_size(project, 2.5)
    assert project.study.gravity == (0.0, 0.0, -9810.0)
    assert project.study.mesh_size == 2.5
    edit.set_gravity(project, None)
    edit.set_mesh_size(project, None)
    assert project.study.gravity is None and project.study.mesh_size is None
    with pytest.raises(InputError, match="mesh size"):
        edit.set_mesh_size(project, 0.0)


def test_remove_by_id(project: Project) -> None:
    edit.add_support(project, "box", [1])
    edit.add_force_load(project, "box", [2], (1.0, 0, 0))
    assert edit.remove(project, "load-1") == "load"
    assert edit.remove(project, "support-1") == "support"
    with pytest.raises(InputError, match="nothing called 'ghost'"):
        edit.remove(project, "ghost")


def test_removing_a_part_that_is_still_used_is_refused(project: Project) -> None:
    edit.add_support(project, "box", [1])
    with pytest.raises(InputError, match="support-1"):
        edit.remove(project, "box")
    edit.remove(project, "support-1")
    assert edit.remove(project, "box") == "part"
    assert project.assembly.parts == []


def test_non_finite_numbers_are_refused_and_change_nothing(project: Project) -> None:
    with pytest.raises(InputError):
        edit.set_gravity(project, (0.0, float("nan"), 0.0))
    with pytest.raises(InputError):
        edit.set_mesh_size(project, float("inf"))
    with pytest.raises(InputError):
        edit.add_pressure_load(project, "box", [3], float("inf"))
    with pytest.raises(InputError):
        edit.add_force_load(project, "box", [2], (0.0, 0.0, float("-inf")))
    assert project.study.gravity is None and project.study.mesh_size is None
    assert project.study.loads == []


def test_malformed_gravity_is_refused_without_changing_the_project(project: Project) -> None:
    with pytest.raises(InputError):
        edit.set_gravity(project, (0.0, -9810.0))  # type: ignore[arg-type]
    assert project.study.gravity is None


def test_ids_are_unique_across_parts_supports_and_loads(box_step: Path) -> None:
    p = new_project("p")
    edit.add_part(p, box_step, part_id="support-1")
    support = edit.add_support(p, "support-1", [1])
    assert support.id != "support-1", "generated ids skip names parts already use"
    with pytest.raises(InputError, match="already used"):
        edit.add_force_load(p, "support-1", [2], (1.0, 0, 0), load_id="support-1")
    with pytest.raises(InputError, match="already used"):
        edit.add_part(p, box_step, part_id=support.id)


def test_parts_are_checked_and_stored_with_an_absolute_path(
    box_step: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.chdir(box_step.parent)
    p = new_project("p")
    part = edit.add_part(p, Path(box_step.name))
    assert part.geometry.is_absolute()
    broken = tmp_path / "broken.stl"
    broken.write_text("solid x\nendsolid x\n", encoding="utf-8")
    with pytest.raises(InputError):
        edit.add_part(p, broken)
    assert len(p.assembly.parts) == 1
