from pathlib import Path

import pytest
from typer.testing import CliRunner

from tessera.cli.app import app

runner = CliRunner()


def test_new_then_info(tmp_path: Path) -> None:
    path = tmp_path / "bracket.tessera.json"
    assert runner.invoke(app, ["new", str(path)]).exit_code == 0
    result = runner.invoke(app, ["info", str(path)])
    assert result.exit_code == 0
    assert "bracket: 0 parts, 0 connections" in result.output


def test_user_errors_print_a_message_not_a_traceback(tmp_path: Path) -> None:
    path = tmp_path / "empty.tessera.json"
    runner.invoke(app, ["new", str(path)])
    result = runner.invoke(app, ["faces", str(path), "ghost"])
    assert result.exit_code == 1
    assert "error: no part 'ghost'" in result.output
    assert "Traceback" not in result.output


def test_info_on_a_broken_file_is_a_clean_error(tmp_path: Path) -> None:
    path = tmp_path / "broken.tessera.json"
    path.write_text("{", encoding="utf-8")
    result = runner.invoke(app, ["info", str(path)])
    assert result.exit_code == 1
    assert "not a valid project file" in result.output
    assert "Traceback" not in result.output


def test_solve_without_parts_is_a_clean_error(tmp_path: Path) -> None:
    path = tmp_path / "empty.tessera.json"
    runner.invoke(app, ["new", str(path)])
    result = runner.invoke(app, ["solve", str(path)])
    assert result.exit_code == 1
    assert "single part" in result.output


def test_internal_bugs_are_not_hidden(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    import tessera.cli.app as cli

    def broken(*_args: object, **_kwargs: object) -> None:
        raise ValueError("operands could not be broadcast")  # a bug, not user input

    path = tmp_path / "p.tessera.json"
    runner.invoke(app, ["new", str(path)])
    monkeypatch.setattr(cli, "solve", broken)
    result = runner.invoke(app, ["solve", str(path)])
    assert isinstance(result.exception, ValueError)


def test_new_container_then_info(tmp_path: Path) -> None:
    path = tmp_path / "hinge.tessera"
    assert runner.invoke(app, ["new", str(path)]).exit_code == 0
    assert "hinge: 0 parts" in runner.invoke(app, ["info", str(path)]).output


def test_pack_bundles_the_geometry(tmp_path: Path) -> None:
    import json
    import zipfile

    (tmp_path / "part.step").write_bytes(b"step")
    source = tmp_path / "p.tessera.json"
    part = {"id": "a", "name": "a", "geometry": "part.step"}
    source.write_text(json.dumps({"name": "p", "assembly": {"parts": [part]}}), encoding="utf-8")
    target = tmp_path / "p.tessera"

    result = runner.invoke(app, ["pack", str(source), str(target)])

    assert result.exit_code == 0, result.output
    with zipfile.ZipFile(target) as zf:
        assert zf.read("geometry/a.step") == b"step"


def test_pack_refuses_a_target_that_is_not_a_container(tmp_path: Path) -> None:
    source = tmp_path / "p.tessera.json"
    runner.invoke(app, ["new", str(source)])
    result = runner.invoke(app, ["pack", str(source), str(tmp_path / "out.json")])
    assert result.exit_code == 1
    assert ".tessera" in result.output


def test_pack_with_missing_geometry_is_a_clean_error(tmp_path: Path) -> None:
    import json

    source = tmp_path / "p.tessera.json"
    part = {"id": "a", "name": "a", "geometry": "gone.step"}
    source.write_text(json.dumps({"name": "p", "assembly": {"parts": [part]}}), encoding="utf-8")
    result = runner.invoke(app, ["pack", str(source), str(tmp_path / "p.tessera")])
    assert result.exit_code == 1
    assert "not found" in result.output


def test_missing_geometry_is_a_clear_one_line_error(tmp_path: Path) -> None:
    import json

    path = tmp_path / "p.tessera.json"
    part = {"id": "a", "name": "a", "geometry": "gone.step"}
    path.write_text(json.dumps({"name": "p", "assembly": {"parts": [part]}}), encoding="utf-8")
    result = runner.invoke(app, ["faces", str(path), "a"])
    assert result.exit_code == 1
    assert "geometry of part 'a' not found" in result.output
    assert "Traceback" not in result.output


@pytest.fixture(scope="module")
def box_step(tmp_path_factory: pytest.TempPathFactory) -> Path:
    import gmsh

    path = tmp_path_factory.mktemp("geom") / "bracket.step"
    gmsh.initialize(interruptible=False)
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.occ.addBox(0, 0, 0, 40, 10, 10)
        gmsh.model.occ.synchronize()
        gmsh.write(str(path))
    finally:
        gmsh.finalize()
    return path


def _ok(*args: str) -> str:
    result = runner.invoke(app, list(args))
    assert result.exit_code == 0, result.output
    return result.output


@pytest.mark.parametrize("suffix", [".tessera.json", ".tessera"])
def test_set_up_a_study_from_the_command_line(box_step: Path, tmp_path: Path, suffix: str) -> None:
    from tessera.core import load_project
    from tessera.model import ForceLoad, MomentLoad, PressureLoad

    path = str(tmp_path / f"study{suffix}")
    _ok("new", path)
    assert "added part 'bracket'" in _ok("part", "add", path, str(box_step))
    _ok("part", "material", path, "bracket", "steel-generic")
    assert "support-1" in _ok("support", "add", path, "bracket", "1")
    _ok("support", "add", path, "bracket", "3", "4", "--dofs", "z", "--id", "slide")
    _ok("load", "force", path, "bracket", "2", "--force", "0", "0", "-100")
    _ok("load", "pressure", path, "bracket", "5", "--pressure", "0.5")
    _ok("load", "moment", path, "bracket", "2", "--moment", "-1000", "0", "0", "--id", "twist")
    _ok("study", path, "--gravity", "0", "0", "-9810", "--mesh-size", "2.5")

    project = load_project(Path(path))
    part = project.assembly.parts[0]
    assert (part.id, part.material_id) == ("bracket", "steel-generic")
    assert [(s.id, s.faces.faces, s.dofs) for s in project.study.supports] == [
        ("support-1", (1,), (1, 2, 3)),
        ("slide", (3, 4), (3,)),
    ]
    force, pressure, moment = project.study.loads
    assert isinstance(force, ForceLoad) and force.force == (0.0, 0.0, -100.0)
    assert isinstance(pressure, PressureLoad) and pressure.pressure == 0.5
    assert isinstance(moment, MomentLoad) and moment.moment == (-1000.0, 0.0, 0.0)
    assert project.study.gravity == (0.0, 0.0, -9810.0)
    assert project.study.mesh_size == 2.5

    shown = _ok("info", path)
    for text in ("bracket", "steel-generic", "slide", "faces 3, 4", "holds z", "twist", "-9810"):
        assert text in shown

    _ok("remove", path, "twist")
    _ok("study", path, "--no-gravity", "--auto-mesh")
    project = load_project(Path(path))
    assert [load.id for load in project.study.loads] == ["load-1", "load-2"]
    assert project.study.gravity is None and project.study.mesh_size is None


def test_editing_errors_are_one_line_and_leave_the_file_alone(
    box_step: Path, tmp_path: Path
) -> None:
    path = tmp_path / "p.tessera.json"
    _ok("new", str(path))
    _ok("part", "add", str(path), str(box_step))
    before = path.read_text(encoding="utf-8")
    for args in (
        ["support", "add", str(path), "bracket", "9"],
        ["support", "add", str(path), "bracket", "1", "--dofs", "w"],
        ["part", "material", str(path), "bracket", "unobtainium"],
        ["load", "force", str(path), "ghost", "1", "--force", "1", "0", "0"],
        ["remove", str(path), "nothing"],
    ):
        result = runner.invoke(app, args)
        assert result.exit_code == 1, (args, result.output)
        assert result.output.startswith("error: "), result.output
        assert "Traceback" not in result.output
    assert path.read_text(encoding="utf-8") == before


def test_materials_lists_the_library() -> None:
    shown = _ok("materials")
    assert "steel-generic" in shown and "aluminium-6061-t6" in shown
