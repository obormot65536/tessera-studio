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
    assert [(s.id, s.faces and s.faces.faces, s.dofs) for s in project.study.supports] == [
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
        ["load", "force", str(path), "bracket", "1", "--force", "0", "0", "0"],
        ["load", "pressure", str(path), "bracket", "1", "--pressure", "0"],
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


def test_study_without_options_is_refused(tmp_path: Path) -> None:
    path = tmp_path / "p.tessera.json"
    _ok("new", str(path))
    result = runner.invoke(app, ["study", str(path)])
    assert result.exit_code == 1
    assert "--mesh-size" in result.output


def test_container_refusals_leave_the_file_alone(box_step: Path, tmp_path: Path) -> None:
    path = tmp_path / "p.tessera"
    _ok("new", str(path))
    _ok("part", "add", str(path), str(box_step))
    before = path.read_bytes()
    result = runner.invoke(app, ["support", "add", str(path), "bracket", "9"])
    assert result.exit_code == 1
    assert path.read_bytes() == before


def _listed_edges(path: str, *options: str) -> list[str]:
    return [line.split()[0] for line in _ok("edges", path, "bracket", *options).splitlines()]


def test_edges_lists_face_pairs_and_line_supports_use_them(box_step: Path, tmp_path: Path) -> None:
    from tessera.core import load_project

    path = str(tmp_path / "p.tessera.json")
    _ok("new", path)
    _ok("part", "add", path, str(box_step))
    listed = _listed_edges(path)
    assert len(listed) == 12
    lines = _ok("edges", path, "bracket").splitlines()
    assert all("Line" in line and " mm " in line and "from (" in line for line in lines)
    on_face_1 = _listed_edges(path, "--face", "1")
    assert len(on_face_1) == 4 and all("1" in pair.split("-") for pair in on_face_1)

    first, second = on_face_1[:2]
    added = _ok("support", "add", path, "bracket", "--edge", first, "--edge", second, "--dofs", "z")
    assert "support-1" in added
    support = load_project(Path(path)).study.supports[0]
    assert support.edges is not None and support.dofs == (3,)
    assert [f"{a}-{b}" for a, b in support.edges.edges] == [first, second]
    assert f"support support-1: bracket edges {first}, {second}, holds z" in _ok("info", path)


def test_line_support_errors_are_one_line(box_step: Path, tmp_path: Path) -> None:
    path = tmp_path / "p.tessera.json"
    _ok("new", str(path))
    _ok("part", "add", str(path), str(box_step))
    listed = set(_listed_edges(str(path)))
    apart = next(
        f"{a}-{b}" for a in range(1, 7) for b in range(a + 1, 7) if f"{a}-{b}" not in listed
    )
    before = path.read_text(encoding="utf-8")
    for args, message in (
        (["support", "add", str(path), "bracket", "1", "--edge", "1-3"], "not both"),
        (["support", "add", str(path), "bracket"], "--edge A-B"),
        (["support", "add", str(path), "bracket", "--edge", "1x3"], "like 3-5"),
        (["support", "add", str(path), "bracket", "--edge", "\u00b2-3"], "like 3-5"),
        (["support", "add", str(path), "bracket", "--edge", "1-\u0663"], "like 3-5"),
        (["support", "add", str(path), "bracket", "--edge", "2-2"], "two different faces"),
        (["support", "add", str(path), "bracket", "--edge", apart], "tessera edges"),
        (["support", "add", str(path), "bracket", "--edge", "1-9"], "no face 9"),
        (["edges", str(path), "bracket", "--face", "9"], "faces 1 to 6"),
    ):
        result = runner.invoke(app, args)
        assert result.exit_code == 1, (args, result.output)
        assert result.output.startswith("error: ") and message in result.output, result.output
        assert len(result.output.strip().splitlines()) == 1, result.output
    assert path.read_text(encoding="utf-8") == before


def test_edges_of_a_part_without_edges_says_so(tmp_path: Path) -> None:
    import gmsh

    sphere = tmp_path / "ball.step"
    gmsh.initialize(interruptible=False)
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.occ.addSphere(0, 0, 0, 5)
        gmsh.model.occ.synchronize()
        gmsh.write(str(sphere))
    finally:
        gmsh.finalize()
    path = str(tmp_path / "p.tessera.json")
    _ok("new", path)
    _ok("part", "add", path, str(sphere))
    assert _ok("edges", path, "ball").strip() == "part 'ball' has no edges"
    assert _ok("edges", path, "ball", "--face", "1").strip() == "face 1 of part 'ball' has no edges"
