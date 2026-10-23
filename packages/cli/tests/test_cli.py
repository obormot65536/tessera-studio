import dataclasses
import re
from pathlib import Path

import pytest
from typer.testing import CliRunner

from tessera.cli.app import app
from tessera.model import SolverError, SolverTimeoutError, StudyResult

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


def _solved(
    peak_flags: tuple[str, ...], clear: bool, n_faces: int = 7, all_flagged: bool = False
) -> StudyResult:
    """A hand-made result: the raw peak at (0, 5, 5) on face 1, flagged as asked."""
    import numpy as np

    from tessera.model import FlagKind, Mesh, NodeFlags, SolveResult, StressFlag, StressSpot

    flags = tuple(
        StressFlag(FlagKind.SUPPORT, support=name)
        if name != "edge"
        else StressFlag(FlagKind.SHARP_INSIDE_EDGE)
        for name in peak_flags
    )
    n = 10
    mesh = Mesh(nodes=np.zeros((n, 3)), elements=np.arange(10).reshape(1, 10))
    peak = StressSpot(value=143.79, node=0, location=(0.0, 5.0, 5.0), face=1, flags=flags)
    clear_spot = StressSpot(value=98.25, node=1, location=(12.5, -5.0, 5.0), face=None)
    faces = tuple(
        StressSpot(value=150.0 - 10 * i, node=i, location=(float(i), 0.0, 0.0), face=i + 1)
        for i in range(n_faces)
    )
    return StudyResult(
        mesh=mesh,
        solution=SolveResult(np.zeros((n, 3)), np.zeros((n, 6)), solver="fake"),
        max_displacement=1.52238,
        max_von_mises=143.79,
        workdir=Path("run"),
        peak_von_mises=peak,
        clear_peak_von_mises=clear_spot if clear else None,
        face_peaks=tuple(
            dataclasses.replace(f, flags=flags) if i == 0 else f for i, f in enumerate(faces)
        ),
        node_flags=NodeFlags(
            masks={flags[0]: np.ones(n, dtype=bool)} if all_flagged else {},
            n_nodes=n,
            distance=5.0,
            edge_length=2.5,
        ),
    )


def _solve_output(monkeypatch: pytest.MonkeyPatch, result: StudyResult, *args: str) -> list[str]:
    import tessera.cli.app as cli

    asked: list[bool] = []

    def fake_solve(
        _project: object,
        workdir: Path | None = None,
        check_mesh: bool = False,
        timeout: float | None = None,
    ):
        asked.append(check_mesh)
        return result

    monkeypatch.setattr(cli, "load_project", lambda _path: None)
    monkeypatch.setattr(cli, "solve", fake_solve)
    out = runner.invoke(app, ["solve", "p.tessera.json", *args])
    assert out.exit_code == 0, out.output
    assert asked == ["--check-mesh" in args]
    return out.output.splitlines()


def test_solve_leads_with_the_clear_peak_then_the_flagged_raw_peak(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    lines = _solve_output(monkeypatch, _solved(("root", "edge"), clear=True))
    assert lines[:5] == [
        "nodes: 10, elements: 1",
        "max displacement: 1.52238 mm",
        "von Mises stress away from supports and sharp inside edges: 98.25 MPa"
        " at (12.5, -5, 5) mm, inside the part",
        "max von Mises stress: 143.79 MPa at (0, 5, 5) mm, face 1",
        "  near support 'root' and near a sharp inside edge (within 5 mm): the stress there"
        " may keep rising as the mesh is refined; compare two mesh sizes (--check-mesh) before"
        " trusting either value",
    ]
    assert lines[5] == "max von Mises stress by face (5 of 7):"
    assert lines[6] == (
        "  face 1: 150 MPa at (0, 0, 0) mm, near support 'root' and near a sharp inside edge"
    )
    assert lines[7:10] == [
        "  face 2: 140 MPa at (1, 0, 0) mm",
        "  face 3: 130 MPa at (2, 0, 0) mm",
        "  face 4: 120 MPa at (3, 0, 0) mm",
    ]
    assert lines[10] == "  face 5: 110 MPa at (4, 0, 0) mm"
    assert lines[11:] == ["solver files: run"]


def test_solve_without_flags_or_clear_nodes(monkeypatch: pytest.MonkeyPatch) -> None:
    lines = _solve_output(monkeypatch, _solved((), clear=True, n_faces=2))
    assert "max von Mises stress: 143.79 MPa at (0, 5, 5) mm, face 1" in lines
    assert not any(line.startswith("  near") for line in lines)
    assert "max von Mises stress by face (2 of 2):" in lines

    lines = _solve_output(monkeypatch, _solved(("root",), clear=False, all_flagged=True))
    assert lines[2] == (
        "von Mises stress away from supports and sharp inside edges: none, every node is"
        " within 5 mm of a marked support or sharp inside edge"
    )

    # Unflagged nodes exist, but none of them has a stress result; and no face has one either.
    lines = _solve_output(monkeypatch, _solved(("root",), clear=False, n_faces=0))
    assert lines[2] == (
        "von Mises stress away from supports and sharp inside edges: none, no node there has a"
        " stress result"
    )
    assert not any("by face" in line for line in lines)


ORTHOTROPIC_NOTE = (
    "note: this part uses an orthotropic (layer-wise) material; von Mises stress shows where the"
    " stress is, but it is not a strength check: printed parts are weaker across layers."
    " Layer-wise failure checks are planned for v0.5."
)


def test_solve_notes_an_orthotropic_part_before_the_stress_lines(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    isotropic = _solved(("root",), clear=True)
    plain = _solve_output(monkeypatch, isotropic)
    assert ORTHOTROPIC_NOTE not in plain
    assert not any("orthotropic" in line for line in plain)

    lines = _solve_output(monkeypatch, dataclasses.replace(isotropic, orthotropic=True))
    assert lines[:3] == [
        "nodes: 10, elements: 1",
        "max displacement: 1.52238 mm",
        ORTHOTROPIC_NOTE,
    ]
    # Everything else is unchanged: the note is the only added line.
    assert lines[3:] == plain[2:]

    checked = _checked(isotropic, peak=120.0, faces={1: 120.0})
    plain = _solve_output(monkeypatch, checked, "--check-mesh")
    lines = _solve_output(
        monkeypatch, dataclasses.replace(checked, orthotropic=True), "--check-mesh"
    )
    assert lines == [*plain[:3], ORTHOTROPIC_NOTE, *plain[3:]]


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


def _checked(
    result: StudyResult, peak: float, faces: dict[int, float], peak_refinement: float = 0.7
) -> StudyResult:
    """``result`` with a two-mesh check: the coarse raw peak and face peaks as given; every
    spot refined 0.7 times except the raw peak, refined ``peak_refinement`` times."""
    from tessera.model import MeshCheck
    from tessera.post import spot_check

    clear = result.clear_peak_von_mises
    return dataclasses.replace(
        result,
        mesh_check=MeshCheck(
            coarse_edge_length=3.54,
            fine_edge_length=2.5,
            coarse_nodes=4,
            threshold=0.05,
            refinement_limit=0.85,
            max_displacement=spot_check(1.52, result.max_displacement, 0.7),
            peak_von_mises=spot_check(peak, result.max_von_mises, peak_refinement),
            clear_peak_von_mises=None if clear is None else spot_check(97.0, clear.value, 0.7),
            face_peaks={f: spot_check(v, 150.0 - 10 * (f - 1), 0.7) for f, v in faces.items()},
        ),
    )


def test_check_mesh_gives_each_number_a_verdict(monkeypatch: pytest.MonkeyPatch) -> None:
    result = _checked(_solved(("root",), clear=True), peak=120.0, faces={1: 120.0, 2: 139.0})
    lines = _solve_output(monkeypatch, result, "--check-mesh")
    assert lines[:6] == [
        "nodes: 10, elements: 1",
        "mesh check: solved with median element edge 3.54 mm (4 nodes), then 2.5 mm; the numbers"
        " below are from the finer mesh, settled = changed by at most 5 %",
        "max displacement: 1.52238 mm (settled, changed 0.16 % between the meshes)",
        "von Mises stress away from supports and sharp inside edges: 98.25 MPa"
        " at (12.5, -5, 5) mm, inside the part (settled, changed 1.27 % between the meshes)",
        "max von Mises stress: 143.79 MPa at (0, 5, 5) mm, face 1",
        "  near support 'root' (within 5 mm): mesh-dependent, changed 16.5 % between the meshes;"
        " do not compare it with a strength",
    ]
    assert lines[7:10] == [
        "  face 1: 150 MPa at (0, 0, 0) mm, near support 'root' (mesh-dependent, changed 20 %"
        " between the meshes)",
        "  face 2: 140 MPa at (1, 0, 0) mm (settled, changed 0.71 % between the meshes)",
        "  face 3: 130 MPa at (2, 0, 0) mm",  # no peak on the coarse mesh: no verdict
    ]


def test_check_mesh_verdict_replaces_the_flag_warning(monkeypatch: pytest.MonkeyPatch) -> None:
    settled = _checked(_solved(("root",), clear=True), peak=143.0, faces={})
    lines = _solve_output(monkeypatch, settled, "--check-mesh")
    assert (
        lines[5]
        == "  near support 'root' (within 5 mm): settled, changed 0.55 % between the meshes"
    )
    assert not any("may keep rising" in line for line in lines)

    unflagged = _checked(_solved((), clear=True), peak=100.0, faces={})
    lines = _solve_output(monkeypatch, unflagged, "--check-mesh")
    assert lines[5] == (
        "  mesh-dependent, changed 30.5 % between the meshes; do not compare it with a strength"
    )

    overshoot = _checked(_solved((), clear=True), peak=400.0, faces={})
    lines = _solve_output(monkeypatch, overshoot, "--check-mesh")
    assert lines[5].startswith("  mesh-dependent, changed 178 % between the meshes")


def test_check_mesh_percentages_tell_apart_the_two_sides_of_the_threshold(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # 143.79 vs 136.66 is 4.96 %, settled; vs 136.54 it is 5.04 %, mesh-dependent.
    for coarse, expected in (
        (136.66, "settled, changed 4.96 %"),
        (136.54, "mesh-dependent, changed 5.04 %"),
    ):
        lines = _solve_output(
            monkeypatch, _checked(_solved((), clear=True), peak=coarse, faces={}), "--check-mesh"
        )
        assert lines[5].startswith(f"  {expected} between the meshes")


def test_a_spot_the_check_could_not_refine_keeps_its_flag_warning(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # The peak barely moved, but the mesh there barely got finer: no verdict, the flag stands.
    result = _checked(_solved(("root",), clear=True), peak=143.0, faces={}, peak_refinement=0.95)
    lines = _solve_output(monkeypatch, result, "--check-mesh")
    assert lines[5] == (
        "  near support 'root' (within 5 mm): the stress there may keep rising as the mesh is"
        " refined; not judged: the finer mesh is only 0.95 times the element size here"
    )
    assert not any("settled" in line for line in lines[4:6])

    unflagged = _checked(_solved((), clear=True), peak=143.0, faces={}, peak_refinement=0.95)
    lines = _solve_output(monkeypatch, unflagged, "--check-mesh")
    assert lines[5] == "  not judged: the finer mesh is only 0.95 times the element size here"


def test_solve_passes_the_timeout_to_core(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    import tessera.cli.app as cli

    asked: list[float | None] = []

    def fake_solve(_project: object, workdir: Path | None = None, **kwargs: object):
        asked.append(kwargs.get("timeout"))  # type: ignore[arg-type]
        raise SolverError("stop here")

    monkeypatch.setattr(cli, "load_project", lambda _path: None)
    monkeypatch.setattr(cli, "solve", fake_solve)
    runner.invoke(app, ["solve", "p.tessera.json", "--timeout", "90"])
    runner.invoke(app, ["solve", "p.tessera.json"])
    assert asked == [90.0, None]  # None: core's default (SOLVE_TIMEOUT)


def test_a_solver_timeout_says_how_to_give_it_longer(monkeypatch: pytest.MonkeyPatch) -> None:
    import tessera.cli.app as cli

    def slow_solve(*_args: object, **_kwargs: object) -> None:
        raise SolverTimeoutError("CalculiX did not finish within 90 s and was stopped", 90.0)

    monkeypatch.setattr(cli, "load_project", lambda _path: None)
    monkeypatch.setattr(cli, "solve", slow_solve)
    result = runner.invoke(app, ["solve", "p.tessera.json", "--timeout", "90"])
    assert result.exit_code == 2
    assert result.output.strip() == (
        "solver failed: CalculiX did not finish within 90 s and was stopped; "
        "give it longer with --timeout SECONDS"
    )


def test_a_timeout_that_is_not_positive_is_a_clean_error(tmp_path: Path) -> None:
    path = tmp_path / "empty.tessera.json"
    runner.invoke(app, ["new", str(path)])
    result = runner.invoke(app, ["solve", str(path), "--timeout", "0"])
    assert result.exit_code == 1
    assert "error: the solver timeout must be a positive number of seconds, got 0" in result.output
    assert "Traceback" not in result.output


def test_solve_help_shows_the_default_timeout() -> None:
    from tessera.core import SOLVE_TIMEOUT

    result = runner.invoke(app, ["solve", "--help"])
    # Rich styles the help on CI (GITHUB_ACTIONS forces colour), splitting words with escapes.
    text = re.sub(r"\x1b\[[0-9;]*m", "", result.output)
    assert "--timeout" in text
    assert f"{SOLVE_TIMEOUT:g}" in text
