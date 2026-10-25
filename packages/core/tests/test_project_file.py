"""The ``.tessera`` container (ADR-003 amendment) and relative paths in ``*.tessera.json``."""

import json
import os
import zipfile
from pathlib import Path

import pytest

from tessera.core import load_project, new_project, save_project
from tessera.model import FaceSet, InputError, PartInstance, Project, Support


def _project_with(geometry: Path, *, pid: str = "p1") -> Project:
    project = new_project("bracket")
    project.assembly.parts.append(PartInstance(id=pid, name=pid, geometry=geometry))
    project.study.supports.append(Support(id="s1", faces=FaceSet(part_id=pid, faces=(1,))))
    return project


def _write_zip(path: Path, entries: dict[str, bytes | str]) -> Path:
    with zipfile.ZipFile(path, "w") as zf:
        for name, data in entries.items():
            zf.writestr(name, data)
    return path


def _project_json(geometry: str) -> str:
    return json.dumps(
        {"name": "x", "assembly": {"parts": [{"id": "p", "name": "p", "geometry": geometry}]}}
    )


def test_container_holds_project_json_and_a_copy_of_the_geometry(tmp_path: Path) -> None:
    source = tmp_path / "src" / "Bracket v2.STEP"
    source.parent.mkdir()
    source.write_bytes(b"ISO-10303-21; fake step")
    path = tmp_path / "bracket.tessera"

    save_project(_project_with(source), path)

    with zipfile.ZipFile(path) as zf:
        stored = json.loads(zf.read("project.json"))
        entry = stored["assembly"]["parts"][0]["geometry"]
        assert entry.startswith("geometry/")
        assert zf.read(entry) == source.read_bytes()
    assert str(tmp_path) not in json.dumps(stored), "the sender's folders must not leak"


def test_container_opens_after_the_original_geometry_is_gone(tmp_path: Path) -> None:
    source = tmp_path / "part.stl"
    source.write_bytes(b"solid fake")
    path = tmp_path / "part.tessera"
    original = _project_with(source)
    save_project(original, path)
    source.unlink()

    loaded = load_project(path)

    geometry = loaded.assembly.parts[0].geometry
    assert geometry.is_absolute() and geometry.read_bytes() == b"solid fake"
    assert geometry.suffix == ".stl"
    assert loaded.study == original.study


def test_container_round_trips_through_a_second_save(tmp_path: Path) -> None:
    source = tmp_path / "part.step"
    source.write_bytes(b"step bytes")
    first, second = tmp_path / "a.tessera", tmp_path / "b.tessera"
    save_project(_project_with(source), first)

    save_project(load_project(first), second)

    assert load_project(second).assembly.parts[0].geometry.read_bytes() == b"step bytes"


def test_saving_over_the_file_it_was_loaded_from_keeps_the_geometry(tmp_path: Path) -> None:
    source = tmp_path / "part.step"
    source.write_bytes(b"step bytes")
    path = tmp_path / "a.tessera"
    save_project(_project_with(source), path)
    source.unlink()

    project = load_project(path)
    project.name = "renamed"
    save_project(project, path)

    again = load_project(path)
    assert again.name == "renamed"
    assert again.assembly.parts[0].geometry.read_bytes() == b"step bytes"


def test_parts_sharing_a_source_file_share_one_entry(tmp_path: Path) -> None:
    source = tmp_path / "bolt.step"
    source.write_bytes(b"bolt")
    project = _project_with(source)
    project.assembly.parts.append(PartInstance(id="p2", name="p2", geometry=source))
    path = tmp_path / "pair.tessera"

    save_project(project, path)

    with zipfile.ZipFile(path) as zf:
        assert len([n for n in zf.namelist() if n.startswith("geometry/")]) == 1


def test_part_ids_that_are_not_file_names_still_save(tmp_path: Path) -> None:
    source = tmp_path / "a.step"
    source.write_bytes(b"a")
    other = tmp_path / "b.step"
    other.write_bytes(b"b")
    project = _project_with(source, pid="../../etc/passwd")
    project.assembly.parts.append(PartInstance(id="..\\x:y", name="q", geometry=other))
    path = tmp_path / "odd.tessera"

    save_project(project, path)

    with zipfile.ZipFile(path) as zf:
        names = [n for n in zf.namelist() if n.startswith("geometry/")]
    assert len(names) == 2
    assert all(n.count("/") == 1 and ".." not in n and "\\" not in n for n in names)
    loaded = load_project(path)
    assert [p.geometry.read_bytes() for p in loaded.assembly.parts] == [b"a", b"b"]


def test_saving_with_missing_geometry_is_an_error_and_leaves_the_old_file(tmp_path: Path) -> None:
    path = tmp_path / "keep.tessera"
    save_project(new_project("keep"), path)
    before = path.read_bytes()

    with pytest.raises(FileNotFoundError):
        save_project(_project_with(tmp_path / "missing.step"), path)

    assert path.read_bytes() == before
    assert [p.name for p in tmp_path.iterdir()] == ["keep.tessera"], "no temp files left"


def test_container_without_project_json_is_an_input_error(tmp_path: Path) -> None:
    path = _write_zip(tmp_path / "x.tessera", {"readme.txt": "hi"})
    with pytest.raises(InputError, match=r"project\.json"):
        load_project(path)


def test_container_with_invalid_project_json_is_an_input_error(tmp_path: Path) -> None:
    path = _write_zip(tmp_path / "x.tessera", {"project.json": "{"})
    with pytest.raises(InputError, match="not a valid project"):
        load_project(path)


@pytest.mark.parametrize(
    "geometry",
    [
        "../outside.step",
        "/etc/passwd",
        "geometry/../../outside.step",
        "geometry/sub/part.step",
        "geometry/..",
        "C:/part.step",
        "part.step",
    ],
)
def test_container_geometry_must_be_a_plain_entry_in_the_geometry_folder(
    tmp_path: Path, geometry: str
) -> None:
    path = _write_zip(tmp_path / "x.tessera", {"project.json": _project_json(geometry)})
    with pytest.raises(InputError, match="geometry"):
        load_project(path)


def test_container_geometry_entry_must_exist(tmp_path: Path) -> None:
    path = _write_zip(
        tmp_path / "x.tessera", {"project.json": _project_json("geometry/missing.step")}
    )
    with pytest.raises(InputError, match=r"missing\.step"):
        load_project(path)


def test_oversized_entries_are_refused_before_decompression(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    import tessera.core.project_io as project_io

    path = _write_zip(
        tmp_path / "x.tessera",
        {"project.json": _project_json("geometry/p.step"), "geometry/p.step": b"0" * 4096},
    )
    monkeypatch.setattr(project_io, "MAX_ENTRY_BYTES", 1024)
    with pytest.raises(InputError, match="too large"):
        load_project(path)


def test_container_with_newer_schema_is_rejected(tmp_path: Path) -> None:
    path = _write_zip(
        tmp_path / "x.tessera", {"project.json": '{"schema_version": 999, "name": "x"}'}
    )
    with pytest.raises(InputError, match="schema 999"):
        load_project(path)


def test_project_materials_and_raster_angle_round_trip_through_a_container(
    tmp_path: Path,
) -> None:
    from tessera.model import OrthotropicDefinition, PrintSetup

    source = tmp_path / "part.step"
    source.write_bytes(b"step bytes")
    project = _project_with(source)
    project.assembly.parts[0].print_setup = PrintSetup(process="fdm", raster_angle=90.0)
    project.materials.append(
        OrthotropicDefinition(
            id="pla-xyz",
            e1=3000.0,
            e2=2400.0,
            e3=1800.0,
            nu12=0.35,
            nu13=0.30,
            nu23=0.30,
            g12=1000.0,
            g13=700.0,
            g23=650.0,
        )
    )
    path = tmp_path / "part.tessera"
    save_project(project, path)

    with zipfile.ZipFile(path) as zf:
        stored = json.loads(zf.read("project.json"))
    assert stored["schema_version"] == 7
    assert stored["materials"][0]["kind"] == "orthotropic"
    assert stored["assembly"]["parts"][0]["print_setup"]["raster_angle"] == 90.0
    loaded = load_project(path)
    assert loaded.materials == project.materials
    assert loaded.assembly.parts[0].print_setup == project.assembly.parts[0].print_setup


def test_invalid_material_in_a_container_is_an_input_error(tmp_path: Path) -> None:
    bad = {"kind": "orthotropic", "id": "bad", "e1": 3000.0, "e2": 3000.0, "e3": 1800.0}
    bad |= {"nu12": 1.2, "nu13": 0.3, "nu23": 0.3, "g12": 1000.0, "g13": 700.0, "g23": 650.0}
    text = json.dumps({"schema_version": 6, "name": "x", "materials": [bad]})
    path = _write_zip(tmp_path / "x.tessera", {"project.json": text})
    with pytest.raises(InputError, match="positive definite"):
        load_project(path)


def test_json_geometry_paths_resolve_against_the_project_folder(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    folder = tmp_path / "proj"
    (folder / "geometry").mkdir(parents=True)
    (folder / "geometry" / "p.step").write_bytes(b"step")
    path = folder / "project.json"
    path.write_text(_project_json("geometry/p.step"), encoding="utf-8")
    monkeypatch.chdir(tmp_path)  # not the project folder

    geometry = load_project(path).assembly.parts[0].geometry

    assert geometry == folder / "geometry" / "p.step"


def test_json_save_stores_paths_relative_to_the_project_folder(tmp_path: Path) -> None:
    (tmp_path / "parts").mkdir()
    source = tmp_path / "parts" / "p.step"
    source.write_bytes(b"step")
    path = tmp_path / "p.tessera.json"

    save_project(_project_with(source), path)

    stored = json.loads(path.read_text(encoding="utf-8"))
    assert stored["assembly"]["parts"][0]["geometry"] == "parts/p.step"  # "/" on every OS
    assert load_project(path).assembly.parts[0].geometry == source


def test_unzipped_container_is_an_editable_json_project(tmp_path: Path) -> None:
    source = tmp_path / "part.step"
    source.write_bytes(b"step bytes")
    path = tmp_path / "a.tessera"
    save_project(_project_with(source), path)
    unzipped = tmp_path / "unzipped"
    with zipfile.ZipFile(path) as zf:
        zf.extractall(unzipped)

    loaded = load_project(unzipped / "project.json")

    assert loaded.assembly.parts[0].geometry.read_bytes() == b"step bytes"


def test_container_entries_use_forward_slashes(tmp_path: Path) -> None:
    source = tmp_path / "part.step"
    source.write_bytes(b"step")
    path = tmp_path / "a.tessera"
    save_project(_project_with(source), path)
    with zipfile.ZipFile(path) as zf:
        stored = json.loads(zf.read("project.json"))
    assert stored["assembly"]["parts"][0]["geometry"] == "geometry/p1.step"


def test_part_ids_differing_only_in_case_keep_their_own_geometry(tmp_path: Path) -> None:
    upper, lower = tmp_path / "u.step", tmp_path / "l.step"
    upper.write_bytes(b"upper")
    lower.write_bytes(b"lower")
    project = _project_with(upper, pid="Bracket")
    project.assembly.parts.append(PartInstance(id="bracket", name="b", geometry=lower))
    path = tmp_path / "case.tessera"
    save_project(project, path)

    with zipfile.ZipFile(path) as zf:
        names = [n.casefold() for n in zf.namelist()]
    assert len(names) == len(set(names))
    loaded = load_project(path)
    assert [p.geometry.read_bytes() for p in loaded.assembly.parts] == [b"upper", b"lower"]


def test_crafted_entry_names_cannot_collide_on_extraction(tmp_path: Path) -> None:
    parts = [
        {"id": "a", "name": "a", "geometry": "geometry/X.step"},
        {"id": "b", "name": "b", "geometry": "geometry/x.step"},
    ]
    path = _write_zip(
        tmp_path / "x.tessera",
        {
            "project.json": json.dumps({"name": "x", "assembly": {"parts": parts}}),
            "geometry/X.step": b"upper",
            "geometry/x.step": b"lower",
        },
    )
    loaded = load_project(path)
    assert [p.geometry.read_bytes() for p in loaded.assembly.parts] == [b"upper", b"lower"]


def test_saving_an_opened_container_as_json_is_refused(tmp_path: Path) -> None:
    source = tmp_path / "part.step"
    source.write_bytes(b"step")
    path = tmp_path / "a.tessera"
    save_project(_project_with(source), path)

    with pytest.raises(InputError, match="unzip"):
        save_project(load_project(path), tmp_path / "a.tessera.json")


def test_only_geometry_files_are_packed(tmp_path: Path) -> None:
    # A crafted JSON project must not be able to pack e.g. an SSH key into a shared file.
    secret = tmp_path / "id_ed25519"
    secret.write_bytes(b"PRIVATE KEY")
    path = tmp_path / "x.tessera"
    with pytest.raises(InputError, match="not a STEP or STL"):
        save_project(_project_with(secret), path)
    assert not path.exists()


@pytest.mark.skipif(os.name != "posix", reason="POSIX permissions")
def test_saving_keeps_the_file_permissions(tmp_path: Path) -> None:
    path = tmp_path / "shared.tessera"
    save_project(new_project("shared"), path)
    path.chmod(0o640)
    save_project(new_project("shared"), path)
    assert path.stat().st_mode & 0o777 == 0o640


def test_damaged_container_is_an_input_error(tmp_path: Path) -> None:
    good = tmp_path / "good.tessera"
    save_project(new_project("x"), good)
    data = bytearray(good.read_bytes())
    data[40:60] = b"\xff" * 20  # corrupt the compressed project.json
    bad = tmp_path / "bad.tessera"
    bad.write_bytes(bytes(data))
    with pytest.raises(InputError):
        load_project(bad)


def test_a_folder_is_not_a_project_file(tmp_path: Path) -> None:
    with pytest.raises(InputError, match="not a file"):
        load_project(tmp_path)


def test_total_geometry_size_is_capped(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    import tessera.core.project_io as project_io

    parts = [
        {"id": "a", "name": "a", "geometry": "geometry/a.step"},
        {"id": "b", "name": "b", "geometry": "geometry/b.step"},
    ]
    path = _write_zip(
        tmp_path / "x.tessera",
        {
            "project.json": json.dumps({"name": "x", "assembly": {"parts": parts}}),
            "geometry/a.step": b"0" * 600,
            "geometry/b.step": b"0" * 600,
        },
    )
    monkeypatch.setattr(project_io, "MAX_TOTAL_BYTES", 1000)
    with pytest.raises(InputError, match="in total"):
        load_project(path)
