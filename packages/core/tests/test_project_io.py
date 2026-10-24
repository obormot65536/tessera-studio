from pathlib import Path

import pytest

from tessera.core import load_project, new_project, save_project
from tessera.model import SCHEMA_VERSION


def test_save_and_load(tmp_path: Path) -> None:
    path = tmp_path / "demo.tessera.json"
    save_project(new_project("demo"), path)
    assert load_project(path).name == "demo"


def test_newer_schema_is_rejected(tmp_path: Path) -> None:
    path = tmp_path / "future.tessera.json"
    path.write_text('{"schema_version": 999, "name": "x"}', encoding="utf-8")
    with pytest.raises(ValueError, match="schema 999"):
        load_project(path)


def test_schema_7_file_is_refused(tmp_path: Path) -> None:
    path = tmp_path / "future.tessera.json"
    path.write_text('{"schema_version": 7, "name": "x", "materials": []}', encoding="utf-8")
    with pytest.raises(ValueError, match=r"uses project schema 7; .* reads up to 6"):
        load_project(path)


def test_schema_5_file_loads_and_saves_as_schema_6(tmp_path: Path) -> None:
    old = tmp_path / "old.tessera.json"
    old.write_text('{"schema_version": 5, "name": "old"}', encoding="utf-8")
    project = load_project(old)
    assert project.materials == []
    new = tmp_path / "new.tessera.json"
    save_project(project, new)
    assert '"schema_version": 6' in new.read_text(encoding="utf-8")


def test_invalid_project_file_is_an_input_error(tmp_path: Path) -> None:
    from tessera.model import InputError

    path = tmp_path / "broken.tessera.json"
    path.write_text('{"name": 3', encoding="utf-8")
    with pytest.raises(InputError, match="not a valid project"):
        load_project(path)


def test_files_from_schema_1_still_load(tmp_path: Path) -> None:
    path = tmp_path / "old.tessera.json"
    path.write_text('{"schema_version": 1, "name": "old"}', encoding="utf-8")
    assert load_project(path).name == "old"


def test_saved_files_carry_the_current_schema(tmp_path: Path) -> None:
    # From schema 3 on, builds check the version before the fields (see test below).
    path = tmp_path / "new.tessera.json"
    save_project(new_project("new"), path)
    assert f'"schema_version": {SCHEMA_VERSION}' in path.read_text(encoding="utf-8")


def test_newer_schema_with_unknown_fields_reports_the_schema(tmp_path: Path) -> None:
    path = tmp_path / "future.tessera.json"
    path.write_text('{"schema_version": 999, "name": "x", "hologram": true}', encoding="utf-8")
    with pytest.raises(ValueError, match="schema 999"):
        load_project(path)


def test_saving_an_older_project_writes_the_current_schema(tmp_path: Path) -> None:
    old = tmp_path / "old.tessera.json"
    old.write_text('{"schema_version": 2, "name": "old"}', encoding="utf-8")
    new = tmp_path / "new.tessera.json"
    save_project(load_project(old), new)
    assert f'"schema_version": {SCHEMA_VERSION}' in new.read_text(encoding="utf-8")


def test_saving_an_invalid_project_is_an_input_error(tmp_path: Path) -> None:
    from tessera.model import FaceSet, InputError, Support

    project = new_project("broken")
    project.study.supports.append(Support(id="s", faces=FaceSet(part_id="ghost", faces=(1,))))
    with pytest.raises(InputError, match="unknown part 'ghost'"):
        save_project(project, tmp_path / "broken.tessera.json")
    assert not (tmp_path / "broken.tessera.json").exists()
