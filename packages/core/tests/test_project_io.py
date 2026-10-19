from pathlib import Path

import pytest

from tessera.core import load_project, new_project, save_project


def test_save_and_load(tmp_path: Path) -> None:
    path = tmp_path / "demo.tessera.json"
    save_project(new_project("demo"), path)
    assert load_project(path).name == "demo"


def test_newer_schema_is_rejected(tmp_path: Path) -> None:
    path = tmp_path / "future.tessera.json"
    path.write_text('{"schema_version": 999, "name": "x"}', encoding="utf-8")
    with pytest.raises(ValueError, match="schema 999"):
        load_project(path)
