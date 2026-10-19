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
