from pathlib import Path

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
