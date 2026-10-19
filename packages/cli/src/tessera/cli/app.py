"""``tessera`` command-line entry point. Thin layer over ``tessera.core``."""

from importlib.metadata import version
from pathlib import Path

import typer

from tessera.core import load_project, new_project, save_project

app = typer.Typer(
    help="Tessera Studio: stress analysis for 3D-printed parts.", no_args_is_help=True
)


@app.command("version")
def show_version() -> None:
    """Print the installed Tessera version."""
    typer.echo(version("tessera-cli"))


@app.command()
def new(path: Path, name: str | None = None) -> None:
    """Create an empty project file."""
    save_project(new_project(name or path.name.removesuffix(".tessera.json")), path)
    typer.echo(f"created {path}")


@app.command()
def info(path: Path) -> None:
    """Summarize a project file."""
    project = load_project(path)
    asm = project.assembly
    typer.echo(f"{project.name}: {len(asm.parts)} parts, {len(asm.connections)} connections")


def main() -> None:
    app()
