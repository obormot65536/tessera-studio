"""``tessera`` command-line entry point. Thin layer over ``tessera.core``."""

from collections.abc import Iterator
from contextlib import contextmanager
from importlib.metadata import version
from pathlib import Path

import typer

from tessera.core import describe_faces, load_project, new_project, save_project, solve
from tessera.model import InputError, SolverError

app = typer.Typer(
    help="Tessera Studio: stress analysis for 3D-printed parts.", no_args_is_help=True
)


@contextmanager
def _user_errors() -> Iterator[None]:
    """Report problems in the user's input as one line instead of a traceback."""
    try:
        yield
    except (InputError, FileNotFoundError, NotImplementedError) as exc:
        typer.echo(f"error: {exc}", err=True)
        raise typer.Exit(1) from None
    except SolverError as exc:
        typer.echo(f"solver failed: {exc}", err=True)
        raise typer.Exit(2) from None


@app.command("version")
def show_version() -> None:
    """Print the installed Tessera version."""
    typer.echo(version("tessera-cli"))


@app.command()
def new(path: Path, name: str | None = None) -> None:
    """Create an empty project: a .tessera file, or editable JSON for any other name."""
    default = path.name.removesuffix(".json").removesuffix(".tessera")
    with _user_errors():
        save_project(new_project(name or default), path)
    typer.echo(f"created {path}")


@app.command()
def pack(source: Path, target: Path) -> None:
    """Bundle a project and copies of its geometry into one shareable .tessera file."""
    if not target.name.endswith(".tessera"):
        typer.echo("error: the packed file name must end in .tessera", err=True)
        raise typer.Exit(1)
    with _user_errors():
        save_project(load_project(source), target)
    typer.echo(f"packed {target}")


@app.command()
def info(path: Path) -> None:
    """Summarize a project file."""
    with _user_errors():
        project = load_project(path)
    asm = project.assembly
    typer.echo(f"{project.name}: {len(asm.parts)} parts, {len(asm.connections)} connections")


@app.command()
def faces(path: Path, part: str) -> None:
    """List a part's faces (index, type, area, centroid, normal)."""
    with _user_errors():
        found = describe_faces(load_project(path), part)
    for f in found:
        normal = "curved" if f.normal is None else "({:+.2f}, {:+.2f}, {:+.2f})".format(*f.normal)
        centroid = "({:.2f}, {:.2f}, {:.2f})".format(*f.centroid)
        typer.echo(f"{f.index:4d}  {f.kind:<10} {f.area:12.3f} mm2  {centroid}  {normal}")


@app.command("solve")
def solve_cmd(path: Path, workdir: Path | None = None) -> None:
    """Mesh and solve the project's static study, then print a summary."""
    with _user_errors():
        result = solve(load_project(path), workdir=workdir)
    typer.echo(f"nodes: {result.mesh.n_nodes}, elements: {result.mesh.n_elements}")
    typer.echo(f"max displacement: {result.max_displacement:.6g} mm")
    typer.echo(f"max von Mises stress: {result.max_von_mises:.6g} MPa")
    typer.echo(f"solver files: {result.workdir}")


def main() -> None:
    app()
