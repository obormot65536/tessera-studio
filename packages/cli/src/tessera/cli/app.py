"""``tessera`` command-line entry point. Thin layer over ``tessera.core``."""

from importlib.metadata import version
from pathlib import Path

import typer

from tessera.core import describe_faces, load_project, new_project, save_project, solve

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


@app.command()
def faces(path: Path, part: str) -> None:
    """List a part's faces (index, type, area, centroid, normal)."""
    for f in describe_faces(load_project(path), part):
        normal = "curved" if f.normal is None else "({:+.2f}, {:+.2f}, {:+.2f})".format(*f.normal)
        centroid = "({:.2f}, {:.2f}, {:.2f})".format(*f.centroid)
        typer.echo(f"{f.index:4d}  {f.kind:<10} {f.area:12.3f} mm2  {centroid}  {normal}")


@app.command("solve")
def solve_cmd(path: Path, workdir: Path | None = None) -> None:
    """Mesh and solve the project's static study, then print a summary."""
    result = solve(load_project(path), workdir=workdir)
    typer.echo(f"nodes: {result.mesh.n_nodes}, elements: {result.mesh.n_elements}")
    typer.echo(f"max displacement: {result.max_displacement:.6g} mm")
    typer.echo(f"max von Mises stress: {result.max_von_mises:.6g} MPa")


def main() -> None:
    app()
