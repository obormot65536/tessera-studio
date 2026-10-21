"""``tessera`` command-line entry point. Thin layer over ``tessera.core``."""

from collections.abc import Callable, Iterator
from contextlib import contextmanager
from importlib.metadata import version
from pathlib import Path

import typer

from tessera.core import (
    describe_faces,
    edit,
    is_container_name,
    load_project,
    new_project,
    save_project,
    solve,
)
from tessera.model import ForceLoad, InputError, MomentLoad, Project, SolverError, StressSpot

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
    if not is_container_name(target):
        typer.echo("error: the packed file name must end in .tessera", err=True)
        raise typer.Exit(1)
    with _user_errors():
        save_project(load_project(source), target)
    typer.echo(f"packed {target}")


_AXES = "xyz"


def _vec(v: tuple[float, ...]) -> str:
    return "(" + ", ".join(f"{x:g}" for x in v) + ")"


@app.command()
def info(path: Path) -> None:
    """Summarize a project: parts, supports, loads and study settings."""
    with _user_errors():
        project = load_project(path)
    asm, study = project.assembly, project.study
    typer.echo(f"{project.name}: {len(asm.parts)} parts, {len(asm.connections)} connections")
    for part in asm.parts:
        material = part.material_id or "no material"
        typer.echo(f"  part {part.id}: {material}")
    for s in study.supports:
        held = "fixed" if s.dofs == (1, 2, 3) else "holds " + "".join(_AXES[d - 1] for d in s.dofs)
        typer.echo(f"  support {s.id}: {s.faces.part_id} faces {_faces(s.faces.faces)}, {held}")
    for load in study.loads:
        where = f"{load.faces.part_id} faces {_faces(load.faces.faces)}"
        if isinstance(load, ForceLoad):
            what = f"force {_vec(load.force)} N"
        elif isinstance(load, MomentLoad):
            what = f"moment {_vec(load.moment)} N mm"
        else:
            what = f"pressure {load.pressure:g} MPa"
        typer.echo(f"  load {load.id}: {what} on {where}")
    if study.gravity is not None:
        typer.echo(f"  gravity: {_vec(study.gravity)} mm/s2")
    mesh = "automatic" if study.mesh_size is None else f"{study.mesh_size:g} mm"
    typer.echo(f"  mesh size: {mesh}")


def _faces(faces: tuple[int, ...]) -> str:
    return ", ".join(str(f) for f in faces)


def _edit(path: Path, change: Callable[[Project], str]) -> None:
    """Load the project, apply one change, save it back and report what changed."""
    with _user_errors():
        project = load_project(path)
        message = change(project)
        save_project(project, path)
    typer.echo(message)


def _dofs(letters: str) -> tuple[int, ...]:
    letters = letters.lower()
    if not letters or any(c not in _AXES for c in letters) or len(set(letters)) != len(letters):
        raise InputError(f"--dofs takes the directions to hold, e.g. xyz or z; got {letters!r}")
    return tuple(sorted(_AXES.index(c) + 1 for c in letters))


@app.command()
def materials() -> None:
    """List the built-in materials."""
    for material_id in edit.material_ids():
        typer.echo(material_id)


part_app = typer.Typer(help="Add parts and choose their material.", no_args_is_help=True)
app.add_typer(part_app, name="part")


@part_app.command("add")
def part_add(
    path: Path,
    geometry: Path,
    part_id: str | None = typer.Option(None, "--id", help="Part id (default: file name)."),
    name: str | None = typer.Option(None, help="Display name (default: the id)."),
    material: str | None = typer.Option(None, help="Material id, see `tessera materials`."),
) -> None:
    """Add a STEP or STL part to the project."""
    _edit(
        path,
        lambda p: (
            "added part "
            + repr(edit.add_part(p, geometry, part_id=part_id, name=name, material_id=material).id)
        ),
    )


@part_app.command("material")
def part_material(path: Path, part: str, material: str) -> None:
    """Set a part's material (see `tessera materials`)."""

    def change(p: Project) -> str:
        edit.set_material(p, part, material)
        return f"part {part!r} is {material}"

    _edit(path, change)


_FACES = typer.Argument(..., help="Face numbers, see `tessera faces`.")

support_app = typer.Typer(help="Hold faces in place.", no_args_is_help=True)
app.add_typer(support_app, name="support")


@support_app.command("add")
def support_add(
    path: Path,
    part: str,
    faces: list[int] = _FACES,
    dofs: str = typer.Option("xyz", help="Directions held: xyz is fixed, z alone a roller."),
    support_id: str | None = typer.Option(None, "--id", help="Support id."),
) -> None:
    """Add a support on faces of a part."""

    def change(p: Project) -> str:
        support = edit.add_support(p, part, faces, dofs=_dofs(dofs), support_id=support_id)
        return f"added support {support.id!r}"

    _edit(path, change)


load_app = typer.Typer(help="Add loads on faces.", no_args_is_help=True)
app.add_typer(load_app, name="load")
_LOAD_ID = typer.Option(None, "--id", help="Load id.")


@load_app.command("force")
def load_force(
    path: Path,
    part: str,
    faces: list[int] = _FACES,
    force: tuple[float, float, float] = typer.Option(..., help="Total force FX FY FZ in N."),
    load_id: str | None = _LOAD_ID,
) -> None:
    """Total force (N) spread uniformly over the faces."""
    _edit(
        path,
        lambda p: f"added load {edit.add_force_load(p, part, faces, force, load_id=load_id).id!r}",
    )


@load_app.command("pressure")
def load_pressure(
    path: Path,
    part: str,
    faces: list[int] = _FACES,
    pressure: float = typer.Option(..., help="MPa; positive pushes into the part."),
    load_id: str | None = _LOAD_ID,
) -> None:
    """Pressure (MPa) normal to the faces."""
    _edit(
        path,
        lambda p: (
            "added load "
            + repr(edit.add_pressure_load(p, part, faces, pressure, load_id=load_id).id)
        ),
    )


@load_app.command("moment")
def load_moment(
    path: Path,
    part: str,
    faces: list[int] = _FACES,
    moment: tuple[float, float, float] = typer.Option(..., help="MX MY MZ in N mm."),
    load_id: str | None = _LOAD_ID,
) -> None:
    """Moment (N mm) about the faces' centroid, e.g. a torque on a shaft end."""
    _edit(
        path,
        lambda p: (
            "added load " + repr(edit.add_moment_load(p, part, faces, moment, load_id=load_id).id)
        ),
    )


@app.command()
def study(
    path: Path,
    gravity: tuple[float, float, float] | None = typer.Option(
        None, help="Gravity GX GY GZ in mm/s2, e.g. 0 0 -9810."
    ),
    no_gravity: bool = typer.Option(False, "--no-gravity", help="Turn gravity off."),
    mesh_size: float | None = typer.Option(None, help="Target element size in mm."),
    auto_mesh: bool = typer.Option(False, "--auto-mesh", help="Let the mesher choose."),
) -> None:
    """Change study settings: gravity and mesh size."""

    def change(p: Project) -> str:
        if gravity is not None and no_gravity:
            raise InputError("use either --gravity or --no-gravity")
        if mesh_size is not None and auto_mesh:
            raise InputError("use either --mesh-size or --auto-mesh")
        if gravity is None and not no_gravity and mesh_size is None and not auto_mesh:
            raise InputError(
                "nothing to change: give --gravity, --no-gravity, --mesh-size or --auto-mesh"
            )
        if gravity is not None or no_gravity:
            edit.set_gravity(p, None if no_gravity else gravity)
        if mesh_size is not None or auto_mesh:
            edit.set_mesh_size(p, None if auto_mesh else mesh_size)
        return "study updated"

    _edit(path, change)


@app.command()
def remove(path: Path, item: str) -> None:
    """Remove a support, load or part by its id."""
    _edit(path, lambda p: f"removed {edit.remove(p, item)} {item!r}")


@app.command()
def faces(path: Path, part: str) -> None:
    """List a part's faces (index, type, area, centroid, normal)."""
    with _user_errors():
        found = describe_faces(load_project(path), part)
    for f in found:
        normal = "curved" if f.normal is None else "({:+.2f}, {:+.2f}, {:+.2f})".format(*f.normal)
        centroid = "({:.2f}, {:.2f}, {:.2f})".format(*f.centroid)
        typer.echo(f"{f.index:4d}  {f.kind:<10} {f.area:12.3f} mm2  {centroid}  {normal}")


FACE_PEAKS_SHOWN = 5
"""Faces listed by ``tessera solve`` (ADR-007); the API has them all."""
_CLEAR = "away from supports and sharp inside edges"


def _where(spot: StressSpot) -> str:
    x, y, z = spot.location
    face = "inside the part" if spot.face is None else f"face {spot.face}"
    return f"at ({x:.4g}, {y:.4g}, {z:.4g}) mm, {face}"


def _reasons(spot: StressSpot) -> str:
    return " and ".join(str(flag) for flag in spot.flags)


@app.command("solve")
def solve_cmd(path: Path, workdir: Path | None = None) -> None:
    """Mesh and solve the project's static study, then print a summary."""
    with _user_errors():
        result = solve(load_project(path), workdir=workdir)
    typer.echo(f"nodes: {result.mesh.n_nodes}, elements: {result.mesh.n_elements}")
    typer.echo(f"max displacement: {result.max_displacement:.6g} mm")

    clear, peak = result.clear_peak_von_mises, result.peak_von_mises
    d = result.node_flags.distance
    if clear is None:
        typer.echo(f"von Mises stress {_CLEAR}: none, every node is within {d:.3g} mm of one")
    else:
        typer.echo(f"von Mises stress {_CLEAR}: {clear.value:.6g} MPa {_where(clear)}")
    typer.echo(f"max von Mises stress: {result.max_von_mises:.6g} MPa {_where(peak)}")
    if peak.flagged:
        typer.echo(
            f"  {_reasons(peak)} (within {d:.3g} mm): the stress there keeps rising as the "
            "mesh is refined, so this value depends on the mesh size"
        )
    shown = result.face_peaks[:FACE_PEAKS_SHOWN]
    typer.echo(f"max von Mises stress by face ({len(shown)} of {len(result.face_peaks)}):")
    for spot in shown:
        x, y, z = spot.location
        line = f"  face {spot.face}: {spot.value:.6g} MPa at ({x:.4g}, {y:.4g}, {z:.4g}) mm"
        typer.echo(f"{line}, {_reasons(spot)}" if spot.flagged else line)
    typer.echo(f"solver files: {result.workdir}")


def main() -> None:
    app()
