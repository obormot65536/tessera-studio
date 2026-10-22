"""``tessera`` command-line entry point. Thin layer over ``tessera.core``."""

import math
import re
from collections.abc import Callable, Iterator
from contextlib import contextmanager
from importlib.metadata import version
from pathlib import Path

import typer

from tessera.core import (
    WALL_ELEMENTS_MIN,
    describe_edges,
    describe_faces,
    edit,
    is_container_name,
    load_project,
    new_project,
    save_project,
    solve,
)
from tessera.model import (
    ForceLoad,
    InputError,
    MeshCheck,
    MomentLoad,
    Project,
    SolverError,
    SpotCheck,
    StressSpot,
    Verdict,
    WallSpan,
)

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
        if s.edges is not None:
            where = f"edges {_edges(s.edges.edges)}"
        else:
            assert s.faces is not None  # a support holds faces or edges
            where = f"faces {_faces(s.faces.faces)}"
        typer.echo(f"  support {s.id}: {s.part_id} {where}, {held}")
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


def _edges(edges: tuple[tuple[int, int], ...]) -> str:
    return ", ".join(f"{a}-{b}" for a, b in edges)


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

support_app = typer.Typer(help="Hold faces or edges in place.", no_args_is_help=True)
app.add_typer(support_app, name="support")


_EDGE = re.compile(r"\s*([0-9]+)\s*-\s*([0-9]+)\s*")


def _edge(text: str) -> tuple[int, int]:
    match = _EDGE.fullmatch(text)
    if match is None:
        raise InputError(f"--edge takes two face numbers like 3-5; got {text!r}")
    return int(match[1]), int(match[2])


_SUPPORT_FACES = typer.Argument(None, help="Face numbers, see `tessera faces`; or use --edge.")
_SUPPORT_EDGES = typer.Option(
    None, "--edge", help="Edge between two faces, e.g. 3-5 (see `tessera edges`); repeatable."
)


@support_app.command("add")
def support_add(
    path: Path,
    part: str,
    faces: list[int] | None = _SUPPORT_FACES,
    edges: list[str] | None = _SUPPORT_EDGES,
    dofs: str = typer.Option("xyz", help="Directions held: xyz is fixed, z alone a roller."),
    support_id: str | None = typer.Option(None, "--id", help="Support id."),
) -> None:
    """Add a support on faces of a part, or a line support on edges (--edge)."""

    def change(p: Project) -> str:
        if faces and edges:
            raise InputError("give face numbers or --edge, not both")
        if edges:
            pairs = [_edge(e) for e in edges]
            support = edit.add_edge_support(p, part, pairs, dofs=_dofs(dofs), support_id=support_id)
        elif faces:
            support = edit.add_support(p, part, faces, dofs=_dofs(dofs), support_id=support_id)
        else:
            raise InputError("give face numbers or --edge A-B")
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


@app.command()
def edges(
    path: Path,
    part: str,
    face: int | None = typer.Option(None, help="Only the edges of this face."),
) -> None:
    """List a part's edges by face pair (faces, type, length, a point, ends)."""
    with _user_errors():
        project = load_project(path)
        found = describe_edges(project, part)
        if face is not None:
            known = [f.index for f in describe_faces(project, part)]
            if face not in known:
                raise InputError(f"part {part!r} has faces {min(known)} to {max(known)}")
            found = [e for e in found if face in e.faces]
    if not found:
        where = f"face {face} of part {part!r}" if face is not None else f"part {part!r}"
        typer.echo(f"{where} has no edges")
    for e in found:
        pair = f"{e.faces[0]}-{e.faces[1]}"
        point = "({:.2f}, {:.2f}, {:.2f})".format(*e.point)
        if e.curves > 1:
            ends = f"{e.curves} curves"
        elif e.start is None or e.end is None:
            ends = "closed"
        else:
            ends = "from ({:.2f}, {:.2f}, {:.2f})".format(*e.start)
            ends += " to ({:.2f}, {:.2f}, {:.2f})".format(*e.end)
        typer.echo(f"{pair:>9}  {e.kind:<10} {e.length:12.3f} mm  at {point}  {ends}")


FACE_PEAKS_SHOWN = 5
"""Faces listed by ``tessera solve`` (ADR-007); the API has them all."""
_CLEAR = "away from supports and sharp inside edges"


def _where(spot: StressSpot) -> str:
    x, y, z = spot.location
    face = "inside the part" if spot.face is None else f"face {spot.face}"
    return f"at ({x:.4g}, {y:.4g}, {z:.4g}) mm, {face}"


def _reasons(spot: StressSpot) -> str:
    return " and ".join(str(flag) for flag in spot.flags)


def _percent(fraction: float) -> str:
    """0.029 %, 4.96 %, 17.4 %, 178 %: enough digits to tell 4.96 from 5.04 at the threshold."""
    percent = 100 * fraction
    if percent < 1:
        return f"{percent:.2g}"
    return f"{percent:.3g}" if percent < 100 else f"{percent:.0f}"


def _size_at_most(size: float) -> str:
    """``size`` rounded down to two significant digits: a mesh size that is small enough."""
    step = 10.0 ** (math.floor(math.log10(size)) - 1)
    return f"{math.floor(size / step + 1e-9) * step:.2g}"


def _elements(count: float) -> str:
    return f"{count:.3g}"


def _change(
    check: SpotCheck,
    spot: StressSpot | None = None,
    minimum: float = WALL_ELEMENTS_MIN,
    mesh_size: float | None = None,
) -> str:
    """How a number moved between the two meshes of ``--check-mesh``; ``spot`` is where it is,
    ``minimum`` the elements a flagged spot needs across its wall to be judged and
    ``mesh_size`` the coarse mesh's element size, for the advice."""
    if check.verdict is Verdict.NOT_REFINED:
        return (
            f"not judged: the finer mesh is only {check.refinement:.2f} times the element size here"
        )
    if check.verdict is Verdict.THIN_WALL:
        line = f"not judged: changed {_percent(check.change)} % between the meshes, but the "
        wall = None if spot is None else spot.wall
        if check.elements_across is None or wall is None:
            return line + "wall could not be measured here"
        across = check.elements_across
        advice = (
            ""
            if mesh_size is None
            else f"; try --mesh-size {_size_at_most(mesh_size * across / minimum)} or less"
        )
        return line + (
            f"coarse mesh had only {_elements(across)} elements across the "
            f"{wall.thickness:.3g} mm wall here (at least {minimum:g} needed){advice}"
        )
    if math.isinf(check.change):
        return f"{check.verdict}, was {check.coarse:.6g} on the coarse mesh"
    return f"{check.verdict}, changed {_percent(check.change)} % between the meshes"


def _thin_wall(
    spot: StressSpot, check: SpotCheck | None, minimum: float
) -> tuple[WallSpan, float] | None:
    """The wall at a flagged spot, and its elements across, when it has too few for the
    two-mesh check to judge it and no verdict says so already: without a check (the elements
    of this mesh), or where the check did not refine (the coarse mesh's elements)."""
    wall = spot.wall
    if not spot.flagged or wall is None:
        return None
    if check is None:
        across = wall.elements
    elif check.verdict is Verdict.NOT_REFINED and check.elements_across is not None:
        across = check.elements_across
    else:
        return None
    return (wall, across) if across < minimum else None


def _mesh_check_header(check: MeshCheck) -> str:
    return (
        f"mesh check: solved with median element edge {check.coarse_edge_length:.3g} mm "
        f"({check.coarse_nodes} nodes), then {check.fine_edge_length:.3g} mm; the numbers below "
        f"are from the finer mesh, settled = changed by at most {_percent(check.threshold)} %"
    )


@app.command("solve")
def solve_cmd(
    path: Path,
    workdir: Path | None = None,
    check_mesh: bool = typer.Option(
        False,
        "--check-mesh",
        help="Solve again with elements 1/sqrt(2) as long and say which numbers settled "
        "(takes 3 to 5 times as long).",
    ),
) -> None:
    """Mesh and solve the project's static study, then print a summary."""
    with _user_errors():
        project = load_project(path)
        result = solve(project, workdir=workdir, check_mesh=check_mesh)
    check = result.mesh_check
    typer.echo(f"nodes: {result.mesh.n_nodes}, elements: {result.mesh.n_elements}")
    if check is not None:
        typer.echo(_mesh_check_header(check))

    minimum = WALL_ELEMENTS_MIN if check is None else check.wall_elements_min
    # The size the user set for the (coarse) mesh; with an automatic size, its median edge.
    mesh_size = project.study.mesh_size or (
        result.node_flags.edge_length if check is None else check.coarse_edge_length
    )

    def verdict(spot_check: SpotCheck | None, spot: StressSpot | None = None) -> str:
        if spot_check is None:
            return ""
        return f" ({_change(spot_check, spot, minimum, mesh_size)})"

    typer.echo(
        f"max displacement: {result.max_displacement:.6g} mm"
        + verdict(check and check.max_displacement)
    )

    clear, peak = result.clear_peak_von_mises, result.peak_von_mises
    d = result.node_flags.distance
    if clear is None and result.node_flags.flagged.all():
        typer.echo(
            f"von Mises stress {_CLEAR}: none, every node is within {d:.3g} mm of a marked "
            "support or sharp inside edge"
        )
    elif clear is None:
        typer.echo(f"von Mises stress {_CLEAR}: none, no node there has a stress result")
    else:
        typer.echo(
            f"von Mises stress {_CLEAR}: {clear.value:.6g} MPa {_where(clear)}"
            + verdict(check and check.clear_peak_von_mises, clear)
        )
    typer.echo(f"max von Mises stress: {result.max_von_mises:.6g} MPa {_where(peak)}")
    near = f"{_reasons(peak)} (within {d:.3g} mm)"
    unjudged = (Verdict.NOT_REFINED, Verdict.THIN_WALL)
    if check is not None and check.peak_von_mises.verdict not in unjudged:
        # The measured verdict replaces the flags' guess (ADR-007).
        line = _change(check.peak_von_mises, peak, minimum, mesh_size)
        if check.peak_von_mises.verdict is Verdict.MESH_DEPENDENT:
            line += "; do not compare it with a strength"
        typer.echo(f"  {near}: {line}" if peak.flagged else f"  {line}")
    elif peak.flagged:
        advice = (
            "compare two mesh sizes (--check-mesh) before trusting either value"
            if check is None
            else _change(check.peak_von_mises, peak, minimum, mesh_size)
        )
        typer.echo(f"  {near}: the stress there may keep rising as the mesh is refined; {advice}")
    elif check is not None:
        typer.echo(f"  {_change(check.peak_von_mises, peak, minimum, mesh_size)}")
    if (thin := _thin_wall(peak, check and check.peak_von_mises, minimum)) is not None:
        wall, across = thin
        where = "" if check is None else " on the coarse mesh"
        typer.echo(
            f"  thin wall: {_elements(across)} elements across the {wall.thickness:.3g} mm wall "
            f"here{where}; use --mesh-size {_size_at_most(mesh_size * across / minimum)} or "
            "less before --check-mesh, or the comparison cannot see the rise"
        )
    shown = result.face_peaks[:FACE_PEAKS_SHOWN]
    if shown:
        typer.echo(f"max von Mises stress by face ({len(shown)} of {len(result.face_peaks)}):")
    for spot in shown:
        x, y, z = spot.location
        line = f"  face {spot.face}: {spot.value:.6g} MPa at ({x:.4g}, {y:.4g}, {z:.4g}) mm"
        if spot.flagged:
            line += f", {_reasons(spot)}"
        face_check = None if check is None or spot.face is None else check.face_peaks.get(spot.face)
        thin = None if check and not face_check else _thin_wall(spot, face_check, minimum)
        if thin is not None:
            line += f", thin wall ({_elements(thin[1])} elements across)"
        line += verdict(face_check, spot)
        typer.echo(line)
    typer.echo(f"solver files: {result.workdir}")


def main() -> None:
    app()
