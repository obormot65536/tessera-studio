"""Change a project: parts, materials, supports, loads and study settings.

Every function checks its input against the project (part ids, face numbers of the actual
geometry, material ids, ids unique across the project, finite numbers) and raises
``InputError`` before changing anything, so a refused change leaves the project as it was.
Clients build their editing commands on these.
"""

import re
from collections.abc import Callable, Iterable, Mapping
from pathlib import Path
from typing import Any, Literal

from pydantic import BaseModel, TypeAdapter, ValidationError

import tessera.geometry as geometry
import tessera.materials as materials
import tessera.printsetup as printsetup
from tessera.core._checks import part_named, validated
from tessera.model import (
    EdgeSet,
    FaceSet,
    ForceLoad,
    InputError,
    MaterialDefinition,
    MomentLoad,
    PartInstance,
    PressureLoad,
    PrintSetup,
    Project,
    Study,
    Support,
)

Vec3 = tuple[float, float, float]


def material_ids() -> list[str]:
    """Ids of the built-in materials, for ``material_id``. Project ones are in ``materials``."""
    return materials.material_ids()


def _check_material(project: Project, material_id: str) -> None:
    """The id names a project material or a built-in one (``get_material`` raises if not)."""
    materials.get_material(material_id, project.materials)


def _all_ids(project: Project) -> set[str]:
    """Ids of parts, connections, supports and loads: one namespace, so ``remove`` is clear."""
    asm, study = project.assembly, project.study
    return {item.id for item in [*asm.parts, *asm.connections, *study.supports, *study.loads]}


def _new_id(prefix: str, taken: set[str], start: int = 1) -> str:
    n = start
    while f"{prefix}-{n}" in taken:
        n += 1
    return f"{prefix}-{n}"


def _checked_id(item_id: str | None, prefix: str, taken: set[str]) -> str:
    if item_id is None:
        return _new_id(prefix, taken)
    if not item_id.strip():
        raise InputError("an id cannot be empty")
    if item_id in taken:
        raise InputError(f"the id {item_id!r} is already used in this project")
    return item_id


def _build[M: BaseModel](model: type[M], **fields: Any) -> M:
    try:
        return model(**fields)
    except ValidationError as exc:
        raise InputError(f"invalid {model.__name__}:\n{exc}") from None


def _check_faces_exist(project: Project, part_id: str, faces: Iterable[int]) -> None:
    known = [f.index for f in geometry.describe_faces(part_named(project, part_id).geometry)]
    bad = [f for f in faces if f not in known]
    if bad:
        raise InputError(
            f"part {part_id!r} has faces {min(known)} to {max(known)}; there is no face "
            f"{bad[0]} (list them with `tessera faces`)"
        )


def _face_set(project: Project, part_id: str, faces: Iterable[int]) -> FaceSet:
    chosen = tuple(faces)
    if not chosen:
        raise InputError("choose at least one face")
    _check_faces_exist(project, part_id, chosen)
    if len(set(chosen)) != len(chosen):
        raise InputError(f"a face is listed twice: {list(chosen)}")
    return FaceSet(part_id=part_id, faces=chosen)


def _edge_set(project: Project, part_id: str, edges: Iterable[tuple[int, int]]) -> EdgeSet:
    chosen = [tuple(edge) for edge in edges]
    if not chosen:
        raise InputError("choose at least one edge")
    pairs: list[tuple[int, int]] = []
    for edge in chosen:
        if len(edge) != 2 or edge[0] == edge[1]:
            raise InputError(
                f"an edge is named by two different faces, e.g. 3-5; got {'-'.join(map(str, edge))}"
            )
        pair = (min(edge), max(edge))
        if pair in pairs:
            raise InputError(f"edge {pair[0]}-{pair[1]} is listed twice")
        pairs.append(pair)
    _check_faces_exist(project, part_id, [f for pair in pairs for f in pair])
    known = {e.faces for e in geometry.describe_edges(part_named(project, part_id).geometry)}
    for a, b in pairs:
        if (a, b) not in known:
            raise InputError(
                f"faces {a} and {b} of part {part_id!r} do not meet at an edge "
                "(list the edges with `tessera edges`)"
            )
    return EdgeSet(part_id=part_id, edges=tuple(pairs))


def _check_gravity_densities(project: Project) -> None:
    """Gravity needs every part's material to have a density; say so now, not at solve time."""
    if project.study.gravity is None:
        return
    for part in project.assembly.parts:
        if part.material_id is None:
            continue  # solving refuses a part without a material anyway
        if materials.get_material(part.material_id, project.materials).density is None:
            raise InputError(
                f"part {part.id!r}: material {part.material_id!r} has no density, "
                "which gravity needs; choose another material or turn gravity off"
            )


def _commit(project: Project, change: Callable[[Project], None]) -> None:
    """Apply ``change`` to a copy, validate it, then apply it to ``project``."""
    trial = project.model_copy(deep=True)
    change(trial)
    validated(trial)
    _check_gravity_densities(trial)
    change(project)


def add_part(
    project: Project,
    geometry_path: Path,
    *,
    part_id: str | None = None,
    name: str | None = None,
    material_id: str | None = None,
) -> PartInstance:
    """Add a STEP or STL part. Its id defaults to the file name, made unique.

    The file is read once here, so a file Tessera cannot use is refused straight away.
    """
    if not geometry_path.is_file():
        raise FileNotFoundError(f"geometry file not found: {geometry_path}")
    suffix = geometry_path.suffix.lower()
    if suffix not in (*geometry.STEP_SUFFIXES, *geometry.STL_SUFFIXES):
        raise InputError(f"{geometry_path.name!r} is not a STEP or STL file")
    if material_id is not None:
        _check_material(project, material_id)
    taken = _all_ids(project)
    if part_id is None:
        stem = re.sub(r"[^A-Za-z0-9_-]+", "-", geometry_path.stem).strip("-") or "part"
        part_id = stem if stem not in taken else _new_id(stem, taken, start=2)
    else:
        part_id = _checked_id(part_id, "part", taken)
    geometry.describe_faces(geometry_path)  # raises InputError for an unusable file
    part = _build(
        PartInstance,
        id=part_id,
        name=name or part_id,
        geometry=geometry_path.absolute(),
        material_id=material_id,
    )
    _commit(project, lambda p: p.assembly.parts.append(part.model_copy()))
    return project.assembly.parts[-1]


def set_material(project: Project, part_id: str, material_id: str) -> None:
    """Give the part a built-in material or one defined in the project (``add_material``)."""
    _check_material(project, material_id)
    part_named(project, part_id)

    def change(p: Project) -> None:
        part_named(p, part_id).material_id = material_id

    _commit(project, change)


_DEFINITION: TypeAdapter[MaterialDefinition] = TypeAdapter(MaterialDefinition)


def add_material(
    project: Project, definition: MaterialDefinition | Mapping[str, Any]
) -> MaterialDefinition:
    """Define a material in the project (ADR-009 D3), for ``set_material``.

    ``definition`` is an ``IsotropicDefinition``, an ``OrthotropicDefinition``, or their fields
    as a mapping with ``kind``. It is checked again here (orthotropic constants with
    ``check_orthotropic``). The id may not be a built-in material's or one already defined.
    """
    fields = definition if isinstance(definition, Mapping) else definition.model_dump()
    try:
        material = _DEFINITION.validate_python(fields)
    except ValidationError as exc:
        raise InputError(f"invalid material definition:\n{exc}") from None
    if material.id in materials.material_ids():
        raise InputError(f"{material.id!r} is a built-in material; choose another id")
    if any(m.id == material.id for m in project.materials):
        raise InputError(f"material {material.id!r} is already defined in this project")
    _commit(project, lambda p: p.materials.append(material.model_copy()))
    return project.materials[-1]


def set_print_setup(
    project: Project,
    part_id: str,
    process: Literal["fdm", "resin"],
    build_direction: Vec3 = (0.0, 0.0, 1.0),
    raster_angle: float = 0.0,
) -> PrintSetup:
    """How the part is printed (ADR-009 D2).

    ``build_direction`` is the part-coordinate vector pointing up on the printer (material axis
    3); ``raster_angle`` turns the rasters about it, in degrees counter-clockwise looking down
    it. An orthotropic part needs one to be solved; an isotropic part ignores it.
    """
    part_named(project, part_id)
    setup = _build(
        PrintSetup, process=process, build_direction=build_direction, raster_angle=raster_angle
    )
    try:
        printsetup.material_frame(setup.build_direction, setup.raster_angle)
    except ValueError as exc:
        raise InputError(f"part {part_id!r}: {exc}") from None

    def change(p: Project) -> None:
        part_named(p, part_id).print_setup = setup.model_copy()

    _commit(project, change)
    return setup


def add_support(
    project: Project,
    part_id: str,
    faces: Iterable[int],
    *,
    dofs: tuple[int, ...] = (1, 2, 3),
    support_id: str | None = None,
) -> Support:
    """Hold the faces: all directions by default, or only ``dofs`` (1 = x, 2 = y, 3 = z)."""
    item_id = _checked_id(support_id, "support", _all_ids(project))
    support = _build(
        Support, id=item_id, faces=_face_set(project, part_id, faces), dofs=tuple(dofs)
    )
    _commit(project, lambda p: p.study.supports.append(support.model_copy()))
    return project.study.supports[-1]


def add_edge_support(
    project: Project,
    part_id: str,
    edges: Iterable[tuple[int, int]],
    *,
    dofs: tuple[int, ...] = (1, 2, 3),
    support_id: str | None = None,
) -> Support:
    """Hold the edges where pairs of faces meet (a line support, ADR-006).

    Each edge is named by its two faces, e.g. ``(3, 5)``. All directions held makes a pinned
    line; one direction, e.g. ``dofs=(3,)``, a roller line. The part can still rotate about
    the line.
    """
    item_id = _checked_id(support_id, "support", _all_ids(project))
    support = _build(
        Support, id=item_id, edges=_edge_set(project, part_id, edges), dofs=tuple(dofs)
    )
    _commit(project, lambda p: p.study.supports.append(support.model_copy()))
    return project.study.supports[-1]


def _check_nonzero(kind: str, *values: float) -> None:
    """A load of zero does nothing, so it is almost certainly a typo."""
    if not any(values):
        raise InputError(f"a {kind} of zero does nothing; give a non-zero value")


def _add_load[L: ForceLoad | PressureLoad | MomentLoad](project: Project, load: L) -> L:
    _commit(project, lambda p: p.study.loads.append(load.model_copy()))
    return load


def add_force_load(
    project: Project,
    part_id: str,
    faces: Iterable[int],
    force: Vec3,
    *,
    load_id: str | None = None,
) -> ForceLoad:
    """Total force in N, spread uniformly over the faces."""
    item_id = _checked_id(load_id, "load", _all_ids(project))
    face_set = _face_set(project, part_id, faces)
    load = _build(ForceLoad, id=item_id, faces=face_set, force=force)
    _check_nonzero("force", *load.force)
    return _add_load(project, load)


def add_pressure_load(
    project: Project,
    part_id: str,
    faces: Iterable[int],
    pressure: float,
    *,
    load_id: str | None = None,
) -> PressureLoad:
    """Pressure in MPa normal to the faces; positive pushes into the part."""
    item_id = _checked_id(load_id, "load", _all_ids(project))
    face_set = _face_set(project, part_id, faces)
    load = _build(PressureLoad, id=item_id, faces=face_set, pressure=pressure)
    _check_nonzero("pressure", load.pressure)
    return _add_load(project, load)


def add_moment_load(
    project: Project,
    part_id: str,
    faces: Iterable[int],
    moment: Vec3,
    *,
    load_id: str | None = None,
) -> MomentLoad:
    """Moment in N mm about the faces' area centroid."""
    item_id = _checked_id(load_id, "load", _all_ids(project))
    face_set = _face_set(project, part_id, faces)
    load = _build(MomentLoad, id=item_id, faces=face_set, moment=moment)
    _check_nonzero("moment", *load.moment)
    return _add_load(project, load)


def set_gravity(project: Project, acceleration: Vec3 | None) -> None:
    """Gravity in mm/s^2 for every part, e.g. (0, 0, -9810); ``None`` turns it off."""
    gravity = _build(Study, gravity=acceleration).gravity  # checked before it touches anything

    def change(p: Project) -> None:
        p.study.gravity = gravity

    _commit(project, change)


def set_mesh_size(project: Project, size: float | None) -> None:
    """Target element size in mm; ``None`` lets the mesher choose."""
    if size is not None and not size > 0:
        raise InputError(f"the mesh size must be positive, got {size:g} mm")

    def change(p: Project) -> None:
        p.study.mesh_size = size

    _commit(project, change)


def remove(project: Project, item_id: str) -> str:
    """Remove the support, load or part with this id; return which kind it was."""
    study = project.study
    for kind, items in (("support", study.supports), ("load", study.loads)):
        for i, item in enumerate(items):
            if item.id == item_id:
                del items[i]
                return kind
    for i, part in enumerate(project.assembly.parts):
        if part.id == item_id:
            users = [s.id for s in study.supports if s.part_id == item_id]
            users += [load.id for load in study.loads if load.faces.part_id == item_id]
            users += [
                c.id for c in project.assembly.connections if item_id in (c.a.part_id, c.b.part_id)
            ]
            if users:
                raise InputError(
                    f"part {item_id!r} is still used by {', '.join(users)}; remove those first"
                )
            del project.assembly.parts[i]
            return "part"
    raise InputError(f"nothing called {item_id!r} in project {project.name!r}")
