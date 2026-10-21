"""Change a project: parts, materials, supports, loads and study settings.

Every function checks its input against the project (part ids, face numbers of the actual
geometry, material ids, ids unique across the project, finite numbers) and raises
``InputError`` before changing anything, so a refused change leaves the project as it was.
Clients build their editing commands on these.
"""

import re
from collections.abc import Callable, Iterable
from pathlib import Path
from typing import Any

from pydantic import BaseModel, ValidationError

import tessera.geometry as geometry
import tessera.materials as materials
from tessera.core._checks import part_named, validated
from tessera.model import (
    FaceSet,
    ForceLoad,
    InputError,
    MomentLoad,
    PartInstance,
    PressureLoad,
    Project,
    Study,
    Support,
)

Vec3 = tuple[float, float, float]


def material_ids() -> list[str]:
    """Ids of the built-in materials, for ``material_id``."""
    return materials.material_ids()


def _check_material(material_id: str) -> None:
    if material_id not in materials.material_ids():
        known = ", ".join(materials.material_ids())
        raise InputError(f"unknown material {material_id!r}; known: {known}")


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


def _face_set(project: Project, part_id: str, faces: Iterable[int]) -> FaceSet:
    chosen = tuple(faces)
    if not chosen:
        raise InputError("choose at least one face")
    known = [f.index for f in geometry.describe_faces(part_named(project, part_id).geometry)]
    bad = [f for f in chosen if f not in known]
    if bad:
        raise InputError(
            f"part {part_id!r} has faces {min(known)} to {max(known)}; there is no face "
            f"{bad[0]} (list them with `tessera faces`)"
        )
    if len(set(chosen)) != len(chosen):
        raise InputError(f"a face is listed twice: {list(chosen)}")
    return FaceSet(part_id=part_id, faces=chosen)


def _check_gravity_densities(project: Project) -> None:
    """Gravity needs every part's material to have a density; say so now, not at solve time."""
    if project.study.gravity is None:
        return
    for part in project.assembly.parts:
        if part.material_id is None:
            continue  # solving refuses a part without a material anyway
        if materials.get_isotropic(part.material_id).density is None:
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
        _check_material(material_id)
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
    _check_material(material_id)
    part_named(project, part_id)

    def change(p: Project) -> None:
        part_named(p, part_id).material_id = material_id

    _commit(project, change)


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
    return _add_load(project, _build(ForceLoad, id=item_id, faces=face_set, force=force))


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
    return _add_load(project, _build(PressureLoad, id=item_id, faces=face_set, pressure=pressure))


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
    return _add_load(project, _build(MomentLoad, id=item_id, faces=face_set, moment=moment))


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
            users = [s.id for s in [*study.supports, *study.loads] if s.faces.part_id == item_id]
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
