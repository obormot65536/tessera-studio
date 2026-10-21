"""Checks shared by the core modules: project validity, parts and their geometry files."""

from pydantic import ValidationError

from tessera.model import InputError, PartInstance, Project


def validated(project: Project) -> Project:
    """Check the whole project again: lists edited after it was built skip validation."""
    try:
        return Project.model_validate(project.model_dump())
    except ValidationError as exc:
        raise InputError(f"project {project.name!r} is not valid:\n{exc}") from None


def check_geometry(part: PartInstance) -> PartInstance:
    if not part.geometry.is_file():
        raise FileNotFoundError(f"geometry of part {part.id!r} not found: {part.geometry}")
    return part


def part_named(project: Project, part_id: str) -> PartInstance:
    for part in project.assembly.parts:
        if part.id == part_id:
            return check_geometry(part)
    raise InputError(f"no part {part_id!r} in project {project.name!r}")
