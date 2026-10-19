"""Create, save and load Tessera projects.

v0.1 stores a project as one JSON file (``*.tessera.json``). The final
project file format is ADR-003.
"""

from pathlib import Path

from pydantic import ValidationError

from tessera.model import SCHEMA_VERSION, InputError, Project


def new_project(name: str) -> Project:
    return Project(name=name)


def save_project(project: Project, path: Path) -> None:
    path.write_text(project.model_dump_json(indent=2), encoding="utf-8")


def load_project(path: Path) -> Project:
    try:
        project = Project.model_validate_json(path.read_text(encoding="utf-8"))
    except ValidationError as exc:
        raise InputError(f"{path} is not a valid project file:\n{exc}") from None
    if project.schema_version > SCHEMA_VERSION:
        raise InputError(
            f"{path} uses project schema {project.schema_version}; "
            f"this version of Tessera reads up to {SCHEMA_VERSION}"
        )
    return project
