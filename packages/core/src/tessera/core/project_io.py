"""Create, save and load Tessera projects.

v0.1 stores a project as one JSON file (``*.tessera.json``). The final
project file format is ADR-003.
"""

from pathlib import Path

from tessera.model import SCHEMA_VERSION, Project


def new_project(name: str) -> Project:
    return Project(name=name)


def save_project(project: Project, path: Path) -> None:
    path.write_text(project.model_dump_json(indent=2), encoding="utf-8")


def load_project(path: Path) -> Project:
    project = Project.model_validate_json(path.read_text(encoding="utf-8"))
    if project.schema_version > SCHEMA_VERSION:
        raise ValueError(
            f"{path} uses project schema {project.schema_version}; "
            f"this version of Tessera reads up to {SCHEMA_VERSION}"
        )
    return project
