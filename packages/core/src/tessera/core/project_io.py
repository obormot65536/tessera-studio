"""Create, save and load Tessera projects (ADR-003).

A ``.tessera`` file is a zip container: ``project.json`` plus a copy of every part's geometry
under ``geometry/``. Any other file name is saved as plain project JSON, the hand-editable form
whose geometry paths are relative to the JSON file. An unzipped container is such a project.
"""

import atexit
import os
import re
import shutil
import tempfile
import zipfile
import zlib
from collections.abc import Callable
from pathlib import Path
from typing import Any

from pydantic import ValidationError

from tessera.model import SCHEMA_VERSION, InputError, Project

CONTAINER_SUFFIX = ".tessera"
PROJECT_ENTRY = "project.json"
MAX_ENTRY_BYTES = 1 << 30
"""Largest single entry read from a container, uncompressed (checked before decompressing)."""
MAX_TOTAL_BYTES = 4 << 30
"""Largest total geometry read from one container, uncompressed."""

_GEOMETRY_ENTRY = re.compile(r"geometry/[A-Za-z0-9_-][A-Za-z0-9_.-]*")
_SAFE_SUFFIX = re.compile(r"\.[a-z0-9]{1,10}")


def new_project(name: str) -> Project:
    return Project(name=name)


def save_project(project: Project, path: Path) -> None:
    """Save as a ``.tessera`` container, or as project JSON for any other file name.

    The file is replaced in one step, so a failed save leaves the previous file intact.
    """
    for part in project.assembly.parts:
        if not part.geometry.is_file():
            raise FileNotFoundError(f"geometry of part {part.id!r} not found: {part.geometry}")
    data = project.model_dump(mode="json")
    if path.name.endswith(CONTAINER_SUFFIX):
        _replace(path, lambda tmp: _write_container(project, data, tmp))
    else:
        folder = path.absolute().parent
        for part, stored in zip(project.assembly.parts, data["assembly"]["parts"], strict=True):
            stored["geometry"] = _relative_to(part.geometry.absolute(), folder)
        text = Project.model_validate(data).model_dump_json(indent=2)
        _replace(path, lambda tmp: tmp.write_text(text, encoding="utf-8"))


def load_project(path: Path) -> Project:
    """Load a ``.tessera`` container or project JSON, whatever the file is named.

    Geometry from a container is copied to a private temporary folder that is removed when the
    process exits; part geometry paths point there.
    """
    if zipfile.is_zipfile(path):
        try:
            return _read_container(path)
        except (zipfile.BadZipFile, zlib.error, EOFError) as exc:
            raise InputError(f"{path} is a damaged project file: {exc}") from None
    project = _parse(path.read_text(encoding="utf-8"), path)
    folder = path.absolute().parent
    for part in project.assembly.parts:
        part.geometry = folder / part.geometry  # an absolute path stays as it is
    return project


def _parse(text: str | bytes, path: Path) -> Project:
    try:
        project = Project.model_validate_json(text)
    except ValidationError as exc:
        raise InputError(f"{path} is not a valid project file:\n{exc}") from None
    if project.schema_version > SCHEMA_VERSION:
        raise InputError(
            f"{path} uses project schema {project.schema_version}; "
            f"this version of Tessera reads up to {SCHEMA_VERSION}"
        )
    return project


def _relative_to(target: Path, folder: Path) -> str:
    try:
        return Path(os.path.relpath(target, folder)).as_posix()
    except ValueError:  # Windows: a different drive has no relative path
        return target.as_posix()


def _replace(path: Path, write: Callable[[Path], object]) -> None:
    """Write through a temporary file in the target's folder, then rename it over the target."""
    fd, name = tempfile.mkstemp(dir=path.absolute().parent, prefix=f".{path.name}.", suffix=".tmp")
    os.close(fd)
    tmp = Path(name)
    try:
        write(tmp)
        os.replace(tmp, path)
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise


def _entry_name(part_id: str, suffix: str, used: set[str]) -> str:
    stem = re.sub(r"[^A-Za-z0-9_-]+", "_", part_id).strip("_") or "part"
    suffix = suffix.lower() if _SAFE_SUFFIX.fullmatch(suffix.lower()) else ""
    name, n = f"geometry/{stem}{suffix}", 1
    while name in used:
        n += 1
        name = f"geometry/{stem}-{n}{suffix}"
    used.add(name)
    return name


def _write_container(project: Project, data: dict[str, Any], tmp: Path) -> None:
    entries: dict[Path, str] = {}
    for part, stored in zip(project.assembly.parts, data["assembly"]["parts"], strict=True):
        source = part.geometry.absolute()
        if source not in entries:
            entries[source] = _entry_name(part.id, source.suffix, set(entries.values()))
        stored["geometry"] = entries[source]
    text = Project.model_validate(data).model_dump_json(indent=2)
    with zipfile.ZipFile(tmp, "w", compression=zipfile.ZIP_DEFLATED, strict_timestamps=False) as zf:
        zf.writestr(PROJECT_ENTRY, text)
        for source, entry in entries.items():
            zf.write(source, entry)


def _checked_entry(zf: zipfile.ZipFile, name: str, path: Path) -> zipfile.ZipInfo:
    try:
        info = zf.getinfo(name)
    except KeyError:
        raise InputError(f"{path} lists {name} but does not contain it") from None
    if info.file_size > MAX_ENTRY_BYTES:
        raise InputError(f"{path}: {name} is too large ({info.file_size} bytes)")
    return info


def _read_container(path: Path) -> Project:
    with zipfile.ZipFile(path) as zf:
        if PROJECT_ENTRY not in zf.namelist():
            raise InputError(f"{path} has no {PROJECT_ENTRY}; is it a Tessera project?")
        project = _parse(zf.read(_checked_entry(zf, PROJECT_ENTRY, path)), path)

        names = [part.geometry.as_posix() for part in project.assembly.parts]
        for name in names:
            if not _GEOMETRY_ENTRY.fullmatch(name):
                raise InputError(
                    f"{path}: part geometry {name!r} must be a file in the container's "
                    "geometry/ folder"
                )
        infos = {name: _checked_entry(zf, name, path) for name in set(names)}
        if sum(info.file_size for info in infos.values()) > MAX_TOTAL_BYTES:
            raise InputError(f"{path}: the geometry is too large in total")

        folder = _extraction_folder()
        extracted: dict[str, Path] = {}
        for name, info in infos.items():
            target = folder / name.removeprefix("geometry/")
            with zf.open(info) as src, target.open("wb") as dst:
                shutil.copyfileobj(src, dst)
            extracted[name] = target
    for part, name in zip(project.assembly.parts, names, strict=True):
        part.geometry = extracted[name]
    return project


def _extraction_folder() -> Path:
    folder = Path(tempfile.mkdtemp(prefix="tessera-project-"))
    atexit.register(shutil.rmtree, folder, ignore_errors=True)
    return folder
