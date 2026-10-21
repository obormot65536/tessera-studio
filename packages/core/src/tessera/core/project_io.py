"""Create, save and load Tessera projects (ADR-003).

A ``.tessera`` file is a zip container: ``project.json`` plus a copy of every part's geometry
under ``geometry/``. Any other file name is saved as plain project JSON, the hand-editable form
whose geometry paths are relative to the JSON file. An unzipped container is such a project.
"""

import atexit
import json
import os
import re
import shutil
import stat
import tempfile
import zipfile
import zlib
from collections.abc import Callable
from pathlib import Path, PurePosixPath
from typing import Any

from pydantic import ValidationError

import tessera.geometry as geometry
from tessera.model import SCHEMA_VERSION, InputError, Project

CONTAINER_SUFFIX = ".tessera"
PROJECT_ENTRY = "project.json"
MAX_ENTRY_BYTES = 1 << 30
"""Largest single entry read from a container, uncompressed (checked before decompressing)."""
MAX_TOTAL_BYTES = 4 << 30
"""Largest total geometry read from one container, uncompressed."""
GEOMETRY_SUFFIXES = (*geometry.STEP_SUFFIXES, *geometry.STL_SUFFIXES)
"""Only these files are copied into a container, so a crafted project cannot pack other files."""

_GEOMETRY_ENTRY = re.compile(r"geometry/[A-Za-z0-9_-][A-Za-z0-9_.-]*")
_extraction_folders: list[Path] = []


def new_project(name: str) -> Project:
    return Project(name=name)


def is_container_name(path: Path) -> bool:
    return path.name.lower().endswith(CONTAINER_SUFFIX)


def save_project(project: Project, path: Path) -> None:
    """Save as a ``.tessera`` container, or as project JSON for any other file name.

    The file is replaced in one step, so a failed save leaves the previous file intact.
    """
    for part in project.assembly.parts:
        if not part.geometry.is_file():
            raise FileNotFoundError(f"geometry of part {part.id!r} not found: {part.geometry}")
    data = project.model_dump(mode="json")
    # The file is written in this version's schema, whatever version it was read from.
    data["schema_version"] = SCHEMA_VERSION
    stored_parts = data["assembly"]["parts"]
    if is_container_name(path):
        _replace(path, lambda tmp: _write_container(project, data, tmp))
        return
    folder = path.absolute().parent
    for part, stored in zip(project.assembly.parts, stored_parts, strict=True):
        if any(part.geometry.is_relative_to(f) for f in _extraction_folders):
            raise InputError(
                f"part {part.id!r} uses geometry from an opened .tessera file, which is only "
                "kept until Tessera exits; save as .tessera, or unzip the file to edit it as JSON"
            )
        stored["geometry"] = _relative_to(part.geometry.absolute(), folder)
    _replace(path, lambda tmp: tmp.write_text(_dump(data), encoding="utf-8"))


def load_project(path: Path) -> Project:
    """Load a ``.tessera`` container or project JSON, whatever the file is named.

    Geometry from a container is copied to a private temporary folder that is removed when the
    process exits; part geometry paths point there.
    """
    if path.exists() and not path.is_file():
        raise InputError(f"{path} is not a file")
    if zipfile.is_zipfile(path):
        try:
            return _read_container(path)
        except (zipfile.BadZipFile, zlib.error, EOFError, RuntimeError, NotImplementedError) as exc:
            # RuntimeError: encrypted entry; NotImplementedError: unsupported compression.
            raise InputError(f"{path} is a damaged or unsupported project file: {exc}") from None
    project = _parse(path.read_text(encoding="utf-8"), path)
    folder = path.absolute().parent
    for part in project.assembly.parts:
        part.geometry = folder / part.geometry  # an absolute path stays as it is
    return project


def _parse(text: str | bytes, path: Path) -> Project:
    try:
        raw = json.loads(text)
    except ValueError as exc:  # JSONDecodeError, or bytes that are not UTF-8
        raise InputError(f"{path} is not a valid project file:\n{exc}") from None
    # Check the version first: a newer file may use fields this version does not know.
    version = raw.get("schema_version") if isinstance(raw, dict) else None
    if isinstance(version, int) and version > SCHEMA_VERSION:
        raise InputError(
            f"{path} uses project schema {version}; "
            f"this version of Tessera reads up to {SCHEMA_VERSION}"
        )
    try:
        return Project.model_validate(raw)
    except ValidationError as exc:
        raise InputError(f"{path} is not a valid project file:\n{exc}") from None


def _dump(data: dict[str, Any]) -> str:
    """Validated JSON text of ``data``, with geometry paths written exactly as given.

    Dumping through the model would turn them back into ``Path`` and use ``\\`` on Windows.
    """
    try:
        Project.model_validate(data)
    except ValidationError as exc:
        raise InputError(f"project {data.get('name')!r} is not valid:\n{exc}") from None
    return json.dumps(data, indent=2, ensure_ascii=False)


def _relative_to(target: Path, folder: Path) -> str:
    try:
        return Path(os.path.relpath(target, folder)).as_posix()
    except ValueError:  # Windows: a different drive has no relative path
        return target.as_posix()


def _replace(path: Path, write: Callable[[Path], object]) -> None:
    """Write through a temporary file in the target's folder, then rename it over the target.

    A symlinked target is updated where it points; the file keeps its permissions (new files
    get the usual ones for the user's umask).
    """
    target = path.resolve() if path.is_symlink() else path
    fd, name = tempfile.mkstemp(
        dir=target.absolute().parent, prefix=f".{target.name}.", suffix=".tmp"
    )
    os.close(fd)
    tmp = Path(name)
    try:
        write(tmp)
        if target.exists():
            mode = stat.S_IMODE(target.stat().st_mode)
        else:
            umask = os.umask(0)
            os.umask(umask)
            mode = 0o666 & ~umask
        tmp.chmod(mode)
        with tmp.open("rb+") as f:
            os.fsync(f.fileno())
        os.replace(tmp, target)
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise


def _entry_name(part_id: str, suffix: str, used: set[str]) -> str:
    """A container entry for a part's geometry, unique even on case-insensitive file systems."""
    stem = re.sub(r"[^A-Za-z0-9_-]+", "_", part_id).strip("_") or "part"
    name, n = f"geometry/{stem}{suffix}", 1
    while name.casefold() in used:
        n += 1
        name = f"geometry/{stem}-{n}{suffix}"
    used.add(name.casefold())
    return name


def _write_container(project: Project, data: dict[str, Any], tmp: Path) -> None:
    entries: dict[Path, str] = {}
    used: set[str] = set()
    for part, stored in zip(project.assembly.parts, data["assembly"]["parts"], strict=True):
        source = part.geometry.resolve()
        suffix = source.suffix.lower()
        if suffix not in GEOMETRY_SUFFIXES:
            raise InputError(
                f"part {part.id!r}: {part.geometry.name!r} is not a STEP or STL file, "
                "so it is not packed into the project"
            )
        if source not in entries:
            entries[source] = _entry_name(part.id, suffix, used)
        stored["geometry"] = entries[source]
    text = _dump(data)
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
        infos = {name: _checked_entry(zf, name, path) for name in dict.fromkeys(names)}
        if sum(info.file_size for info in infos.values()) > MAX_TOTAL_BYTES:
            raise InputError(f"{path}: the geometry is too large in total")

        # Files get names Tessera makes up, so entry names cannot collide or address devices.
        folder = _extraction_folder()
        extracted: dict[str, Path] = {}
        for i, (name, info) in enumerate(infos.items()):
            target = folder / f"{i}{PurePosixPath(name).suffix.lower()}"
            with zf.open(info) as src, target.open("wb") as dst:
                shutil.copyfileobj(src, dst)
            extracted[name] = target
    for part, name in zip(project.assembly.parts, names, strict=True):
        part.geometry = extracted[name]
    return project


def _extraction_folder() -> Path:
    folder = Path(tempfile.mkdtemp(prefix="tessera-project-"))
    _extraction_folders.append(folder)
    atexit.register(shutil.rmtree, folder, ignore_errors=True)
    return folder
