"""STEP/STL/3MF import and geometry queries."""

import hashlib
import threading
from collections import OrderedDict
from pathlib import Path

from tessera.geometry.step import describe_step_faces
from tessera.geometry.stl import describe_surface, read_stl
from tessera.model import FaceInfo, InputError

STEP_SUFFIXES = (".step", ".stp")
STL_SUFFIXES = (".stl",)
FACE_CACHE_SIZE = 32
"""How many files' faces are kept in memory."""

_face_cache: OrderedDict[tuple[str, bytes], tuple[FaceInfo, ...]] = OrderedDict()
_face_cache_lock = threading.Lock()


def describe_faces(path: Path) -> list[FaceInfo]:
    """Faces of a single-solid STEP or STL file, ordered by face index.

    Faces are remembered by file content, so asking again for an unchanged file (for every
    support and load added to a part) does not parse it again, and a changed file is re-read.
    """
    suffix = path.suffix.lower()
    if suffix in STL_SUFFIXES:
        kind = "stl"
    elif suffix in STEP_SUFFIXES:
        kind = "step"
    else:
        raise InputError(
            f"unsupported geometry file {path.name!r}: use STEP (.step/.stp) or STL (.stl)"
        )
    if not path.is_file():
        raise FileNotFoundError(path)
    try:
        digest = _digest(path)
    except OSError as exc:
        raise InputError(f"cannot read {path}: {exc}") from None
    # STL faces also depend on read_stl's feature angle; once that becomes a setting, it must
    # be part of the key, or faces and mesh would be split at different edges (ADR-005).
    key = (kind, digest)
    with _face_cache_lock:
        if key in _face_cache:
            _face_cache.move_to_end(key)
            return list(_face_cache[key])
    faces = describe_surface(read_stl(path)) if kind == "stl" else describe_step_faces(path)
    # Not cached if the file changed while it was being parsed. A change that is undone before
    # the second hash goes unnoticed; that needs two writes within one parse.
    try:
        unchanged = _digest(path) == digest
    except OSError:
        unchanged = False
    if unchanged:
        with _face_cache_lock:
            _face_cache[key] = tuple(faces)
            while len(_face_cache) > FACE_CACHE_SIZE:
                _face_cache.popitem(last=False)
    return list(faces)


def clear_face_cache() -> None:
    """Forget all remembered faces, so the next ``describe_faces`` reads its file again."""
    with _face_cache_lock:
        _face_cache.clear()


def _digest(path: Path) -> bytes:
    with path.open("rb") as f:
        return hashlib.file_digest(f, "sha256").digest()


__all__ = [
    "FACE_CACHE_SIZE",
    "STEP_SUFFIXES",
    "STL_SUFFIXES",
    "clear_face_cache",
    "describe_faces",
    "describe_step_faces",
    "describe_surface",
    "read_stl",
]
