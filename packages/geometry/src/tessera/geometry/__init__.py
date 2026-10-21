"""STEP/STL/3MF import and geometry queries."""

import hashlib
import threading
from collections import OrderedDict
from pathlib import Path

from tessera.geometry.edges import describe_surface_edges
from tessera.geometry.step import describe_step
from tessera.geometry.stl import describe_surface, read_stl
from tessera.model import EdgeInfo, FaceInfo, InputError

STEP_SUFFIXES = (".step", ".stp")
STL_SUFFIXES = (".stl",)
FACE_CACHE_SIZE = 32
"""How many files' faces and edges are kept in memory."""

_Description = tuple[tuple[FaceInfo, ...], tuple[EdgeInfo, ...]]
_face_cache: OrderedDict[tuple[str, bytes], _Description] = OrderedDict()
_face_cache_lock = threading.Lock()


def describe_faces(path: Path) -> list[FaceInfo]:
    """Faces of a single-solid STEP or STL file, ordered by face index.

    Faces are remembered by file content, so asking again for an unchanged file (for every
    support and load added to a part) does not parse it again, and a changed file is re-read.
    """
    return list(_describe(path)[0])


def describe_edges(path: Path) -> list[EdgeInfo]:
    """Edges of a single-solid STEP or STL file, each named by its two faces (ADR-006).

    Ordered by face pair. Seams and zero-length curves are left out. Remembered together with
    the faces, by file content.
    """
    return list(_describe(path)[1])


def _describe(path: Path) -> _Description:
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
            return _face_cache[key]
    if kind == "stl":
        surface = read_stl(path)
        found = (tuple(describe_surface(surface)), tuple(describe_surface_edges(surface)))
    else:
        faces, edges = describe_step(path)
        found = (tuple(faces), tuple(edges))
    # Not cached if the file changed while it was being parsed. A change that is undone before
    # the second hash goes unnoticed; that needs two writes within one parse.
    try:
        unchanged = _digest(path) == digest
    except OSError:
        unchanged = False
    if unchanged:
        with _face_cache_lock:
            _face_cache[key] = found
            while len(_face_cache) > FACE_CACHE_SIZE:
                _face_cache.popitem(last=False)
    return found


def clear_face_cache() -> None:
    """Forget all remembered faces and edges, so the next query reads its file again."""
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
    "describe_edges",
    "describe_faces",
    "describe_step",
    "describe_surface",
    "describe_surface_edges",
    "read_stl",
]
