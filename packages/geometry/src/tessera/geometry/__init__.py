"""STEP/STL/3MF import and geometry queries."""

from pathlib import Path

from tessera.geometry.step import describe_step_faces
from tessera.geometry.stl import describe_surface, read_stl
from tessera.model import FaceInfo, InputError

STEP_SUFFIXES = (".step", ".stp")
STL_SUFFIXES = (".stl",)


def describe_faces(path: Path) -> list[FaceInfo]:
    """Faces of a single-solid STEP or STL file, ordered by face index."""
    suffix = path.suffix.lower()
    if suffix in STL_SUFFIXES:
        return describe_surface(read_stl(path))
    if suffix in STEP_SUFFIXES:
        return describe_step_faces(path)
    raise InputError(
        f"unsupported geometry file {path.name!r}: use STEP (.step/.stp) or STL (.stl)"
    )


__all__ = [
    "STEP_SUFFIXES",
    "STL_SUFFIXES",
    "describe_faces",
    "describe_step_faces",
    "describe_surface",
    "read_stl",
]
