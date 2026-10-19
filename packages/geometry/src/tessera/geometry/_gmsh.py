"""Scoped access to gmsh's global state."""

from collections.abc import Iterator, Mapping
from contextlib import contextmanager
from itertools import count

import gmsh

_models = count()

# STEP lengths are converted to mm on import (ADR-002), whatever the file's own unit.
_BASE_OPTIONS: dict[str, float | str] = {"General.Terminal": 0, "Geometry.OCCTargetUnit": "MM"}


def _get(name: str, like: float | str) -> float | str:
    return gmsh.option.getString(name) if isinstance(like, str) else gmsh.option.getNumber(name)


def _set(name: str, value: float | str) -> None:
    if isinstance(value, str):
        gmsh.option.setString(name, value)
    else:
        gmsh.option.setNumber(name, value)


@contextmanager
def gmsh_model(options: Mapping[str, float | str] | None = None) -> Iterator[None]:
    """Run a block in a fresh gmsh model with ``options`` set.

    gmsh options are process-global, so every option set here is restored afterwards and a
    long-lived process (GUI, MCP server) is never left with this call's settings.
    """
    owner = not gmsh.isInitialized()
    if owner:
        gmsh.initialize(interruptible=False)
    wanted = {**_BASE_OPTIONS, **(options or {})}
    saved = {name: _get(name, value) for name, value in wanted.items()}
    previous = None if owner else gmsh.model.getCurrent()
    name = f"tessera-geometry-{next(_models)}"
    try:
        for option, value in wanted.items():
            _set(option, value)
        gmsh.model.add(name)
        yield
    finally:
        if name in gmsh.model.list():
            gmsh.model.setCurrent(name)
            gmsh.model.remove()
        if previous:
            gmsh.model.setCurrent(previous)
        for option, value in saved.items():
            _set(option, value)
        if owner:
            gmsh.finalize()
