"""Scoped access to gmsh's global state."""

from collections.abc import Iterator
from contextlib import contextmanager
from itertools import count

import gmsh

_models = count()


@contextmanager
def gmsh_model() -> Iterator[None]:
    """Run a block in a fresh gmsh model; initialize/finalize gmsh only if we started it."""
    owner = not gmsh.isInitialized()
    if owner:
        gmsh.initialize(interruptible=False)
    gmsh.option.setNumber("General.Terminal", 0)
    previous = gmsh.model.getCurrent() if not owner else None
    name = f"tessera-mesh-{next(_models)}"
    gmsh.model.add(name)
    try:
        yield
    finally:
        gmsh.model.setCurrent(name)
        gmsh.model.remove()
        if previous:
            gmsh.model.setCurrent(previous)
        if owner:
            gmsh.finalize()
