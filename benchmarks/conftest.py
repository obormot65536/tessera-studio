"""Shared helpers for the analytical validation benchmarks.

Benchmarks may use gmsh directly to build their geometry; the code under test only ever sees
the STEP file, exactly like a user import.
"""

from __future__ import annotations

import tomllib
from collections.abc import Callable
from pathlib import Path
from typing import Any

import pytest

REFERENCES_DIR = Path(__file__).parent / "references"


def write_box_step(
    path: Path, origin: tuple[float, float, float], size: tuple[float, float, float]
) -> Path:
    """Write an axis-aligned box (mm) as a STEP file with the gmsh OpenCASCADE kernel."""
    import gmsh

    owns_session = not gmsh.isInitialized()
    if owns_session:
        gmsh.initialize(interruptible=False)
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add(f"benchmark-{path.stem}")
    try:
        gmsh.model.occ.addBox(*origin, *size)
        gmsh.model.occ.synchronize()
        gmsh.write(str(path))
    finally:
        gmsh.model.remove()
        if owns_session:
            gmsh.finalize()
    return path


@pytest.fixture(scope="session")
def box_step(tmp_path_factory: pytest.TempPathFactory) -> Callable[..., Path]:
    """Factory: ``box_step(name, origin, size)`` returns the path of a box STEP file."""
    directory = tmp_path_factory.mktemp("geometry")

    def make(
        name: str, origin: tuple[float, float, float], size: tuple[float, float, float]
    ) -> Path:
        path = directory / f"{name}.step"
        if not path.exists():
            write_box_step(path, origin, size)
        return path

    return make


@pytest.fixture(scope="session")
def load_reference() -> Callable[[str], dict[str, Any]]:
    """Factory: ``load_reference(name)`` reads ``references/<name>.toml``.

    Skips the test while the maintainer has not yet approved (committed) the reference.
    """

    def load(name: str) -> dict[str, Any]:
        path = REFERENCES_DIR / f"{name}.toml"
        if not path.is_file():
            pytest.skip("reference awaiting maintainer approval")
        with path.open("rb") as f:
            return tomllib.load(f)

    return load
