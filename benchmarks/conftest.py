"""Shared helpers for the analytical validation benchmarks.

Benchmarks may use gmsh directly to build their geometry; the code under test only ever sees
the STEP file, exactly like a user import.
"""

from __future__ import annotations

import itertools
import tomllib
from collections.abc import Callable
from pathlib import Path
from typing import Any

import numpy as np
import pytest

REFERENCES_DIR = Path(__file__).parent / "references"

# Mesh-irregularity guards for the patch tests (the mesh must not be a structured grid).
MAX_AXIS_ALIGNED_EDGE_FRACTION = 0.1
MAX_SMALLEST_DIHEDRAL_DEG = 30.0

_TET4_EDGES = ((0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3))
_TET4_FACES = ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3))


def check_unstructured(nodes: np.ndarray, elements: np.ndarray) -> None:
    """The mesh must be irregular: few axis-aligned edges and some strongly skewed elements."""
    corners = nodes[elements[:, :4]]  # (M, 4, 3)
    edges = np.concatenate([corners[:, j] - corners[:, i] for i, j in _TET4_EDGES])
    edges /= np.linalg.norm(edges, axis=1)[:, None]
    aligned = float((np.abs(edges).max(axis=1) > 1.0 - 1e-9).mean())
    assert aligned < MAX_AXIS_ALIGNED_EDGE_FRACTION, (
        f"{aligned:.1%} of the element edges are axis-aligned: the mesh looks structured"
    )

    normals = []
    for a, b, c in _TET4_FACES:
        n = np.cross(corners[:, b] - corners[:, a], corners[:, c] - corners[:, a])
        normals.append(n / np.linalg.norm(n, axis=1)[:, None])
    cosines = [
        np.einsum("ij,ij->i", normals[i], normals[j])
        for i, j in itertools.combinations(range(4), 2)
    ]
    # Dihedral angle between two faces = 180 deg - angle between their (consistent) normals.
    smallest = float(180.0 - np.degrees(np.arccos(np.clip(np.min(cosines), -1.0, 1.0))))
    assert smallest < MAX_SMALLEST_DIHEDRAL_DEG, (
        f"smallest dihedral angle {smallest:.1f} deg: the mesh is not distorted"
    )


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


def write_box_stl(
    path: Path, origin: tuple[float, float, float], size: tuple[float, float, float]
) -> Path:
    """Write an axis-aligned box (mm) as a binary STL surface triangulation."""
    import gmsh

    owns_session = not gmsh.isInitialized()
    if owns_session:
        gmsh.initialize(interruptible=False)
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add(f"benchmark-{path.stem}-stl")
    try:
        gmsh.model.occ.addBox(*origin, *size)
        gmsh.model.occ.synchronize()
        gmsh.option.setNumber("Mesh.MeshSizeMax", max(size))
        gmsh.option.setNumber("Mesh.Binary", 1)
        gmsh.model.mesh.generate(2)
        gmsh.write(str(path))
    finally:
        gmsh.model.remove()
        if owns_session:
            gmsh.finalize()
    return path


@pytest.fixture(scope="session")
def box_stl(tmp_path_factory: pytest.TempPathFactory) -> Callable[..., Path]:
    """Factory: ``box_stl(name, origin, size)`` returns the path of a box STL file."""
    directory = tmp_path_factory.mktemp("geometry-stl")

    def make(
        name: str, origin: tuple[float, float, float], size: tuple[float, float, float]
    ) -> Path:
        path = directory / f"{name}.stl"
        if not path.exists():
            write_box_stl(path, origin, size)
        return path

    return make


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


@pytest.fixture(scope="session")
def assert_unstructured() -> Callable[[np.ndarray, np.ndarray], None]:
    """``assert_unstructured(nodes, elements)``: fail unless the tet mesh is irregular.

    Shared by the patch tests (TS-03, TO-01), so a future structured mesher cannot turn them
    into trivial tests.
    """
    return check_unstructured
