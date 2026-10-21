"""Run a project's study end to end: geometry → mesh → physics → solver → post."""

import tempfile
from pathlib import Path

import tessera.geometry as geometry
import tessera.materials as materials
import tessera.mesh as mesh
import tessera.physics as physics
import tessera.post as post
from tessera.core._checks import check_geometry, part_named, validated
from tessera.model import (
    AnalysisKind,
    EdgeInfo,
    FaceInfo,
    InputError,
    Mesh,
    PartInstance,
    Project,
    SolverAdapter,
    StudyResult,
    Transform,
)


def describe_faces(project: Project, part_id: str) -> list[FaceInfo]:
    """Faces of a part, for choosing where supports and loads go."""
    return geometry.describe_faces(part_named(validated(project), part_id).geometry)


def describe_edges(project: Project, part_id: str) -> list[EdgeInfo]:
    """Edges of a part, each named by its two faces, for choosing line supports (ADR-006)."""
    return geometry.describe_edges(part_named(validated(project), part_id).geometry)


def _mesh_part(part: PartInstance, size: float | None) -> Mesh:
    suffix = part.geometry.suffix.lower()
    if suffix in geometry.STL_SUFFIXES:
        surface = geometry.read_stl(part.geometry)
        return mesh.mesh_surface(surface, part_id=part.id, size=size)
    if suffix in geometry.STEP_SUFFIXES:
        return mesh.mesh_step(part.geometry, part_id=part.id, size=size)
    raise InputError(f"unsupported geometry file {part.geometry.name!r}: use STEP or STL")


def _default_adapter() -> SolverAdapter:
    # v0.1 has one solver; ADR-004's registry by analysis kind arrives with the second one.
    from tessera.solvers import CalculixAdapter

    return CalculixAdapter()


def solve(
    project: Project, workdir: Path | None = None, adapter: SolverAdapter | None = None
) -> StudyResult:
    """Mesh and solve the project's static study. v0.1: one part, no connections."""
    project = validated(project)
    parts = project.assembly.parts
    if len(parts) != 1 or project.assembly.connections:
        raise NotImplementedError("v0.1 solves a single part; assemblies arrive in v0.5")
    part = check_geometry(parts[0])
    if part.transform != Transform():
        raise NotImplementedError("v0.1 does not apply part placements yet")
    if part.material_id is None:
        raise InputError(f"part {part.id!r} has no material")

    adapter = adapter or _default_adapter()
    analysis = AnalysisKind(project.study.analysis)
    if analysis not in adapter.capabilities():
        raise InputError(f"solver {adapter.name!r} cannot run a {analysis} analysis")

    part_mesh = _mesh_part(part, project.study.mesh_size)
    job = physics.build_static_job(
        project.study, part, part_mesh, materials.get_isotropic(part.material_id)
    )
    if workdir is None:
        workdir = Path(tempfile.mkdtemp(prefix="tessera-"))
    workdir.mkdir(parents=True, exist_ok=True)
    solution = adapter.run(job, workdir)
    return post.summarize(job, solution, workdir=workdir)
