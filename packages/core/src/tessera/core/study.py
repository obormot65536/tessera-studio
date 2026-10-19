"""Run a project's study end to end: geometry → mesh → physics → solver → post."""

import tempfile
from pathlib import Path

import tessera.geometry as geometry
import tessera.materials as materials
import tessera.mesh as mesh
import tessera.physics as physics
import tessera.post as post
from tessera.model import (
    AnalysisKind,
    FaceInfo,
    InputError,
    PartInstance,
    Project,
    SolverAdapter,
    StudyResult,
    Transform,
)


def _part(project: Project, part_id: str) -> PartInstance:
    for part in project.assembly.parts:
        if part.id == part_id:
            return part
    raise InputError(f"no part {part_id!r} in project {project.name!r}")


def describe_faces(project: Project, part_id: str) -> list[FaceInfo]:
    """Faces of a part, for choosing where supports and loads go."""
    return geometry.describe_faces(_part(project, part_id).geometry)


def _default_adapter() -> SolverAdapter:
    # v0.1 has one solver; ADR-004's registry by analysis kind arrives with the second one.
    from tessera.solvers import CalculixAdapter

    return CalculixAdapter()


def solve(
    project: Project, workdir: Path | None = None, adapter: SolverAdapter | None = None
) -> StudyResult:
    """Mesh and solve the project's static study. v0.1: one part, no connections."""
    parts = project.assembly.parts
    if len(parts) != 1 or project.assembly.connections:
        raise NotImplementedError("v0.1 solves a single part; assemblies arrive in v0.5")
    part = parts[0]
    if part.transform != Transform():
        raise NotImplementedError("v0.1 does not apply part placements yet")
    if part.material_id is None:
        raise InputError(f"part {part.id!r} has no material")

    adapter = adapter or _default_adapter()
    analysis = AnalysisKind(project.study.analysis)
    if analysis not in adapter.capabilities():
        raise InputError(f"solver {adapter.name!r} cannot run a {analysis} analysis")

    part_mesh = mesh.mesh_step(part.geometry, part_id=part.id, size=project.study.mesh_size)
    job = physics.build_static_job(
        project.study, part, part_mesh, materials.get_isotropic(part.material_id)
    )
    if workdir is None:
        workdir = Path(tempfile.mkdtemp(prefix="tessera-"))
    workdir.mkdir(parents=True, exist_ok=True)
    solution = adapter.run(job, workdir)
    return post.summarize(job.mesh, solution, workdir=workdir)
