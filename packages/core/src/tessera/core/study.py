"""Run a project's study end to end: geometry → mesh → physics → solver → post."""

import dataclasses
import math
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


MESH_CHECK_REFINEMENT = 1.0 / math.sqrt(2.0)
"""Element size of the fine mesh in a two-mesh check, as a fraction of the coarse one (ADR-007,
approved by the maintainer)."""


def solve(
    project: Project,
    workdir: Path | None = None,
    adapter: SolverAdapter | None = None,
    *,
    check_mesh: bool = False,
) -> StudyResult:
    """Mesh and solve the project's static study. v0.1: one part, no connections.

    With ``check_mesh`` the study is solved twice (ADR-007): at the project's mesh size, then
    with elements ``MESH_CHECK_REFINEMENT`` times as long. The result is the fine solve, with
    ``mesh_check`` saying which of its numbers changed by more than
    ``tessera.post.MESH_CHECK_THRESHOLD`` (5 %) between the two. The
    coarse solve's files go to ``workdir / "coarse"``.
    """
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
    material = materials.get_isotropic(part.material_id)
    if workdir is None:
        workdir = Path(tempfile.mkdtemp(prefix="tessera-"))

    def run(size: float | None, where: Path) -> StudyResult:
        job = physics.build_static_job(project.study, part, _mesh_part(part, size), material)
        where.mkdir(parents=True, exist_ok=True)
        return post.summarize(job, adapter.run(job, where), workdir=where)

    size = project.study.mesh_size
    if not check_mesh:
        return run(size, workdir)
    coarse = run(size, workdir / "coarse")
    # With an automatic size, refine what the mesher chose: its median element edge.
    fine = run((size or coarse.node_flags.edge_length) * MESH_CHECK_REFINEMENT, workdir)
    return dataclasses.replace(fine, mesh_check=post.compare_meshes(coarse, fine))
