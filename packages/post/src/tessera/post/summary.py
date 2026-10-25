"""Summaries of a solved study: largest displacement and where the stress peaks (ADR-007)."""

from pathlib import Path

import numpy as np
from numpy.typing import NDArray

from tessera.model import (
    FlagKind,
    Mesh,
    OrthotropicMaterial,
    SolveJob,
    SolveResult,
    StrengthResult,
    StrengthSpot,
    StressFlag,
    StressSpot,
    StudyResult,
    WallSpan,
)
from tessera.post.failure import criterion, failure_mode, strength_ratio
from tessera.post.flags import SUPPORT_PREFIX, flag_nodes
from tessera.post.points import element_nodal_stresses
from tessera.post.stress import von_mises
from tessera.post.wall import OuterSurface, outer_surface, samples_on, wall_span

_FACE_PREFIX = "face:"


def _faces(mesh: Mesh) -> dict[int, NDArray[np.int64]]:
    """Node indices of each geometry face, from the mesh's ``face:<i>`` node sets."""
    faces = {}
    for name, nodes in mesh.node_sets.items():
        index = name.removeprefix(_FACE_PREFIX)
        if name.startswith(_FACE_PREFIX) and index.isdigit():
            faces[int(index)] = nodes
    return dict(sorted(faces.items()))


def _peak(values: NDArray[np.float64], candidates: NDArray[np.int64]) -> int | None:
    """The candidate node with the largest value, skipping NaN; ``None`` if there is none."""
    finite = candidates[~np.isnan(values[candidates])]
    return int(finite[np.argmax(values[finite])]) if finite.size else None


def _lowest(values: NDArray[np.float64], candidates: NDArray[np.int64]) -> int | None:
    """The candidate node with the smallest value, skipping NaN; ``None`` if there is none."""
    finite = candidates[~np.isnan(values[candidates])]
    return int(finite[np.argmin(values[finite])]) if finite.size else None


def summarize(job: SolveJob, solution: SolveResult, workdir: Path | None = None) -> StudyResult:
    """Summarize a solution of ``job``. Nodes no element uses may have no result (NaN) and are
    skipped.

    When a material of the job has a strength and the solution has point stresses, the result
    carries strength ratios (ADR-010): from each element's own stresses at its nodes
    (``element_nodal_stresses``), in its material axes, the lowest over the elements at each
    node. Without point stresses there is no strength result; averaged nodal stresses are never
    used for it."""
    mesh = job.mesh
    stress = von_mises(solution.stresses)
    flags = flag_nodes(job)
    faces = _faces(mesh)
    lowest_face = np.full(mesh.n_nodes, -1, dtype=np.int64)
    for index in reversed(faces):
        lowest_face[faces[index]] = index

    surface: OuterSurface | None = None
    walls: dict[int, WallSpan | None] = {}

    def held_faces(node_flags: tuple[StressFlag, ...]) -> NDArray[np.int64]:
        """Nodes of the face supports that flag a node: the part continues into the fixture
        there, so the wall is measured as if they were not surface (line supports hold no
        face)."""
        ids = {f.support for f in node_flags if f.kind is FlagKind.SUPPORT}
        sets = [
            mesh.node_sets[s.node_set]
            for s in job.supports
            if not s.line and s.node_set.removeprefix(SUPPORT_PREFIX) in ids
        ]
        return np.concatenate(sets) if sets else np.zeros(0, dtype=np.int64)

    def wall(node: int) -> WallSpan | None:
        """The wall at a node, measured once (ADR-007 Amendment 2)."""
        nonlocal surface
        if node not in walls:
            surface = surface or outer_surface(mesh)
            held = held_faces(flags.at(node))
            ignore = samples_on(surface, held) if held.size else None
            walls[node] = wall_span(surface, node, ignore)
        return walls[node]

    def spot(node: int, face: int | None, measure: bool = False) -> StressSpot:
        """``measure``: give the wall even if the node is not flagged."""
        x, y, z = (float(c) for c in mesh.nodes[node])
        node_flags = flags.at(node)
        return StressSpot(
            value=float(stress[node]),
            node=node,
            location=(x, y, z),
            face=face,
            flags=node_flags,
            wall=wall(node) if measure or node_flags else None,
        )

    def on_lowest_face(node: int) -> StressSpot:
        face = int(lowest_face[node])
        return spot(node, face if face >= 0 else None, measure=True)

    all_nodes = np.arange(mesh.n_nodes)
    peak = _peak(stress, all_nodes)
    if peak is None:
        raise ValueError("no node has a stress result")
    clear = _peak(stress, np.flatnonzero(~flags.flagged))
    face_peaks = [
        spot(node, face)
        for face, nodes in faces.items()
        if (node := _peak(stress, nodes)) is not None
    ]
    face_peaks.sort(key=lambda s: -s.value)

    strength = None
    materials = {m.name: m for m in job.materials}
    checked = [sec for sec in job.sections if materials[sec.material].strength is not None]
    if checked and solution.point_stresses is not None:
        nodal = element_nodal_stresses(mesh, solution.point_stresses)
        section_of = np.full(mesh.n_elements, -1, dtype=np.int64)
        element_ratio = np.full(mesh.elements.shape, np.nan)
        for index, sec in enumerate(checked):
            elements = mesh.element_sets[sec.element_set]
            section_of[elements] = index
            values = strength_ratio(
                nodal[elements].reshape(-1, 6), materials[sec.material], sec.orientation
            )
            element_ratio[elements] = values.reshape(len(elements), -1)
        ratio = np.full(mesh.n_nodes, np.nan)
        np.fmin.at(ratio, mesh.elements.ravel(), element_ratio.ravel())

        def strength_spot(node: int, face: int | None, measure: bool = False) -> StrengthSpot:
            """The spot at ``node``, from the element that gave its lowest ratio."""
            at_node = (mesh.elements == node) & (element_ratio == ratio[node])
            element, corner = np.argwhere(at_node)[0]
            sec = checked[section_of[element]]
            mode = failure_mode(nodal[element, corner], materials[sec.material], sec.orientation)
            x, y, z = (float(c) for c in mesh.nodes[node])
            node_flags = flags.at(node)
            return StrengthSpot(
                ratio=float(ratio[node]),
                node=node,
                element=int(element),
                location=(x, y, z),
                face=face,
                mode=mode,
                flags=node_flags,
                wall=wall(node) if measure or node_flags else None,
            )

        def strength_on_lowest_face(node: int) -> StrengthSpot:
            face = int(lowest_face[node])
            return strength_spot(node, face if face >= 0 else None, measure=True)

        lowest = _lowest(ratio, all_nodes)
        if lowest is not None:
            clear_lowest = _lowest(ratio, np.flatnonzero(~flags.flagged))
            face_lowest = [
                strength_spot(node, face)
                for face, nodes in faces.items()
                if (node := _lowest(ratio, nodes)) is not None
            ]
            face_lowest.sort(key=lambda s: s.ratio)
            used = {materials[sec.material] for sec in checked}
            strength = StrengthResult(
                criteria={m.name: criterion(m) for m in sorted(used, key=lambda m: m.name)},
                lowest=strength_on_lowest_face(lowest),
                clear_lowest=(
                    None if clear_lowest is None else strength_on_lowest_face(clear_lowest)
                ),
                face_lowest=tuple(face_lowest),
                assumed={
                    m.name: m.strength.assumed for m in used if m.strength and m.strength.assumed
                },
                unchecked=tuple(sorted({m.name for m in job.materials if m.strength is None})),
            )

    return StudyResult(
        mesh=mesh,
        solution=solution,
        max_displacement=float(np.nanmax(np.linalg.norm(solution.displacements, axis=1))),
        max_von_mises=float(stress[peak]),
        workdir=workdir,
        peak_von_mises=on_lowest_face(peak),
        clear_peak_von_mises=None if clear is None else on_lowest_face(clear),
        face_peaks=tuple(face_peaks),
        node_flags=flags,
        orthotropic=any(isinstance(m, OrthotropicMaterial) for m in job.materials),
        strength=strength,
    )
