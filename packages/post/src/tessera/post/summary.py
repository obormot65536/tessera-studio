"""Summaries of a solved study: largest displacement and where the stress peaks (ADR-007)."""

from pathlib import Path

import numpy as np
from numpy.typing import NDArray

from tessera.model import (
    Mesh,
    OrthotropicMaterial,
    SolveJob,
    SolveResult,
    StressSpot,
    StudyResult,
)
from tessera.post.flags import flag_nodes

_FACE_PREFIX = "face:"


def von_mises(stresses: NDArray[np.float64]) -> NDArray[np.float64]:
    """Von Mises stress from (..., 6) components xx, yy, zz, xy, yz, zx (MPa), shape (...).

    (N, 6) nodal stresses give (N,); (M, Q, 6) integration-point stresses give (M, Q).
    """
    if stresses.shape[-1:] != (6,):
        raise ValueError(f"stresses must have shape (..., 6), got {stresses.shape}")
    sxx, syy, szz, sxy, syz, szx = (stresses[..., k] for k in range(6))
    return np.sqrt(
        0.5 * ((sxx - syy) ** 2 + (syy - szz) ** 2 + (szz - sxx) ** 2)
        + 3.0 * (sxy**2 + syz**2 + szx**2)
    )


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


def summarize(job: SolveJob, solution: SolveResult, workdir: Path | None = None) -> StudyResult:
    """Summarize a solution of ``job``. Nodes no element uses may have no result (NaN) and are
    skipped."""
    mesh = job.mesh
    stress = von_mises(solution.stresses)
    flags = flag_nodes(job)
    faces = _faces(mesh)
    lowest_face = np.full(mesh.n_nodes, -1, dtype=np.int64)
    for index in reversed(faces):
        lowest_face[faces[index]] = index

    def spot(node: int, face: int | None) -> StressSpot:
        x, y, z = (float(c) for c in mesh.nodes[node])
        return StressSpot(
            value=float(stress[node]),
            node=node,
            location=(x, y, z),
            face=face,
            flags=flags.at(node),
        )

    def on_lowest_face(node: int) -> StressSpot:
        face = int(lowest_face[node])
        return spot(node, face if face >= 0 else None)

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
    )
