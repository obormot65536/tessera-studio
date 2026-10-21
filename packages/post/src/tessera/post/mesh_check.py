"""Two-mesh check (ADR-007): compare what a study reports at two mesh sizes."""

import math

from tessera.model import MeshCheck, SpotCheck, StudyResult, Verdict

MESH_CHECK_THRESHOLD = 0.05
"""Largest relative change between the two meshes that still counts as settled (ADR-007,
approved by the maintainer)."""


def spot_check(coarse: float, fine: float, threshold: float = MESH_CHECK_THRESHOLD) -> SpotCheck:
    """Compare one number at two meshes; the change is relative to the fine value."""
    change = abs(fine - coarse) / abs(fine) if fine else (math.inf if coarse else 0.0)
    verdict = Verdict.SETTLED if change <= threshold else Verdict.MESH_DEPENDENT
    return SpotCheck(coarse=coarse, fine=fine, change=change, verdict=verdict)


def compare_meshes(
    coarse: StudyResult, fine: StudyResult, threshold: float = MESH_CHECK_THRESHOLD
) -> MeshCheck:
    """The two-mesh check of ``fine`` against ``coarse``: the same study at two mesh sizes.

    Each reported number is compared with the same kind of number on the other mesh (the raw
    peak with the raw peak, face 3's peak with face 3's peak), wherever on the face it sits.
    """

    def check(a: float, b: float) -> SpotCheck:
        return spot_check(a, b, threshold)

    clear = None
    if coarse.clear_peak_von_mises is not None and fine.clear_peak_von_mises is not None:
        clear = check(coarse.clear_peak_von_mises.value, fine.clear_peak_von_mises.value)
    coarse_faces = {spot.face: spot.value for spot in coarse.face_peaks}
    faces = {
        spot.face: check(coarse_faces[spot.face], spot.value)
        for spot in fine.face_peaks
        if spot.face is not None and spot.face in coarse_faces
    }
    return MeshCheck(
        coarse_edge_length=coarse.node_flags.edge_length,
        fine_edge_length=fine.node_flags.edge_length,
        coarse_nodes=coarse.mesh.n_nodes,
        threshold=threshold,
        max_displacement=check(coarse.max_displacement, fine.max_displacement),
        peak_von_mises=check(coarse.max_von_mises, fine.max_von_mises),
        clear_peak_von_mises=clear,
        face_peaks=faces,
    )
