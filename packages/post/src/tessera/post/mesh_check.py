"""Two-mesh check (ADR-007): compare what a study reports at two mesh sizes."""

import math

import numpy as np
from scipy.spatial import KDTree

from tessera.model import TET10_EDGES, Mesh, MeshCheck, SpotCheck, StressSpot, StudyResult, Verdict

MESH_CHECK_THRESHOLD = 0.05
"""Largest relative change between the two meshes that still counts as settled (ADR-007,
approved by the maintainer)."""

MESH_CHECK_REFINEMENT_LIMIT = 0.85
"""A number is only judged where the fine mesh's local element edge is at most this fraction of
the coarse one's (ADR-007, approved by the maintainer); elsewhere it is ``not-refined``."""

WALL_ELEMENTS_MIN = 5.0
"""Fewest coarse-mesh elements across the wall at which a flagged spot can read settled (ADR-007
Amendment 2, approved by the maintainer); with fewer it is ``thin-wall``, because a mesh that
coarse may not show the stress rise at all."""


def local_edge_length(mesh: Mesh, point: tuple[float, float, float]) -> float:
    """Median corner-to-corner edge length of the elements at the node nearest ``point``, mm."""
    _, node = KDTree(mesh.nodes).query(point)
    return node_edge_length(mesh, int(node))


def node_edge_length(mesh: Mesh, node: int) -> float:
    """Median corner-to-corner edge length of the elements at ``node``, mm."""
    elements = mesh.elements[
        np.flatnonzero(mesh.elements.ravel() == node) // mesh.elements.shape[1]
    ]
    edges = np.array(TET10_EDGES)
    a, b = elements[:, edges[:, 0]].ravel(), elements[:, edges[:, 1]].ravel()
    return float(np.median(np.linalg.norm(mesh.nodes[a] - mesh.nodes[b], axis=1)))


def spot_check(
    coarse: float,
    fine: float,
    refinement: float,
    threshold: float = MESH_CHECK_THRESHOLD,
    refinement_limit: float = MESH_CHECK_REFINEMENT_LIMIT,
    *,
    elements_across: float | None = None,
    guard_thin_wall: bool = False,
    wall_elements_min: float = WALL_ELEMENTS_MIN,
) -> SpotCheck:
    """Compare one number at two meshes; the change is relative to the fine value.

    ``refinement`` is the local edge length ratio, fine / coarse, where the number was taken.
    ``elements_across`` is how many coarse elements span the wall there. With
    ``guard_thin_wall`` (a flagged spot), a change within the threshold on fewer than
    ``wall_elements_min`` elements across, or with no measured wall (``None``), is
    ``thin-wall``, not settled: the guard fails closed.
    """
    change = abs(fine - coarse) / abs(fine) if fine else (math.inf if coarse else 0.0)
    if not refinement <= refinement_limit:  # NaN too: no measured refinement, no verdict
        verdict = Verdict.NOT_REFINED
    elif change <= threshold:
        thin = elements_across is None or elements_across < wall_elements_min
        verdict = Verdict.THIN_WALL if guard_thin_wall and thin else Verdict.SETTLED
    else:
        verdict = Verdict.MESH_DEPENDENT
    return SpotCheck(
        coarse=coarse,
        fine=fine,
        change=change,
        refinement=refinement,
        verdict=verdict,
        elements_across=elements_across,
    )


def compare_meshes(
    coarse: StudyResult,
    fine: StudyResult,
    threshold: float = MESH_CHECK_THRESHOLD,
    refinement_limit: float = MESH_CHECK_REFINEMENT_LIMIT,
    wall_elements_min: float = WALL_ELEMENTS_MIN,
) -> MeshCheck:
    """The two-mesh check of ``fine`` against ``coarse``: the same study at two mesh sizes.

    Each reported number is compared with the same kind of number on the other mesh (the raw
    peak with the raw peak, face 3's peak with face 3's peak), wherever on the face it sits. A
    stress is judged only if the mesh got finer where the fine mesh has it; the largest
    displacement, a whole-part number, uses the median edge length of each mesh.

    Elements across the wall are the fine spot's wall thickness over the coarse mesh's local
    edge length there. Only flagged spots are guarded against a thin wall (ADR-007 Amendment 2).
    """

    def check(a: float, spot: StressSpot) -> SpotCheck:
        h_coarse = local_edge_length(coarse.mesh, spot.location)
        ratio = local_edge_length(fine.mesh, spot.location) / h_coarse
        across = None if spot.wall is None else spot.wall.thickness / h_coarse
        return spot_check(
            a,
            spot.value,
            ratio,
            threshold,
            refinement_limit,
            elements_across=across,
            guard_thin_wall=spot.flagged,
            wall_elements_min=wall_elements_min,
        )

    clear = None
    if coarse.clear_peak_von_mises is not None and fine.clear_peak_von_mises is not None:
        clear = check(coarse.clear_peak_von_mises.value, fine.clear_peak_von_mises)
    coarse_faces = {spot.face: spot.value for spot in coarse.face_peaks}
    faces = {
        spot.face: check(coarse_faces[spot.face], spot)
        for spot in fine.face_peaks
        if spot.face is not None and spot.face in coarse_faces
    }
    h_coarse, h_fine = coarse.node_flags.edge_length, fine.node_flags.edge_length
    return MeshCheck(
        coarse_edge_length=h_coarse,
        fine_edge_length=h_fine,
        coarse_nodes=coarse.mesh.n_nodes,
        threshold=threshold,
        refinement_limit=refinement_limit,
        max_displacement=spot_check(
            coarse.max_displacement,
            fine.max_displacement,
            h_fine / h_coarse,
            threshold,
            refinement_limit,
        ),
        peak_von_mises=check(coarse.max_von_mises, fine.peak_von_mises),
        clear_peak_von_mises=clear,
        face_peaks=faces,
        wall_elements_min=wall_elements_min,
    )
