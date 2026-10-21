"""Two-mesh check (ADR-007): each reported number compared between a coarse and a fine solve."""

import itertools
import math

import numpy as np
import pytest

from tessera.model import Mesh, NodeFlags, SolveResult, StressSpot, StudyResult, Verdict
from tessera.post import (
    MESH_CHECK_REFINEMENT_LIMIT,
    MESH_CHECK_THRESHOLD,
    compare_meshes,
    local_edge_length,
    spot_check,
)

SPOT = (1.0, 1.0, 1.0)
"""Where every hand-made peak sits: the middle of the 2 mm cube."""


def _cube(cells: int) -> Mesh:
    """The 2 mm cube split into ``cells``^3 cubes of six tets each; only the corner nodes are
    real (the midside columns repeat corners, which the edge lengths do not read)."""
    step = 2.0 / cells
    axis = np.arange(cells + 1)
    grid = np.stack(np.meshgrid(axis, axis, axis, indexing="ij"), axis=-1).reshape(-1, 3)

    def index(i: int, j: int, k: int) -> int:
        return (i * (cells + 1) + j) * (cells + 1) + k

    tets = []
    for i, j, k in itertools.product(range(cells), repeat=3):
        for perm in itertools.permutations(range(3)):
            v = [i, j, k]
            path = [index(*v)]
            for a in perm:
                v[a] += 1
                path.append(index(*v))
            tets.append(path + path[:4] + path[:2])
    return Mesh(nodes=(grid * step).astype(np.float64), elements=np.array(tets, dtype=np.int64))


def _result(
    cells: int,
    displacement: float,
    peak: float,
    clear: float | None,
    faces: dict[int, float],
) -> StudyResult:
    mesh = _cube(cells)
    n_nodes = mesh.n_nodes

    def spot(value: float, face: int | None) -> StressSpot:
        return StressSpot(value=value, node=0, location=SPOT, face=face)

    return StudyResult(
        mesh=mesh,
        solution=SolveResult(np.zeros((n_nodes, 3)), np.zeros((n_nodes, 6)), solver="fake"),
        max_displacement=displacement,
        max_von_mises=peak,
        peak_von_mises=spot(peak, None),
        clear_peak_von_mises=None if clear is None else spot(clear, None),
        face_peaks=tuple(spot(v, f) for f, v in sorted(faces.items(), key=lambda kv: -kv[1])),
        node_flags=NodeFlags(
            masks={}, n_nodes=n_nodes, distance=4.0 / cells, edge_length=2.0 / cells
        ),
    )


def test_the_approved_defaults() -> None:
    assert MESH_CHECK_THRESHOLD == 0.05
    assert MESH_CHECK_REFINEMENT_LIMIT == 0.85


def test_local_edge_length_is_the_mesh_size_at_the_point() -> None:
    # Kuhn tets of cubes with side s have edges s, s*sqrt(2) and s*sqrt(3).
    coarse = local_edge_length(_cube(2), SPOT)
    assert 1.0 <= coarse <= math.sqrt(3)
    assert local_edge_length(_cube(4), SPOT) == pytest.approx(coarse / 2)
    # The node nearest the point is used, here (0, 0, 0) of the coarse cube.
    assert local_edge_length(_cube(2), (0.1, 0.0, 0.2)) == local_edge_length(_cube(2), (0, 0, 0))


@pytest.mark.parametrize(
    ("coarse", "fine", "change", "verdict"),
    [
        (95.0, 100.0, 0.05, Verdict.SETTLED),  # exactly at the threshold still settles
        (94.9, 100.0, 0.051, Verdict.MESH_DEPENDENT),
        (105.0, 100.0, 0.05, Verdict.SETTLED),  # falling counts like rising
        (80.0, 100.0, 0.20, Verdict.MESH_DEPENDENT),
        (-1.0, -1.02, 0.02 / 1.02, Verdict.SETTLED),  # relative to the fine value's size
        (0.0, 0.0, 0.0, Verdict.SETTLED),
        (1.0, 0.0, math.inf, Verdict.MESH_DEPENDENT),
    ],
)
def test_spot_check_is_relative_to_the_fine_value(
    coarse: float, fine: float, change: float, verdict: Verdict
) -> None:
    check = spot_check(coarse, fine, refinement=0.7)
    assert (check.coarse, check.fine, check.refinement) == (coarse, fine, 0.7)
    assert check.change == pytest.approx(change)
    assert check.verdict is verdict


@pytest.mark.parametrize("refinement", [0.86, 1.0, math.nan])
def test_no_verdict_where_the_mesh_did_not_get_finer(refinement: float) -> None:
    # Unchanged, but untested: a mesh that did not refine there cannot show a singularity.
    check = spot_check(100.0, 100.0, refinement=refinement)
    assert check.verdict is Verdict.NOT_REFINED
    assert spot_check(100.0, 100.0, refinement=0.85).verdict is Verdict.SETTLED


def test_compare_meshes_matches_like_with_like() -> None:
    coarse = _result(2, 1.00, peak=120.0, clear=80.0, faces={1: 120.0, 2: 60.0, 4: 10.0})
    fine = _result(4, 1.01, peak=150.0, clear=81.0, faces={1: 150.0, 2: 61.0, 3: 5.0})
    check = compare_meshes(coarse, fine)
    assert (check.coarse_edge_length, check.fine_edge_length) == (1.0, 0.5)
    assert check.coarse_nodes == 27
    assert (check.threshold, check.refinement_limit) == (0.05, 0.85)
    assert check.max_displacement.verdict is Verdict.SETTLED
    assert check.max_displacement.refinement == pytest.approx(0.5)
    assert (check.peak_von_mises.coarse, check.peak_von_mises.fine) == (120.0, 150.0)
    assert check.peak_von_mises.refinement == pytest.approx(0.5)
    assert check.peak_von_mises.verdict is Verdict.MESH_DEPENDENT
    assert check.clear_peak_von_mises is not None
    assert check.clear_peak_von_mises.verdict is Verdict.SETTLED
    # Faces 3 and 4 have a peak on one mesh only, so they get no verdict.
    assert set(check.face_peaks) == {1, 2}
    assert check.face_peaks[1].verdict is Verdict.MESH_DEPENDENT
    assert (check.face_peaks[2].coarse, check.face_peaks[2].fine) == (60.0, 61.0)
    assert check.face_peaks[2].verdict is Verdict.SETTLED


def test_compare_meshes_judges_nothing_on_the_same_mesh() -> None:
    first = _result(2, 1.0, peak=100.0, clear=80.0, faces={1: 100.0})
    again = _result(2, 1.0, peak=100.0, clear=80.0, faces={1: 100.0})
    check = compare_meshes(first, again)
    spots = [check.max_displacement, check.peak_von_mises, check.clear_peak_von_mises]
    assert all(s is not None and s.verdict is Verdict.NOT_REFINED for s in spots)
    assert check.face_peaks[1].verdict is Verdict.NOT_REFINED


def test_no_clear_peak_on_either_mesh_means_no_clear_verdict() -> None:
    with_clear = _result(2, 1.0, peak=10.0, clear=5.0, faces={})
    without = _result(4, 1.0, peak=10.0, clear=None, faces={})
    assert compare_meshes(with_clear, without).clear_peak_von_mises is None
    assert compare_meshes(without, with_clear).clear_peak_von_mises is None


def test_a_custom_threshold_is_used_and_reported() -> None:
    coarse = _result(2, 1.0, peak=97.0, clear=None, faces={})
    fine = _result(4, 1.0, peak=100.0, clear=None, faces={})
    check = compare_meshes(coarse, fine, threshold=0.02)
    assert check.threshold == 0.02
    assert check.peak_von_mises.verdict is Verdict.MESH_DEPENDENT
