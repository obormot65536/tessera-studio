"""Two-mesh check (ADR-007): each reported number compared between a coarse and a fine solve."""

import math

import numpy as np
import pytest

from tessera.model import Mesh, NodeFlags, SolveResult, StressSpot, StudyResult, Verdict
from tessera.post import MESH_CHECK_THRESHOLD, compare_meshes, spot_check


def _result(
    n_nodes: int,
    edge: float,
    displacement: float,
    peak: float,
    clear: float | None,
    faces: dict[int, float],
) -> StudyResult:
    def spot(value: float, face: int | None) -> StressSpot:
        return StressSpot(value=value, node=0, location=(0.0, 0.0, 0.0), face=face)

    return StudyResult(
        mesh=Mesh(nodes=np.zeros((n_nodes, 3)), elements=np.zeros((1, 10), dtype=np.int64)),
        solution=SolveResult(np.zeros((n_nodes, 3)), np.zeros((n_nodes, 6)), solver="fake"),
        max_displacement=displacement,
        max_von_mises=peak,
        peak_von_mises=spot(peak, None),
        clear_peak_von_mises=None if clear is None else spot(clear, None),
        face_peaks=tuple(spot(v, f) for f, v in sorted(faces.items(), key=lambda kv: -kv[1])),
        node_flags=NodeFlags(masks={}, n_nodes=n_nodes, distance=2 * edge, edge_length=edge),
    )


def test_the_approved_threshold_is_five_percent() -> None:
    assert MESH_CHECK_THRESHOLD == 0.05


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
    check = spot_check(coarse, fine)
    assert (check.coarse, check.fine) == (coarse, fine)
    assert check.change == pytest.approx(change)
    assert check.verdict is verdict


def test_compare_meshes_matches_like_with_like() -> None:
    coarse = _result(100, 3.0, 1.00, peak=120.0, clear=80.0, faces={1: 120.0, 2: 60.0, 4: 10.0})
    fine = _result(280, 2.1, 1.01, peak=150.0, clear=81.0, faces={1: 150.0, 2: 61.0, 3: 5.0})
    check = compare_meshes(coarse, fine)
    assert (check.coarse_edge_length, check.fine_edge_length) == (3.0, 2.1)
    assert check.coarse_nodes == 100
    assert check.threshold == MESH_CHECK_THRESHOLD
    assert check.max_displacement.verdict is Verdict.SETTLED
    assert (check.peak_von_mises.coarse, check.peak_von_mises.fine) == (120.0, 150.0)
    assert check.peak_von_mises.verdict is Verdict.MESH_DEPENDENT
    assert check.clear_peak_von_mises is not None
    assert check.clear_peak_von_mises.verdict is Verdict.SETTLED
    # Faces 3 and 4 have a peak on one mesh only, so they get no verdict.
    assert set(check.face_peaks) == {1, 2}
    assert check.face_peaks[1].verdict is Verdict.MESH_DEPENDENT
    assert (check.face_peaks[2].coarse, check.face_peaks[2].fine) == (60.0, 61.0)
    assert check.face_peaks[2].verdict is Verdict.SETTLED


def test_no_clear_peak_on_either_mesh_means_no_clear_verdict() -> None:
    with_clear = _result(10, 1.0, 1.0, peak=10.0, clear=5.0, faces={})
    without = _result(10, 1.0, 1.0, peak=10.0, clear=None, faces={})
    assert compare_meshes(with_clear, without).clear_peak_von_mises is None
    assert compare_meshes(without, with_clear).clear_peak_von_mises is None


def test_a_custom_threshold_is_used_and_reported() -> None:
    coarse = _result(10, 1.0, 1.0, peak=97.0, clear=None, faces={})
    fine = _result(20, 0.7, 1.0, peak=100.0, clear=None, faces={})
    check = compare_meshes(coarse, fine, threshold=0.02)
    assert check.threshold == 0.02
    assert check.peak_von_mises.verdict is Verdict.MESH_DEPENDENT
