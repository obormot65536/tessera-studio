"""Two-mesh check (ADR-007): each reported number compared between a coarse and a fine solve."""

import itertools
import math

import numpy as np
import pytest

from tessera.model import (
    FlagKind,
    Mesh,
    NodeFlags,
    SolveResult,
    StressFlag,
    StressSpot,
    StudyResult,
    Verdict,
    WallSpan,
)
from tessera.post import (
    MESH_CHECK_REFINEMENT_LIMIT,
    MESH_CHECK_THRESHOLD,
    WALL_ELEMENTS_MIN,
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
    flagged: bool = False,
    wall: float | None = None,
) -> StudyResult:
    """Every spot at ``SPOT``; the raw peak and the face peaks flagged if asked, and every spot
    on a wall ``wall`` mm thick if given (the clear peak is never flagged)."""
    mesh = _cube(cells)
    n_nodes = mesh.n_nodes
    span = None if wall is None else WallSpan(wall, wall * cells / 2, (0.0, 0.0, -1.0))
    flags = (StressFlag(FlagKind.SUPPORT, support="root"),) if flagged else ()

    def spot(value: float, face: int | None, flagged: bool = True) -> StressSpot:
        return StressSpot(
            value=value, node=0, location=SPOT, face=face, flags=flags if flagged else (), wall=span
        )

    return StudyResult(
        mesh=mesh,
        solution=SolveResult(np.zeros((n_nodes, 3)), np.zeros((n_nodes, 6)), solver="fake"),
        max_displacement=displacement,
        max_von_mises=peak,
        peak_von_mises=spot(peak, None),
        clear_peak_von_mises=None if clear is None else spot(clear, None, flagged=False),
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


# --- Thin walls (ADR-007 Amendment 2) ----------------------------------------------------


def test_the_approved_wall_minimum() -> None:
    assert WALL_ELEMENTS_MIN == 5.0


@pytest.mark.parametrize(
    ("coarse", "refinement", "across", "guard", "verdict"),
    [
        (99.0, 0.7, 2.0, True, Verdict.THIN_WALL),  # settled-looking, but too few elements
        (99.0, 0.7, 4.99, True, Verdict.THIN_WALL),
        (99.0, 0.7, 5.0, True, Verdict.SETTLED),  # exactly the minimum is enough
        (99.0, 0.7, 2.0, False, Verdict.SETTLED),  # not guarded: unflagged spots
        (99.0, 0.7, None, True, Verdict.SETTLED),  # no wall measured
        (80.0, 0.7, 2.0, True, Verdict.MESH_DEPENDENT),  # a measured failure needs no guard
        (99.0, 0.9, 2.0, True, Verdict.NOT_REFINED),  # no refinement comes first
    ],
)
def test_thin_wall_guard_in_spot_check(
    coarse: float, refinement: float, across: float | None, guard: bool, verdict: Verdict
) -> None:
    check = spot_check(coarse, 100.0, refinement, elements_across=across, guard_thin_wall=guard)
    assert check.verdict is verdict
    assert check.elements_across == across


def test_spot_check_without_a_wall_is_unchanged() -> None:
    check = spot_check(99.0, 100.0, 0.7)
    assert (check.verdict, check.elements_across) == (Verdict.SETTLED, None)


def test_compare_meshes_guards_flagged_spots_on_thin_walls() -> None:
    # A 3 mm wall: about two coarse elements across (Kuhn cubes of 1 mm, median edge ~1.4 mm).
    coarse = _result(2, 1.0, peak=99.0, clear=80.0, faces={1: 99.0, 2: 60.0}, flagged=True)
    fine = _result(
        4, 1.0, peak=100.0, clear=80.5, faces={1: 100.0, 2: 61.0}, flagged=True, wall=3.0
    )
    check = compare_meshes(coarse, fine)
    across = 3.0 / local_edge_length(coarse.mesh, SPOT)
    assert across < WALL_ELEMENTS_MIN
    assert check.wall_elements_min == WALL_ELEMENTS_MIN
    assert check.peak_von_mises.verdict is Verdict.THIN_WALL
    assert check.peak_von_mises.elements_across == pytest.approx(across)
    assert check.face_peaks[1].verdict is Verdict.THIN_WALL
    assert check.face_peaks[2].verdict is Verdict.THIN_WALL
    # The clear peak is not flagged, so it is never guarded.
    assert check.clear_peak_von_mises is not None
    assert check.clear_peak_von_mises.verdict is Verdict.SETTLED
    assert check.clear_peak_von_mises.elements_across == pytest.approx(across)
    assert check.max_displacement.elements_across is None


def test_compare_meshes_judges_flagged_spots_on_walls_thick_enough() -> None:
    coarse = _result(2, 1.0, peak=99.0, clear=None, faces={}, flagged=True)
    fine = _result(4, 1.0, peak=100.0, clear=None, faces={}, flagged=True, wall=8.0)
    assert compare_meshes(coarse, fine).peak_von_mises.verdict is Verdict.SETTLED
    # A custom minimum is used and reported.
    check = compare_meshes(coarse, fine, wall_elements_min=7.0)
    assert check.wall_elements_min == 7.0
    assert check.peak_von_mises.verdict is Verdict.THIN_WALL


def test_compare_meshes_does_not_guard_unflagged_or_unmeasured_spots() -> None:
    coarse = _result(2, 1.0, peak=99.0, clear=None, faces={})
    unflagged = _result(4, 1.0, peak=100.0, clear=None, faces={}, wall=1.0)
    assert compare_meshes(coarse, unflagged).peak_von_mises.verdict is Verdict.SETTLED
    no_wall = _result(4, 1.0, peak=100.0, clear=None, faces={}, flagged=True)
    check = compare_meshes(coarse, no_wall).peak_von_mises
    assert (check.verdict, check.elements_across) == (Verdict.SETTLED, None)
