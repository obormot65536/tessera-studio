"""Peak stress reporting (ADR-007): node flags, raw and clear peaks, face peaks.

The meshes are small structured tet10 grids built here with numpy, so no mesher or solver runs.
"""

import itertools
from collections.abc import Callable

import numpy as np
import pytest
from numpy.typing import NDArray

from tessera.model import (
    TET10_EDGES,
    FixedSupport,
    FlagKind,
    IsotropicMaterial,
    Mesh,
    Section,
    SolveJob,
    SolveResult,
    StressFlag,
)
from tessera.post import flag_nodes, summarize

Mapping = Callable[[NDArray[np.float64]], NDArray[np.float64]]
ROOT = StressFlag(FlagKind.SUPPORT, support="root")
INSIDE_EDGE = StressFlag(FlagKind.SHARP_INSIDE_EDGE)


def tet10_grid(
    cells: tuple[int, int, int],
    keep: NDArray[np.bool_] | None = None,
    mapping: Mapping | None = None,
    curved: bool = True,
) -> Mesh:
    """Quadratic tets on a grid of unit cells, each cell split into six (Kuhn) tets.

    ``keep`` selects cells; ``mapping`` maps grid coordinates to mm (default: identity). With
    ``curved`` the midside nodes are mapped too, so they lie on the curved surface; otherwise
    they sit halfway between the mapped corners.
    """
    nx, ny, nz = cells
    keep = np.ones(cells, dtype=bool) if keep is None else keep
    mapping = mapping or (lambda p: p)
    grid = np.stack(
        np.meshgrid(np.arange(nx + 1), np.arange(ny + 1), np.arange(nz + 1), indexing="ij"),
        axis=-1,
    ).reshape(-1, 3)

    def index(i: int, j: int, k: int) -> int:
        return (i * (ny + 1) + j) * (nz + 1) + k

    tets = []
    for i, j, k in zip(*np.nonzero(keep), strict=True):
        for perm in itertools.permutations(range(3)):
            v = [int(i), int(j), int(k)]
            path = [index(*v)]
            for axis in perm:
                v[axis] += 1
                path.append(index(*v))
            tets.append(path)
    corners = np.array(tets, dtype=np.int64)
    used, corners = np.unique(corners, return_inverse=True)
    corners = corners.reshape(-1, 4)
    param = grid[used].astype(np.float64)
    xyz = mapping(param)
    a, b, c, d = (xyz[corners[:, n]] for n in range(4))
    negative = np.einsum("ij,ij->i", np.cross(b - a, c - a), d - a) < 0
    corners[negative] = corners[negative][:, [0, 2, 1, 3]]

    pairs = np.sort(corners[:, np.array(TET10_EDGES)], axis=2)  # (M, 6, 2)
    edges, inverse = np.unique(pairs.reshape(-1, 2), axis=0, return_inverse=True)
    mid = mapping(param[edges].mean(axis=1)) if curved else xyz[edges].mean(axis=1)
    nodes = np.vstack([xyz, mid])
    elements = np.hstack([corners, len(xyz) + inverse.reshape(-1, 6)])
    return Mesh(
        nodes=nodes,
        elements=elements,
        element_sets={"part:p": np.arange(len(elements))},
    )


def with_sets(mesh: Mesh, **node_sets: NDArray[np.bool_]) -> Mesh:
    """``mesh`` with node sets given as per-node masks, e.g. ``face_1=mask``."""
    named = {name.replace("_", ":"): np.flatnonzero(m) for name, m in node_sets.items()}
    return Mesh(
        nodes=mesh.nodes,
        elements=mesh.elements,
        element_sets=mesh.element_sets,
        node_sets=named,
    )


def job(mesh: Mesh, *supports: FixedSupport) -> SolveJob:
    return SolveJob(
        mesh=mesh,
        materials=(IsotropicMaterial(name="m", youngs_modulus=1000.0, poisson_ratio=0.3),),
        sections=(Section(element_set="part:p", material="m"),),
        supports=supports,
    )


def uniaxial(values: NDArray[np.float64]) -> SolveResult:
    """A solution whose von Mises stress equals ``values`` (uniaxial xx)."""
    stresses = np.zeros((values.size, 6))
    stresses[:, 0] = values
    return SolveResult(
        displacements=np.zeros((values.size, 3)), stresses=stresses, solver="synthetic"
    )


def box_faces(mesh: Mesh) -> dict[str, NDArray[np.bool_]]:
    x = mesh.nodes
    lo, hi = x.min(axis=0), x.max(axis=0)
    faces = {}
    for axis in range(3):
        faces[f"face_{2 * axis + 1}"] = np.isclose(x[:, axis], lo[axis])
        faces[f"face_{2 * axis + 2}"] = np.isclose(x[:, axis], hi[axis])
    return faces


@pytest.fixture
def bar() -> Mesh:
    """10 x 3 x 3 mm bar, clamped (node set ``support:root``) at x = 0, faces 1-6."""
    mesh = tet10_grid((10, 3, 3))
    return with_sets(mesh, support_root=np.isclose(mesh.nodes[:, 0], 0.0), **box_faces(mesh))


# --- Flag distance -----------------------------------------------------------------------


def test_flag_distance_is_twice_the_median_corner_edge(bar: Mesh) -> None:
    # Kuhn grid of unit cells: edges 1 (axis), sqrt 2 (face diagonals), sqrt 3; the face
    # diagonals hold the median.
    flags = flag_nodes(job(bar))
    assert flags.edge_length == pytest.approx(np.sqrt(2.0))
    assert flags.distance == pytest.approx(2.0 * np.sqrt(2.0))


# --- Supports ----------------------------------------------------------------------------


@pytest.mark.parametrize("dofs", [(1, 2, 3), (1, 3), (2, 3)])
def test_support_holding_two_or_three_directions_flags_nodes_within_2h(
    bar: Mesh, dofs: tuple[int, ...]
) -> None:
    flags = flag_nodes(job(bar, FixedSupport(node_set="support:root", dofs=dofs)))
    x = bar.nodes[:, 0]
    d = flags.distance  # 2.83 mm
    np.testing.assert_array_equal(flags.flagged, x <= d)
    assert flags.flagged[np.isclose(x, 2.5)].all()
    assert not flags.flagged[np.isclose(x, 3.0)].any()
    near = int(np.flatnonzero(np.isclose(x, 1.0))[0])
    assert flags.at(near) == (ROOT,)
    assert flags.at(int(np.flatnonzero(np.isclose(x, 5.0))[0])) == ()


@pytest.mark.parametrize("dofs", [(1,), (2,), (3,)])
def test_support_holding_one_direction_is_never_flagged(bar: Mesh, dofs: tuple[int, ...]) -> None:
    flags = flag_nodes(job(bar, FixedSupport(node_set="support:root", dofs=dofs)))
    assert not flags.flagged.any()
    assert flags.masks == {}


def test_support_without_the_support_prefix_keeps_its_node_set_name(bar: Mesh) -> None:
    mesh = with_sets(bar, clamp=np.isclose(bar.nodes[:, 0], 0.0))
    flags = flag_nodes(job(mesh, FixedSupport(node_set="clamp")))
    assert list(flags.masks) == [StressFlag(FlagKind.SUPPORT, support="clamp")]


# --- Sharp inside edges ------------------------------------------------------------------


def l_block() -> Mesh:
    """12 x 2 x 12 mm block minus the 6 x 6 corner x > 6, z > 6: one inside edge along y."""
    keep = np.ones((12, 2, 12), dtype=bool)
    keep[6:, :, 6:] = False
    return tet10_grid((12, 2, 12), keep=keep)


def test_l_block_flags_only_near_the_inside_edge() -> None:
    mesh = l_block()
    flags = flag_nodes(job(mesh))
    x, z = mesh.nodes[:, 0], mesh.nodes[:, 2]
    # Every node has a node of the edge (x = z = 6) at the same y, so the distance to the
    # nearest edge node is the in-plane distance to the edge line.
    to_edge = np.hypot(x - 6.0, z - 6.0)
    np.testing.assert_array_equal(flags.flagged, to_edge <= flags.distance + 1e-9)
    assert list(flags.masks) == [INSIDE_EDGE]
    for cx, cz in [(0, 0), (12, 0), (12, 6), (6, 12), (0, 12)]:  # outside (convex) corners
        corner = np.isclose(x, cx) & np.isclose(z, cz)
        assert corner.any() and not flags.flagged[corner].any()


def test_box_has_no_sharp_inside_edge(bar: Mesh) -> None:
    assert not flag_nodes(job(bar)).flagged.any()


def annular_sector(segments: int) -> Callable[[NDArray[np.float64]], NDArray[np.float64]]:
    """Grid (r, theta, z) -> quarter of a thick ring, inner radius 5 mm: a concave bore."""

    def mapping(p: NDArray[np.float64]) -> NDArray[np.float64]:
        r = 5.0 + p[:, 0]
        theta = p[:, 1] * (np.pi / 2) / segments
        return np.stack([r * np.cos(theta), r * np.sin(theta), p[:, 2]], axis=1)

    return mapping


def test_smooth_curved_bore_is_not_flagged() -> None:
    # 30 degrees per element around the bore: the corners alone would turn 30 degrees, but the
    # quadratic faces follow the curve, so the surface is smooth where elements meet.
    mesh = tet10_grid((3, 3, 2), mapping=annular_sector(3))
    assert not flag_nodes(job(mesh)).flagged.any()


def test_faceted_bore_is_a_sharp_inside_edge() -> None:
    # Straight elements, 45 degrees per facet: the surface really turns inward at each facet
    # edge, as on a coarsely faceted STL hole.
    mesh = tet10_grid((3, 2, 2), mapping=annular_sector(2), curved=False)
    flags = flag_nodes(job(mesh))
    radius = np.hypot(mesh.nodes[:, 0], mesh.nodes[:, 1])
    middle_edge = np.isclose(radius, 5.0) & np.isclose(
        np.arctan2(mesh.nodes[:, 1], mesh.nodes[:, 0]), np.pi / 4
    )
    assert flags.flagged[middle_edge].all()
    assert list(flags.masks) == [INSIDE_EDGE]


# --- Peaks -------------------------------------------------------------------------------


def node_at(mesh: Mesh, point: tuple[float, float, float]) -> int:
    (found,) = np.flatnonzero(np.isclose(mesh.nodes, point).all(axis=1))
    return int(found)


def test_raw_and_clear_peak(bar: Mesh) -> None:
    near, away = node_at(bar, (0.0, 1.0, 3.0)), node_at(bar, (8.0, 1.0, 3.0))
    values = np.ones(bar.n_nodes)
    values[near], values[away] = 100.0, 50.0
    result = summarize(job(bar, FixedSupport(node_set="support:root")), uniaxial(values))

    peak = result.peak_von_mises
    assert (peak.value, peak.node, peak.flags) == (pytest.approx(100.0), near, (ROOT,))
    assert peak.location == pytest.approx((0.0, 1.0, 3.0))
    assert peak.face == 1  # on faces 1 (x = 0) and 6 (z = 3): the lowest number
    assert peak.flagged
    assert result.max_von_mises == pytest.approx(100.0)  # raw peak, as before

    clear = result.clear_peak_von_mises
    assert clear is not None
    assert (clear.value, clear.node, clear.flags, clear.flagged) == (
        pytest.approx(50.0),
        away,
        (),
        False,
    )
    assert clear.face == 6


def test_clear_peak_is_none_when_every_node_is_flagged(bar: Mesh) -> None:
    mesh = with_sets(bar, support_all=np.ones(bar.n_nodes, dtype=bool))
    result = summarize(job(mesh, FixedSupport(node_set="support:all")), uniaxial(bar.nodes[:, 0]))
    assert result.node_flags.flagged.all()
    assert result.clear_peak_von_mises is None
    assert result.peak_von_mises.value == pytest.approx(10.0)


def test_interior_peak_has_no_face(bar: Mesh) -> None:
    inside = node_at(bar, (5.0, 1.0, 1.0))
    values = np.zeros(bar.n_nodes)
    values[inside] = 7.0
    result = summarize(job(bar), uniaxial(values))
    assert result.peak_von_mises.node == inside
    assert result.peak_von_mises.face is None


def test_face_peaks_are_sorted_and_cover_every_face(bar: Mesh) -> None:
    values = np.random.default_rng(7).uniform(0.0, 100.0, bar.n_nodes)
    result = summarize(job(bar, FixedSupport(node_set="support:root")), uniaxial(values))
    peaks = result.face_peaks
    assert sorted(p.face for p in peaks if p.face is not None) == [1, 2, 3, 4, 5, 6]
    assert [p.value for p in peaks] == sorted((p.value for p in peaks), reverse=True)
    for p in peaks:
        on_face = bar.node_sets[f"face:{p.face}"]
        assert p.value == pytest.approx(values[on_face].max())
        assert p.node in on_face
        assert p.location == pytest.approx(tuple(bar.nodes[p.node]))
        assert p.flags == result.node_flags.at(p.node)
    on_surface = np.concatenate([bar.node_sets[f"face:{i}"] for i in range(1, 7)])
    assert peaks[0].value == pytest.approx(values[on_surface].max())


def test_nan_nodes_are_skipped(bar: Mesh) -> None:
    values = np.ones(bar.n_nodes)
    top = node_at(bar, (9.0, 0.0, 0.0))
    values[top] = 30.0
    values[node_at(bar, (8.0, 0.0, 0.0))] = np.nan  # a node with no result
    values[node_at(bar, (10.0, 0.0, 0.0))] = np.nan
    result = summarize(job(bar), uniaxial(values))
    assert result.max_von_mises == pytest.approx(30.0)
    assert result.peak_von_mises.node == top
    clear = result.clear_peak_von_mises
    assert clear is not None and clear.node == top
    assert all(np.isfinite(p.value) for p in result.face_peaks)
    assert np.isfinite(result.max_displacement)


def test_face_without_any_result_is_left_out(bar: Mesh) -> None:
    values = np.ones(bar.n_nodes)
    values[bar.node_sets["face:2"]] = np.nan  # x = 10
    result = summarize(job(bar), uniaxial(values))
    assert sorted(p.face for p in result.face_peaks if p.face is not None) == [1, 3, 4, 5, 6]


def test_no_result_at_all_is_an_error(bar: Mesh) -> None:
    with pytest.raises(ValueError, match="no node"):
        summarize(job(bar), uniaxial(np.full(bar.n_nodes, np.nan)))
