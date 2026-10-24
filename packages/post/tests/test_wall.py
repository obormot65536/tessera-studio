"""Wall thickness at a node (ADR-007 Amendment 2): the largest ball tangent at the node.

The meshes are small structured tet10 grids built here with numpy, so no mesher or solver runs.
"""

import itertools
from collections.abc import Callable

import numpy as np
import pytest
from numpy.typing import NDArray

from tessera.model import TET10_EDGES, Mesh, WallSpan
from tessera.post import (
    WALL_ELEMENTS_MIN,
    local_edge_length,
    outer_surface,
    samples_on,
    wall_span,
)
from tessera.post.flags import _outer_triangles  # pyright: ignore[reportPrivateUsage]
from tessera.post.wall import _surface  # pyright: ignore[reportPrivateUsage]

Mapping = Callable[[NDArray[np.float64]], NDArray[np.float64]]


def tet10_grid(
    cells: tuple[int, int, int],
    keep: NDArray[np.bool_] | None = None,
    mapping: Mapping | None = None,
) -> Mesh:
    """Quadratic tets on a grid of unit cells, six (Kuhn) tets per cell; ``keep`` selects cells
    and ``mapping`` maps grid coordinates to mm, midside nodes included (curved faces)."""
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
    used, corners = np.unique(np.array(tets, dtype=np.int64), return_inverse=True)
    corners = corners.reshape(-1, 4)
    param = grid[used].astype(np.float64)
    xyz = mapping(param)
    a, b, c, d = (xyz[corners[:, n]] for n in range(4))
    negative = np.einsum("ij,ij->i", np.cross(b - a, c - a), d - a) < 0
    corners[negative] = corners[negative][:, [0, 2, 1, 3]]
    pairs = np.sort(corners[:, np.array(TET10_EDGES)], axis=2)
    edges, inverse = np.unique(pairs.reshape(-1, 2), axis=0, return_inverse=True)
    nodes = np.vstack([xyz, mapping(param[edges].mean(axis=1))])
    elements = np.hstack([corners, len(xyz) + inverse.reshape(-1, 6)])
    return Mesh(nodes=nodes, elements=elements, element_sets={"part:p": np.arange(len(elements))})


def scaled(factor: float) -> Mapping:
    return lambda p: p * factor


def node_at(mesh: Mesh, point: tuple[float, float, float]) -> int:
    (found,) = np.flatnonzero(np.isclose(mesh.nodes, point).all(axis=1))
    return int(found)


def span_at(mesh: Mesh, point: tuple[float, float, float]) -> WallSpan | None:
    return wall_span(outer_surface(mesh), node_at(mesh, point))


def test_the_approved_default() -> None:
    assert WALL_ELEMENTS_MIN == 5.0


def test_thin_plate_is_its_thickness_across() -> None:
    plate = tet10_grid((6, 6, 2))  # 6 x 6 x 2 mm, 1 mm cells
    span = span_at(plate, (3.0, 3.0, 2.0))
    assert span is not None
    assert span.thickness == pytest.approx(2.0, rel=1e-6)
    assert span.direction == pytest.approx((0.0, 0.0, -1.0))
    assert span.elements == pytest.approx(2.0 / local_edge_length(plate, (3.0, 3.0, 2.0)))
    # From the bottom too, and from a midside node.
    bottom = span_at(plate, (3.0, 3.0, 0.0))
    assert bottom is not None and bottom.thickness == pytest.approx(2.0, rel=1e-6)
    assert bottom.direction == pytest.approx((0.0, 0.0, 1.0))
    middle = span_at(plate, (3.5, 3.0, 2.0))
    assert middle is not None and middle.thickness == pytest.approx(2.0, rel=1e-6)


def test_elements_across_scale_with_the_layers() -> None:
    # The same 6 x 6 x 2 mm plate with half-size elements: twice the elements, same thickness.
    coarse = span_at(tet10_grid((6, 6, 2)), (3.0, 3.0, 2.0))
    fine = span_at(tet10_grid((12, 12, 4), mapping=scaled(0.5)), (3.0, 3.0, 2.0))
    assert coarse is not None and fine is not None
    assert fine.thickness == pytest.approx(coarse.thickness, rel=1e-6)
    assert fine.elements == pytest.approx(2.0 * coarse.elements, rel=1e-6)


@pytest.mark.parametrize(
    "point",
    [
        (0.0, 0.0, 0.0),  # corner: 8 mm along x, 4 along y, 6 along z
        (4.0, 0.0, 0.0),  # edge between the y = 0 and z = 0 faces
        (0.0, 0.0, 3.0),  # edge between the x = 0 and y = 0 faces
        (8.0, 4.0, 3.0),
        (0.0, 2.0, 3.0),  # middle of the x = 0 end face: the wall, not the 8 mm length
    ],
)
def test_edges_corners_and_end_faces_take_the_thinnest_wall(
    point: tuple[float, float, float],
) -> None:
    span = span_at(tet10_grid((8, 4, 6)), point)
    assert span is not None
    assert span.thickness == pytest.approx(4.0, rel=1e-6)


def test_corner_takes_the_smallest_ball_over_its_faces() -> None:
    # 8 x 4 x 2: the largest ball tangent to the x or the y face at the corner is 4 mm across,
    # the one tangent to the z face 2 mm; the wall is the smallest, along that face's normal.
    span = span_at(tet10_grid((8, 4, 2)), (0.0, 0.0, 0.0))
    assert span is not None
    assert span.thickness == pytest.approx(2.0, rel=1e-6)
    assert span.direction == pytest.approx((0.0, 0.0, 1.0))


def test_clamped_end_face_of_a_long_bar_is_the_section_not_the_length() -> None:
    # 40 x 10 x 10 mm bar: a ray cast from the end face reads the 40 mm length.
    bar = tet10_grid((20, 5, 5), mapping=scaled(2.0))
    surface = outer_surface(bar)
    for point in [(0.0, 5.0, 5.0), (0.0, 0.0, 0.0), (0.0, 0.0, 5.0)]:
        span = wall_span(surface, node_at(bar, point))
        assert span is not None
        assert span.thickness == pytest.approx(10.0, rel=1e-6), point
    # A ball tangent 2 mm from the bottom edge touches the bottom face first.
    off_centre = wall_span(surface, node_at(bar, (0.0, 4.0, 2.0)))
    assert off_centre is not None and off_centre.thickness == pytest.approx(4.0, rel=1e-6)


def test_corner_of_a_patch_in_the_middle_of_an_end_face() -> None:
    # A 20 x 20 x 60 mm post: at (14, 14, 0), where an 8 mm patch in the middle of the end face
    # would have its corner, the largest ball tangent there touches the sides 6 mm away. A ray
    # cast reads the 60 mm length.
    post = tet10_grid((10, 10, 30), mapping=scaled(2.0))
    patch_corner = span_at(post, (14.0, 14.0, 0.0))
    assert patch_corner is not None
    assert patch_corner.thickness == pytest.approx(12.0, rel=1e-6)
    centre = span_at(post, (10.0, 10.0, 0.0))
    assert centre is not None and centre.thickness == pytest.approx(20.0, rel=1e-6)


def test_l_section_inside_corner_is_the_thinner_arm() -> None:
    # 12 x 8 x 10 mm minus x > 4, z > 6: an upright arm 4 mm thick (x < 4) standing on a base
    # 6 mm thick (z < 6). The inside edge runs along y at x = 4, z = 6; it must not collapse.
    keep = np.ones((12, 8, 10), dtype=bool)
    keep[4:, :, 6:] = False
    mesh = tet10_grid((12, 8, 10), keep=keep)
    surface = outer_surface(mesh)
    corner = wall_span(surface, node_at(mesh, (4.0, 4.0, 6.0)))
    assert corner is not None
    assert corner.thickness == pytest.approx(4.0, rel=1e-6)
    assert corner.direction == pytest.approx((-1.0, 0.0, 0.0))
    base = wall_span(surface, node_at(mesh, (8.0, 4.0, 6.0)))
    assert base is not None and base.thickness == pytest.approx(6.0, rel=1e-6)


def annular_sector(segments: int) -> Mapping:
    """Grid (r, theta, z) -> quarter of a tube, bore radius 5 mm, wall 3 mm."""

    def mapping(p: NDArray[np.float64]) -> NDArray[np.float64]:
        r = 5.0 + p[:, 0]
        theta = p[:, 1] * (np.pi / 2) / segments
        return np.stack([r * np.cos(theta), r * np.sin(theta), p[:, 2]], axis=1)

    return mapping


@pytest.mark.parametrize("theta", [45.0, 37.5, 52.5])  # at a corner node, at midside nodes
def test_tube_bore_is_the_wall_within_one_percent(theta: float) -> None:
    mesh = tet10_grid((3, 6, 8), mapping=annular_sector(6))
    t = np.radians(theta)
    span = span_at(mesh, (5.0 * np.cos(t), 5.0 * np.sin(t), 4.0))
    assert span is not None
    assert span.thickness == pytest.approx(3.0, rel=0.01)
    assert span.direction == pytest.approx((np.cos(t), np.sin(t), 0.0), abs=0.02)


def test_interior_node_is_measured_at_the_nearest_surface_node() -> None:
    # 6 x 6 x 4 block: the node (3, 3, 1) is 1 mm under the bottom face, whose wall is 4 mm.
    block = tet10_grid((6, 6, 4))
    span = span_at(block, (3.0, 3.0, 1.0))
    assert span is not None
    assert span.thickness == pytest.approx(4.0, rel=1e-6)
    assert span.direction == pytest.approx((0.0, 0.0, 1.0))
    assert span.elements == pytest.approx(4.0 / local_edge_length(block, (3.0, 3.0, 0.0)))


def test_no_wall_where_the_surface_is_open() -> None:
    # Only the top face of a plate: nothing for the ball to touch below, as on an open STL.
    plate = tet10_grid((6, 6, 2))
    tri6 = _outer_triangles(plate)
    top = np.isclose(plate.nodes[tri6[:, :3], 2], 2.0).all(axis=1)
    surface = _surface(plate, tri6[top])
    assert wall_span(surface, node_at(plate, (3.0, 3.0, 2.0))) is None


# --- Ignoring the faces of the support that flags a spot ---------------------------------


def test_a_clamped_face_is_left_out_so_the_wall_is_the_depth() -> None:
    # 40 x 10 x 10 mm bar clamped on x = 0. At (2, 5, 10), 2 mm from the clamp on the top face,
    # the ball touches the clamped face and reads 4 mm; the part continues into the fixture there,
    # so with the clamp's samples left out the wall is the 10 mm depth.
    bar = tet10_grid((20, 5, 5), mapping=scaled(2.0))
    surface = outer_surface(bar)
    node = node_at(bar, (2.0, 5.0, 10.0))
    near = wall_span(surface, node)
    assert near is not None and near.thickness == pytest.approx(4.0, rel=1e-6)
    clamp = samples_on(surface, np.flatnonzero(np.isclose(bar.nodes[:, 0], 0.0)))
    assert clamp.sum() == 4 * 2 * 5 * 5  # 4 flat triangles per TRI6, 2 per 2 mm square
    span = wall_span(surface, node, ignore=clamp)
    assert span is not None
    assert span.thickness == pytest.approx(10.0, rel=1e-6)
    assert span.direction == pytest.approx((0.0, 0.0, -1.0))


def test_a_clamped_patch_leaves_the_sides_in() -> None:
    # The 20 x 20 x 60 mm post clamped on an 8 x 8 mm patch in the middle of its end face: only
    # the patch is left out, so the ball at the patch corner still touches the sides 6 mm away.
    post = tet10_grid((10, 10, 30), mapping=scaled(2.0))
    x = post.nodes
    on_patch = np.isclose(x[:, 2], 0.0) & (np.abs(x[:, 0] - 10) <= 4) & (np.abs(x[:, 1] - 10) <= 4)
    surface = outer_surface(post)
    patch = samples_on(surface, np.flatnonzero(on_patch))
    assert patch.sum() == 4 * 2 * 4 * 4
    span = wall_span(surface, node_at(post, (14.0, 14.0, 0.0)), ignore=patch)
    assert span is not None
    assert span.thickness == pytest.approx(12.0, rel=1e-6)
