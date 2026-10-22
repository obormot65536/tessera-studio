"""Wall thickness at a surface node (ADR-007 Amendment 2): ray cast against the outer surface.

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
    IsotropicMaterial,
    Mesh,
    Section,
    SolveJob,
    WallSpan,
)
from tessera.post import WALL_ELEMENTS_MIN, local_edge_length, outer_surface, wall_span

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


def job(mesh: Mesh, *supports: FixedSupport) -> SolveJob:
    return SolveJob(
        mesh=mesh,
        materials=(IsotropicMaterial(name="m", youngs_modulus=1000.0, poisson_ratio=0.3),),
        sections=(Section(element_set="part:p", material="m"),),
        supports=supports,
    )


def node_at(mesh: Mesh, point: tuple[float, float, float]) -> int:
    (found,) = np.flatnonzero(np.isclose(mesh.nodes, point).all(axis=1))
    return int(found)


def span_at(solve_job: SolveJob, point: tuple[float, float, float]) -> WallSpan | None:
    return wall_span(outer_surface(solve_job), node_at(solve_job.mesh, point))


def test_the_approved_default() -> None:
    assert WALL_ELEMENTS_MIN == 5.0


def test_thin_plate_is_its_thickness_across() -> None:
    plate = tet10_grid((6, 6, 2))  # 6 x 6 x 2 mm, 1 mm cells
    span = span_at(job(plate), (3.0, 3.0, 2.0))
    assert span is not None
    assert span.thickness == pytest.approx(2.0, rel=1e-6)
    assert span.direction == pytest.approx((0.0, 0.0, -1.0))
    assert not span.mirrored
    assert span.elements == pytest.approx(2.0 / local_edge_length(plate, (3.0, 3.0, 2.0)))
    # From the bottom too, and from a midside node.
    bottom = span_at(job(plate), (3.0, 3.0, 0.0))
    assert bottom is not None and bottom.thickness == pytest.approx(2.0, rel=1e-6)
    assert bottom.direction == pytest.approx((0.0, 0.0, 1.0))
    middle = span_at(job(plate), (3.5, 3.0, 2.0))
    assert middle is not None and middle.thickness == pytest.approx(2.0, rel=1e-6)


def test_elements_across_scale_with_the_layers() -> None:
    # The same 6 x 6 x 2 mm plate with half-size elements: twice the elements, same thickness.
    coarse = span_at(job(tet10_grid((6, 6, 2))), (3.0, 3.0, 2.0))
    fine = span_at(job(tet10_grid((12, 12, 4), mapping=scaled(0.5))), (3.0, 3.0, 2.0))
    assert coarse is not None and fine is not None
    assert fine.thickness == pytest.approx(coarse.thickness, rel=1e-6)
    assert fine.elements == pytest.approx(2.0 * coarse.elements, rel=1e-6)


@pytest.mark.parametrize(
    ("point", "thickness", "direction"),
    [
        ((0.0, 0.0, 0.0), 2.0, (0.0, 0.0, 1.0)),  # corner: x 8, y 4, z 2 mm
        ((4.0, 0.0, 0.0), 2.0, (0.0, 0.0, 1.0)),  # edge between the y = 0 and z = 0 faces
        ((0.0, 0.0, 1.0), 4.0, (0.0, 1.0, 0.0)),  # edge between x = 0 (8 mm) and y = 0 (4 mm)
        ((8.0, 4.0, 1.0), 4.0, (0.0, -1.0, 0.0)),
    ],
)
def test_edges_and_corners_take_the_thinnest_direction(
    point: tuple[float, float, float], thickness: float, direction: tuple[float, float, float]
) -> None:
    span = span_at(job(tet10_grid((8, 4, 2))), point)
    assert span is not None
    assert span.thickness == pytest.approx(thickness, rel=1e-6)
    assert span.direction == pytest.approx(direction)


def test_l_section_inside_corner_is_the_thinner_arm() -> None:
    # 8 x 2 x 10 mm minus x > 4, z > 6: an upright arm 4 mm thick (x < 4) standing on a base
    # 6 mm thick (z < 6). The inside edge runs along y at x = 4, z = 6.
    keep = np.ones((8, 2, 10), dtype=bool)
    keep[4:, :, 6:] = False
    mesh = tet10_grid((8, 2, 10), keep=keep)
    corner = span_at(job(mesh), (4.0, 1.0, 6.0))
    assert corner is not None
    assert corner.thickness == pytest.approx(4.0, rel=1e-6)
    assert corner.direction == pytest.approx((-1.0, 0.0, 0.0))
    base = span_at(job(mesh), (6.0, 1.0, 6.0))
    assert base is not None and base.thickness == pytest.approx(6.0, rel=1e-6)


def _plate_on(*supports: FixedSupport) -> SolveJob:
    """The 6 x 6 x 2 mm plate with node set ``support:sym`` on its bottom face z = 0."""
    plate = tet10_grid((6, 6, 2))
    bottom = np.flatnonzero(np.isclose(plate.nodes[:, 2], 0.0))
    mesh = Mesh(
        nodes=plate.nodes,
        elements=plate.elements,
        element_sets=plate.element_sets,
        node_sets={"support:sym": bottom},
    )
    return job(mesh, *supports)


def test_a_symmetry_plane_doubles_the_wall() -> None:
    # The bottom holds z only: an exact mirror plane, so the real plate is 4 mm thick.
    span = span_at(_plate_on(FixedSupport(node_set="support:sym", dofs=(3,))), (3.0, 3.0, 2.0))
    assert span is not None
    assert span.thickness == pytest.approx(4.0, rel=1e-6)
    assert span.mirrored
    assert span.elements == pytest.approx(
        4.0 / local_edge_length(_plate_on().mesh, (3.0, 3.0, 2.0))
    )


@pytest.mark.parametrize(
    "support",
    [
        FixedSupport(node_set="support:sym", dofs=(1,)),  # a roller along x: no mirror in z
        FixedSupport(node_set="support:sym", dofs=(1, 3)),  # two directions: a clamp
        FixedSupport(node_set="support:sym", dofs=(1, 2, 3)),
        FixedSupport(node_set="support:sym", dofs=(3,), line=True),
    ],
)
def test_other_supports_do_not_mirror(support: FixedSupport) -> None:
    span = span_at(_plate_on(support), (3.0, 3.0, 2.0))
    assert span is not None
    assert span.thickness == pytest.approx(2.0, rel=1e-6)
    assert not span.mirrored


def annular_sector(segments: int) -> Mapping:
    """Grid (r, theta, z) -> quarter of a tube, bore radius 5 mm, wall 3 mm."""

    def mapping(p: NDArray[np.float64]) -> NDArray[np.float64]:
        r = 5.0 + p[:, 0]
        theta = p[:, 1] * (np.pi / 2) / segments
        return np.stack([r * np.cos(theta), r * np.sin(theta), p[:, 2]], axis=1)

    return mapping


@pytest.mark.parametrize("theta", [45.0, 37.5, 52.5])  # at a corner node, at midside nodes
def test_tube_bore_is_the_wall_within_one_percent(theta: float) -> None:
    mesh = tet10_grid((3, 6, 2), mapping=annular_sector(6))
    t = np.radians(theta)
    span = span_at(job(mesh), (5.0 * np.cos(t), 5.0 * np.sin(t), 1.0))
    assert span is not None
    assert span.thickness == pytest.approx(3.0, rel=0.01)
    assert span.direction == pytest.approx((np.cos(t), np.sin(t), 0.0), abs=0.02)


def test_interior_node_has_no_wall() -> None:
    cube = tet10_grid((4, 4, 4))
    assert span_at(job(cube), (2.0, 2.0, 2.0)) is None
