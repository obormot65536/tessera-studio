from pathlib import Path

import gmsh
import numpy as np
import pytest

from tessera.mesh import mesh_step
from tessera.model import TET10_EDGES, TRI6_EDGES


@pytest.fixture
def box_step(tmp_path: Path) -> Path:
    path = tmp_path / "box.step"
    gmsh.initialize(interruptible=False)
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.occ.addBox(0, 0, 0, 40, 10, 10)
    gmsh.model.occ.synchronize()
    gmsh.write(str(path))
    gmsh.finalize()
    return path


def test_c3d10_midside_nodes_sit_on_their_edges(box_step: Path) -> None:
    mesh = mesh_step(box_step, part_id="box", size=5.0)
    x = mesh.nodes
    for k, (i, j) in enumerate(TET10_EDGES, start=4):
        mid = 0.5 * (x[mesh.elements[:, i]] + x[mesh.elements[:, j]])
        np.testing.assert_allclose(x[mesh.elements[:, k]], mid, atol=1e-9)


def test_tetrahedra_have_positive_volume(box_step: Path) -> None:
    mesh = mesh_step(box_step, part_id="box", size=5.0)
    c = mesh.nodes[mesh.elements[:, :4]]
    vol = np.einsum("ij,ij->i", np.cross(c[:, 1] - c[:, 0], c[:, 2] - c[:, 0]), c[:, 3] - c[:, 0])
    assert (vol > 0).all()
    assert vol.sum() / 6 == pytest.approx(40 * 10 * 10)


def test_face_surfaces_cover_the_boundary(box_step: Path) -> None:
    mesh = mesh_step(box_step, part_id="box", size=5.0)
    assert len(mesh.surfaces) == 6
    x = mesh.nodes
    total = 0.0
    for tris in mesh.surfaces.values():
        for k, (i, j) in enumerate(TRI6_EDGES, start=3):
            np.testing.assert_allclose(
                x[tris[:, k]], 0.5 * (x[tris[:, i]] + x[tris[:, j]]), atol=1e-9
            )
        a = np.cross(x[tris[:, 1]] - x[tris[:, 0]], x[tris[:, 2]] - x[tris[:, 0]])
        total += 0.5 * np.linalg.norm(a, axis=1).sum()
    assert total == pytest.approx(2 * (40 * 10 + 40 * 10 + 10 * 10))
    assert mesh.element_sets["part:box"].size == mesh.n_elements


def test_mesh_options_do_not_leak_into_a_running_gmsh_session(box_step: Path) -> None:
    fresh = mesh_step(box_step, part_id="box").n_nodes
    gmsh.initialize(interruptible=False)  # e.g. a long-lived GUI or MCP process
    try:
        before = gmsh.option.getNumber("Mesh.MeshSizeMax")
        mesh_step(box_step, part_id="box", size=2.0)
        assert gmsh.option.getNumber("Mesh.MeshSizeMax") == before
        assert mesh_step(box_step, part_id="box").n_nodes == fresh
    finally:
        gmsh.finalize()


@pytest.fixture(scope="module")
def plate_with_small_hole(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """100 x 40 x 3 plate with a 1 mm radius hole at (50, 20): gmsh sizes the hole by its
    points per circle, not by the size cap."""
    path = tmp_path_factory.mktemp("geom") / "plate.step"
    gmsh.initialize(interruptible=False)
    gmsh.option.setNumber("General.Terminal", 0)
    plate = gmsh.model.occ.addBox(0, 0, 0, 100, 40, 3)
    hole = gmsh.model.occ.addCylinder(50, 20, -1, 0, 0, 5, 1.0)
    gmsh.model.occ.cut([(3, plate)], [(3, hole)])
    gmsh.model.occ.synchronize()
    gmsh.write(str(path))
    gmsh.finalize()
    return path


def _edge_lengths(mesh) -> tuple[float, float]:
    """(median corner edge length of the mesh, the same for the edges on the hole's wall), mm.

    The hole's wall is the steady measure of the size there: the median of all edges near the
    hole mixes in the graded elements around it and varies by +-0.08 with tiny changes of the
    geometry (measured).
    """
    el = mesh.elements
    pairs = np.array(TET10_EDGES)
    a, b = el[:, pairs[:, 0]].ravel(), el[:, pairs[:, 1]].ravel()
    length = np.linalg.norm(mesh.nodes[a] - mesh.nodes[b], axis=1)
    radius = np.hypot(mesh.nodes[:, 0] - 50.0, mesh.nodes[:, 1] - 20.0)
    on_wall = np.abs(radius - 1.0) < 1e-3
    return float(np.median(length)), float(np.median(length[on_wall[a] & on_wall[b]]))


@pytest.mark.parametrize("size", [None, 2.5])
def test_size_factor_refines_everywhere_including_small_holes(
    plate_with_small_hole: Path, size: float | None
) -> None:
    # A two-mesh check (ADR-007) is only a check where the mesh got finer. Capping the size
    # alone leaves the hole's elements unchanged (measured 1.00x on its wall): the factor must
    # reach them (measured 0.69-0.72x for holes of radius 0.97-1.03 mm).
    coarse = _edge_lengths(mesh_step(plate_with_small_hole, part_id="p", size=size))
    fine = _edge_lengths(
        mesh_step(plate_with_small_hole, part_id="p", size=size, size_factor=2**-0.5)
    )
    assert fine[0] / coarse[0] <= 0.85, f"the mesh did not get finer: {fine[0] / coarse[0]:.3f}"
    assert fine[1] / coarse[1] <= 0.8, f"the hole did not get finer: {fine[1] / coarse[1]:.3f}"


@pytest.mark.parametrize("factor", [0.0, -0.5, 1.5])
def test_size_factor_must_refine(box_step: Path, factor: float) -> None:
    with pytest.raises(ValueError, match="size_factor"):
        mesh_step(box_step, part_id="box", size_factor=factor)
