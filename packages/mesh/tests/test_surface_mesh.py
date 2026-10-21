"""Meshing a TriangleSurface (ADR-005): faces must match describe_surface."""

import math
from pathlib import Path

import gmsh
import numpy as np
import pytest

from tessera.geometry import describe_surface, read_stl
from tessera.mesh import mesh_surface
from tessera.model import TET10_EDGES, TriangleSurface


def _write_stl(path: Path, build, size: float | None = None) -> Path:
    gmsh.initialize(interruptible=False)
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        build()
        gmsh.model.occ.synchronize()
        if size is not None:
            gmsh.option.setNumber("Mesh.MeshSizeMax", size)
        gmsh.model.mesh.generate(2)
        gmsh.write(str(path))
    finally:
        gmsh.finalize()
    return path


def _box() -> None:
    gmsh.model.occ.addBox(0, 0, 0, 40, 10, 10)


def _plate_with_hole() -> None:
    plate = gmsh.model.occ.addBox(0, 0, 0, 50, 50, 4)
    hole = gmsh.model.occ.addCylinder(25, 25, -1, 0, 0, 6, 10)
    gmsh.model.occ.cut([(3, plate)], [(3, hole)])


@pytest.fixture(scope="module")
def box_surface(tmp_path_factory: pytest.TempPathFactory) -> TriangleSurface:
    return read_stl(_write_stl(tmp_path_factory.mktemp("stl") / "box.stl", _box))


@pytest.fixture(scope="module")
def plate_surface(tmp_path_factory: pytest.TempPathFactory) -> TriangleSurface:
    path = _write_stl(tmp_path_factory.mktemp("stl") / "plate.stl", _plate_with_hole, size=2.0)
    return read_stl(path)


def _tet_volume(mesh) -> float:
    c = mesh.nodes[mesh.elements[:, :4]]
    vol = np.einsum("ij,ij->i", np.cross(c[:, 1] - c[:, 0], c[:, 2] - c[:, 0]), c[:, 3] - c[:, 0])
    assert (vol > 0).all()
    return float(vol.sum() / 6)


def test_midside_nodes_sit_on_their_edges(box_surface: TriangleSurface) -> None:
    mesh = mesh_surface(box_surface, part_id="box", size=5.0)
    x = mesh.nodes
    assert mesh.elements.shape[1] == 10
    for k, (i, j) in enumerate(TET10_EDGES, start=4):
        mid = 0.5 * (x[mesh.elements[:, i]] + x[mesh.elements[:, j]])
        np.testing.assert_allclose(x[mesh.elements[:, k]], mid, atol=1e-9)


def test_box_volume_and_element_set(box_surface: TriangleSurface) -> None:
    mesh = mesh_surface(box_surface, part_id="box", size=5.0)
    assert _tet_volume(mesh) == pytest.approx(4000.0, rel=1e-6)
    assert mesh.element_sets["part:box"].size == mesh.n_elements


def test_plate_volume(plate_surface: TriangleSurface) -> None:
    mesh = mesh_surface(plate_surface, part_id="plate", size=3.0)
    assert _tet_volume(mesh) == pytest.approx(50 * 50 * 4 - math.pi * 100 * 4, rel=2e-2)


def _check_surfaces_match(surface: TriangleSurface, size: float) -> None:
    faces = describe_surface(surface)
    mesh = mesh_surface(surface, part_id="p", size=size)
    assert set(mesh.surfaces) == {f"face:{f.index}" for f in faces}
    assert len(mesh.surfaces) == surface.n_faces
    x = mesh.nodes
    for f in faces:
        tris = mesh.surfaces[f"face:{f.index}"]
        cross = np.cross(x[tris[:, 1]] - x[tris[:, 0]], x[tris[:, 2]] - x[tris[:, 0]])
        areas = 0.5 * np.linalg.norm(cross, axis=1)
        centroid = (areas[:, None] * x[tris[:, :3]].mean(axis=1)).sum(axis=0) / areas.sum()
        assert areas.sum() == pytest.approx(f.area, rel=2e-2)
        assert centroid == pytest.approx(f.centroid, abs=0.5)
        nodes = mesh.node_sets[f"face:{f.index}"]
        assert set(np.unique(tris)) <= set(nodes.tolist())


def test_box_face_surfaces_match_describe_surface(box_surface: TriangleSurface) -> None:
    _check_surfaces_match(box_surface, 4.0)


def test_plate_face_surfaces_match_describe_surface(plate_surface: TriangleSurface) -> None:
    _check_surfaces_match(plate_surface, 3.0)


def test_face_node_sets_lie_on_the_face(box_surface: TriangleSurface) -> None:
    mesh = mesh_surface(box_surface, part_id="box", size=5.0)
    for f in describe_surface(box_surface):
        assert f.normal is not None
        nodes = mesh.nodes[mesh.node_sets[f"face:{f.index}"]]
        offsets = (nodes - np.array(f.centroid)) @ np.array(f.normal)
        np.testing.assert_allclose(offsets, 0.0, atol=1e-6)


def test_size_changes_the_node_count(box_surface: TriangleSurface) -> None:
    coarse = mesh_surface(box_surface, part_id="box", size=8.0).n_nodes
    fine = mesh_surface(box_surface, part_id="box", size=3.0).n_nodes
    assert fine > coarse


def test_mesh_options_do_not_leak_into_a_running_gmsh_session(
    box_surface: TriangleSurface,
) -> None:
    fresh = mesh_surface(box_surface, part_id="box").n_nodes
    gmsh.initialize(interruptible=False)  # e.g. a long-lived GUI or MCP process
    try:
        before = gmsh.option.getNumber("Mesh.MeshSizeMax")
        mesh_surface(box_surface, part_id="box", size=2.0)
        assert gmsh.option.getNumber("Mesh.MeshSizeMax") == before
        assert mesh_surface(box_surface, part_id="box").n_nodes == fresh
    finally:
        gmsh.finalize()


def test_inverted_stl_meshes_like_the_original(
    tmp_path: Path, box_surface: TriangleSurface
) -> None:
    flipped = TriangleSurface(
        vertices=box_surface.vertices,
        triangles=box_surface.triangles[:, ::-1].copy(),
        face_ids=box_surface.face_ids,
    )
    # read_stl fixes winding; mesh_surface must still produce a valid solid from either.
    mesh = mesh_surface(flipped, part_id="box", size=5.0)
    assert _tet_volume(mesh) == pytest.approx(4000.0, rel=1e-6)
    assert len(mesh.surfaces) == box_surface.n_faces


def test_size_factor_refines_an_stl_mesh(plate_surface: TriangleSurface) -> None:
    coarse = mesh_surface(plate_surface, part_id="plate", size=4.0)
    fine = mesh_surface(plate_surface, part_id="plate", size=4.0, size_factor=2**-0.5)
    # Measured 2.0x on this thin plate (a solid would be nearer 2.8x); 1.5x shows it refined.
    assert fine.n_nodes > 1.5 * coarse.n_nodes
