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
