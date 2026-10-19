from pathlib import Path

import gmsh
import pytest

from tessera.geometry import describe_faces


@pytest.fixture
def box_step(tmp_path: Path) -> Path:
    path = tmp_path / "box.step"
    gmsh.initialize(interruptible=False)
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.occ.addBox(0, 0, 0, 100, 10, 20)
    gmsh.model.occ.synchronize()
    gmsh.write(str(path))
    gmsh.finalize()
    return path


def test_box_has_six_planar_faces_with_outward_normals(box_step: Path) -> None:
    faces = describe_faces(box_step)
    assert len(faces) == 6
    assert all(f.kind == "Plane" for f in faces)
    for f in faces:
        assert f.normal is not None
        # Outward: moving from the box centre to the face centroid goes along the normal.
        d = [c - m for c, m in zip(f.centroid, (50.0, 5.0, 10.0), strict=True)]
        assert sum(a * b for a, b in zip(d, f.normal, strict=True)) > 0
    root = next(f for f in faces if f.centroid[0] == pytest.approx(0.0))
    assert root.normal == pytest.approx((-1.0, 0.0, 0.0))
    assert root.area == pytest.approx(200.0)


def test_missing_file(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        describe_faces(tmp_path / "nope.step")
