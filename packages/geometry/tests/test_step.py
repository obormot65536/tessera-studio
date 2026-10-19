from pathlib import Path

import gmsh
import pytest

from tessera.geometry import describe_faces
from tessera.model import InputError


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


def test_plate_with_hole_normals_point_outward(tmp_path: Path) -> None:
    # The top/bottom face centroids lie inside the hole, off the face itself.
    path = tmp_path / "plate.step"
    gmsh.initialize(interruptible=False)
    gmsh.option.setNumber("General.Terminal", 0)
    plate = gmsh.model.occ.addBox(0, 0, 0, 50, 50, 4)
    hole = gmsh.model.occ.addCylinder(25, 25, -1, 0, 0, 6, 20)
    gmsh.model.occ.cut([(3, plate)], [(3, hole)])
    gmsh.model.occ.synchronize()
    gmsh.write(str(path))
    gmsh.finalize()

    faces = describe_faces(path)
    top = next(f for f in faces if f.kind == "Plane" and f.centroid[2] == pytest.approx(4.0))
    bottom = next(f for f in faces if f.kind == "Plane" and f.centroid[2] == pytest.approx(0.0))
    assert top.normal == pytest.approx((0.0, 0.0, 1.0))
    assert bottom.normal == pytest.approx((0.0, 0.0, -1.0))
    assert any(f.kind == "Cylinder" and f.normal is None for f in faces)


def test_unreadable_step_is_an_input_error(tmp_path: Path) -> None:
    path = tmp_path / "junk.step"
    path.write_text("this is not a STEP file")
    with pytest.raises(InputError, match="cannot read"):
        describe_faces(path)
