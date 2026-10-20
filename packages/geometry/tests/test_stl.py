"""STL import and Tessera-defined faces (ADR-005)."""

import math
from collections.abc import Callable
from pathlib import Path

import gmsh
import numpy as np
import pytest

from tessera.geometry import describe_faces, describe_surface, read_stl
from tessera.model import InputError, TriangleSurface

BOX_CENTRE = (20.0, 5.0, 5.0)


def _write_stl(
    path: Path, build: Callable[[], None], *, binary: bool = False, size: float | None = None
) -> Path:
    gmsh.initialize(interruptible=False)
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        build()
        gmsh.model.occ.synchronize()
        if size is not None:
            gmsh.option.setNumber("Mesh.MeshSizeMax", size)
        gmsh.option.setNumber("Mesh.Binary", 1 if binary else 0)
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


def _hex_prism_triangles() -> np.ndarray:
    """Hexagonal prism (circumradius 10, height 10), counter-clockwise seen from outside."""
    ang = [k * math.pi / 3 for k in range(6)]
    lo = [np.array([10 * math.cos(a), 10 * math.sin(a), 0.0]) for a in ang]
    hi = [p + np.array([0.0, 0.0, 10.0]) for p in lo]
    c_lo, c_hi = np.array([0.0, 0.0, 0.0]), np.array([0.0, 0.0, 10.0])
    tris = []
    for k in range(6):
        n = (k + 1) % 6
        tris.append([c_lo, lo[n], lo[k]])  # bottom, normal -z
        tris.append([c_hi, hi[k], hi[n]])  # top, normal +z
        tris.append([lo[k], lo[n], hi[n]])
        tris.append([lo[k], hi[n], hi[k]])
    return np.array(tris)


def _cylinder() -> None:
    gmsh.model.occ.addCylinder(0, 0, 0, 0, 0, 20, 10)


def _two_boxes() -> None:
    gmsh.model.occ.addBox(0, 0, 0, 10, 10, 10)
    gmsh.model.occ.addBox(50, 0, 0, 10, 10, 10)


def _read_ascii(path: Path) -> np.ndarray:
    """Triangles (M, 3, 3) from an ASCII STL."""
    v = [
        [float(t) for t in line.split()[1:4]]
        for line in path.read_text().splitlines()
        if line.strip().startswith("vertex")
    ]
    return np.array(v).reshape(-1, 3, 3)


def _write_ascii(path: Path, tris: np.ndarray) -> Path:
    lines = ["solid x"]
    for t in tris:
        n = np.cross(t[1] - t[0], t[2] - t[0])
        n = n / (np.linalg.norm(n) or 1.0)
        lines.append(f"facet normal {n[0]:.9g} {n[1]:.9g} {n[2]:.9g}")
        lines.append("outer loop")
        lines += [f"vertex {p[0]:.9g} {p[1]:.9g} {p[2]:.9g}" for p in t]
        lines += ["endloop", "endfacet"]
    lines.append("endsolid x")
    path.write_text("\n".join(lines) + "\n")
    return path


@pytest.fixture(scope="module")
def box_ascii(tmp_path_factory: pytest.TempPathFactory) -> Path:
    return _write_stl(tmp_path_factory.mktemp("stl") / "box.stl", _box)


@pytest.fixture(scope="module")
def box_binary(tmp_path_factory: pytest.TempPathFactory) -> Path:
    return _write_stl(tmp_path_factory.mktemp("stl") / "box_bin.stl", _box, binary=True)


@pytest.fixture(scope="module")
def plate_stl(tmp_path_factory: pytest.TempPathFactory) -> Path:
    return _write_stl(tmp_path_factory.mktemp("stl") / "plate.stl", _plate_with_hole, size=2.0)


@pytest.fixture(scope="module")
def cylinder_stl(tmp_path_factory: pytest.TempPathFactory) -> Path:
    return _write_stl(tmp_path_factory.mktemp("stl") / "cyl.stl", _cylinder, size=1.0)


def test_read_stl_returns_a_triangle_surface(box_ascii: Path) -> None:
    surface = read_stl(box_ascii)
    assert isinstance(surface, TriangleSurface)
    assert surface.n_faces == 6
    assert set(np.unique(surface.face_ids)) == set(range(1, 7))


def test_box_has_six_planar_faces_with_outward_normals(box_ascii: Path) -> None:
    faces = describe_surface(read_stl(box_ascii))
    assert [f.index for f in faces] == [1, 2, 3, 4, 5, 6]
    assert all(f.kind == "Plane" for f in faces)
    for f in faces:
        assert f.normal is not None
        assert np.linalg.norm(f.normal) == pytest.approx(1.0)
        d = np.array(f.centroid) - np.array(BOX_CENTRE)
        assert float(d @ np.array(f.normal)) > 0
    assert sorted(round(f.area, 6) for f in faces) == pytest.approx([100, 100, 400, 400, 400, 400])
    assert sum(f.area for f in faces) == pytest.approx(1800.0)
    root = next(f for f in faces if f.centroid[0] == pytest.approx(0.0))
    assert root.normal == pytest.approx((-1.0, 0.0, 0.0))
    assert root.area == pytest.approx(100.0)
    assert root.centroid == pytest.approx((0.0, 5.0, 5.0))


def test_binary_and_ascii_give_identical_faces(box_ascii: Path, box_binary: Path) -> None:
    a = describe_surface(read_stl(box_ascii))
    b = describe_surface(read_stl(box_binary))
    assert len(a) == len(b) == 6
    for fa, fb in zip(a, b, strict=True):
        assert fa.index == fb.index
        assert fa.kind == fb.kind
        assert fa.area == pytest.approx(fb.area, rel=1e-5)
        assert fa.centroid == pytest.approx(fb.centroid, abs=1e-4)
        assert fa.normal == pytest.approx(fb.normal, abs=1e-5)


def test_inverted_winding_is_fixed(box_ascii: Path, tmp_path: Path) -> None:
    tris = _read_ascii(box_ascii)
    inverted = _write_ascii(tmp_path / "inv.stl", tris[:, ::-1, :])
    original = describe_surface(read_stl(box_ascii))
    fixed = describe_surface(read_stl(inverted))
    assert len(fixed) == 6
    for f in fixed:
        assert f.normal is not None
        d = np.array(f.centroid) - np.array(BOX_CENTRE)
        assert float(d @ np.array(f.normal)) > 0
    # Same set of (normal, area), whatever the face numbering.
    key = lambda f: (*np.round(f.normal, 6), round(f.area, 6))  # noqa: E731
    assert sorted(key(f) for f in fixed) == pytest.approx(sorted(key(f) for f in original))


def test_inverted_triangles_point_outward(box_ascii: Path, tmp_path: Path) -> None:
    tris = _read_ascii(box_ascii)
    inverted = _write_ascii(tmp_path / "inv.stl", tris[:, ::-1, :])
    s = read_stl(inverted)
    t = s.vertices[s.triangles]
    normals = np.cross(t[:, 1] - t[:, 0], t[:, 2] - t[:, 0])
    centroids = t.mean(axis=1)
    assert (np.einsum("ij,ij->i", normals, centroids - np.array(BOX_CENTRE)) > 0).all()


def test_plate_with_hole_has_six_planes_and_a_curved_wall(plate_stl: Path) -> None:
    faces = describe_surface(read_stl(plate_stl))
    assert len(faces) == 7
    assert sum(f.kind == "Plane" for f in faces) == 6
    (wall,) = [f for f in faces if f.kind == "Curved"]
    assert wall.normal is None
    assert wall.area == pytest.approx(math.pi * 20 * 4, rel=1e-2)
    assert wall.centroid[:2] == pytest.approx((25.0, 25.0), abs=0.1)
    big = [f for f in faces if f.kind == "Plane" and f.area > 1000]
    assert len(big) == 2
    for f in big:
        assert f.area == pytest.approx(50 * 50 - math.pi * 10**2, rel=1e-2)
    top = next(f for f in big if f.centroid[2] == pytest.approx(4.0))
    bottom = next(f for f in big if f.centroid[2] == pytest.approx(0.0))
    assert top.normal == pytest.approx((0.0, 0.0, 1.0))
    assert bottom.normal == pytest.approx((0.0, 0.0, -1.0))


def test_hexagonal_prism_splits_into_eight_faces(tmp_path: Path) -> None:
    path = _write_ascii(tmp_path / "hex.stl", _hex_prism_triangles())
    faces = describe_surface(read_stl(path, feature_angle_deg=40.0))
    assert len(faces) == 8
    assert all(f.kind == "Plane" for f in faces)
    side = 10 * 10  # edge length 10 (circumradius 10), height 10
    assert sorted(round(f.area, 3) for f in faces)[:6] == pytest.approx([side] * 6)


def test_hexagonal_prism_merges_into_three_faces_above_60_degrees(tmp_path: Path) -> None:
    path = _write_ascii(tmp_path / "hex.stl", _hex_prism_triangles())
    faces = describe_surface(read_stl(path, feature_angle_deg=70.0))
    # The 60 degree wall edges no longer separate faces; the 90 degree caps still do.
    assert len(faces) == 3
    assert sorted(f.kind for f in faces) == ["Curved", "Plane", "Plane"]


def test_fine_cylinder_wall_is_one_curved_face(cylinder_stl: Path) -> None:
    faces = describe_surface(read_stl(cylinder_stl))
    assert len(faces) == 3
    (wall,) = [f for f in faces if f.kind == "Curved"]
    assert wall.area == pytest.approx(2 * math.pi * 10 * 20, rel=1e-2)


def test_lower_feature_angle_splits_the_cylinder_wall(cylinder_stl: Path) -> None:
    default = read_stl(cylinder_stl).n_faces
    split = read_stl(cylinder_stl, feature_angle_deg=1.0).n_faces
    assert default == 3
    assert split > default


def test_face_ids_are_deterministic(plate_stl: Path) -> None:
    a = read_stl(plate_stl)
    b = read_stl(plate_stl)
    np.testing.assert_array_equal(a.face_ids, b.face_ids)
    np.testing.assert_array_equal(a.triangles, b.triangles)
    np.testing.assert_allclose(a.vertices, b.vertices)


def test_describe_faces_dispatches_on_stl_suffix(box_ascii: Path, tmp_path: Path) -> None:
    upper = tmp_path / "BOX.STL"
    upper.write_bytes(box_ascii.read_bytes())
    for path in (box_ascii, upper):
        faces = describe_faces(path)
        assert len(faces) == 6
        assert faces == describe_surface(read_stl(path))


def test_describe_faces_still_handles_step(tmp_path: Path) -> None:
    path = tmp_path / "box.step"
    gmsh.initialize(interruptible=False)
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.occ.addBox(0, 0, 0, 40, 10, 10)
    gmsh.model.occ.synchronize()
    gmsh.write(str(path))
    gmsh.finalize()
    assert len(describe_faces(path)) == 6


def test_missing_file(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        read_stl(tmp_path / "nope.stl")
    with pytest.raises(FileNotFoundError):
        describe_faces(tmp_path / "nope.stl")


@pytest.mark.parametrize(
    "content", [b"", b"this is not an STL file at all", b"solid x\nendsolid x\n"]
)
def test_unreadable_or_empty_file_is_an_input_error(tmp_path: Path, content: bytes) -> None:
    path = tmp_path / "junk.stl"
    path.write_bytes(content)
    with pytest.raises(InputError):
        read_stl(path)


def test_open_surface_is_rejected(box_ascii: Path, tmp_path: Path) -> None:
    tris = _read_ascii(box_ascii)
    path = _write_ascii(tmp_path / "open.stl", tris[1:])  # drop one triangle
    with pytest.raises(InputError, match=r"closed|watertight"):
        read_stl(path)


def test_two_bodies_are_rejected(tmp_path: Path) -> None:
    path = _write_stl(tmp_path / "two.stl", _two_boxes)
    with pytest.raises(InputError, match=r"bodies|solid|body"):
        read_stl(path)


def _split_with_sliver(tris: np.ndarray) -> np.ndarray:
    """Split triangle (a, b, c) at m, the midpoint of a-b, and close the gap with the
    collinear triangle (a, b, m). The result is closed and consistently wound, with one
    zero-area triangle whose three corners are distinct."""
    a, b, c = tris[0]
    m = 0.5 * (a + b)
    pieces = np.array([[a, m, c], [m, b, c], [a, b, m]])
    return np.concatenate([pieces, tris[1:]])


def test_zero_area_sliver_triangles_are_rejected_clearly(box_ascii: Path, tmp_path: Path) -> None:
    path = _write_ascii(tmp_path / "sliver.stl", _split_with_sliver(_read_ascii(box_ascii)))
    with pytest.raises(InputError, match="zero-area"):
        read_stl(path)


def test_needle_triangles_with_repeated_corners_are_dropped(
    box_ascii: Path, tmp_path: Path
) -> None:
    tris = _read_ascii(box_ascii)
    a, b = tris[0][0], tris[0][1]
    needle = np.array([[a, a, b]])  # two identical corners: no area, no topology
    surface = read_stl(_write_ascii(tmp_path / "needle.stl", np.concatenate([tris, needle])))
    assert surface.n_faces == 6
    faces = describe_surface(surface)
    assert all(f.area > 0 and f.kind == "Plane" for f in faces)


def test_inconsistent_winding_is_rejected(box_ascii: Path, tmp_path: Path) -> None:
    tris = _read_ascii(box_ascii)
    tris[0] = tris[0][::-1]  # flip one triangle only
    with pytest.raises(InputError, match="winding"):
        read_stl(_write_ascii(tmp_path / "mixed.stl", tris))


@pytest.mark.parametrize("angle", [0.0, 0.5, 181.0, -10.0])
def test_feature_angle_outside_1_to_180_degrees_is_rejected(box_ascii: Path, angle: float) -> None:
    with pytest.raises(InputError, match="feature angle"):
        read_stl(box_ascii, feature_angle_deg=angle)


def test_open_surface_is_reported_as_open_not_as_several_bodies(
    box_ascii: Path, tmp_path: Path
) -> None:
    tris = _read_ascii(box_ascii)
    with pytest.raises(InputError, match="not a closed solid"):
        read_stl(_write_ascii(tmp_path / "open2.stl", tris[2:]))


def test_large_strip_of_triangles_is_segmented_quickly(tmp_path: Path) -> None:
    import time

    # A fine cylinder written in shuffled triangle order: one long chain of wall triangles.
    path = _write_stl(tmp_path / "fine_cyl.stl", _cylinder, size=0.4)
    tris = _read_ascii(path)
    rng = np.random.default_rng(0)
    shuffled = _write_ascii(tmp_path / "shuffled.stl", tris[rng.permutation(len(tris))])
    start = time.perf_counter()
    surface = read_stl(shuffled)
    assert time.perf_counter() - start < 3.0, f"{len(tris)} triangles took too long"
    assert surface.n_faces == 3
