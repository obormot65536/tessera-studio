"""describe_edges: where two faces meet, for STEP and STL alike (ADR-006)."""

from collections import Counter
from collections.abc import Callable
from pathlib import Path

import gmsh
import pytest

import tessera.geometry as geometry
from tessera.geometry import describe_edges, describe_faces
from tessera.geometry.edges import chains
from tessera.model import EdgeInfo, FaceInfo


def _write(path: Path, build: Callable[[], object], *, stl: bool = False) -> Path:
    gmsh.initialize(interruptible=False)
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        build()
        gmsh.model.occ.synchronize()
        if stl:
            gmsh.model.mesh.generate(2)
        gmsh.write(str(path))
    finally:
        gmsh.finalize()
    return path


def _box() -> object:
    return gmsh.model.occ.addBox(0, 0, 0, 100, 10, 20)


def _by_centroid(faces: list[FaceInfo]) -> dict[tuple[float, float, float], int]:
    return {tuple(round(c, 6) for c in f.centroid): f.index for f in faces}  # type: ignore[misc]


def test_a_step_box_has_twelve_straight_edges(tmp_path: Path) -> None:
    path = _write(tmp_path / "box.step", _box)
    edges = describe_edges(path)
    assert len(edges) == 12
    assert all(e.kind == "Line" and e.curves == 1 for e in edges)
    assert Counter(round(e.length, 9) for e in edges) == {100.0: 4, 10.0: 4, 20.0: 4}
    assert [e.faces for e in edges] == sorted(e.faces for e in edges)
    assert all(a < b for a, b in (e.faces for e in edges))
    # Every face of a box meets four others.
    assert Counter(f for e in edges for f in e.faces) == {i: 4 for i in range(1, 7)}
    for e in edges:
        assert e.start is not None and e.end is not None
        middle = [(s + t) / 2 for s, t in zip(e.start, e.end, strict=True)]
        assert e.point == pytest.approx(middle)


def test_the_edge_between_two_box_faces_lies_where_they_meet(tmp_path: Path) -> None:
    path = _write(tmp_path / "box.step", _box)
    faces = describe_faces(path)
    bottom = next(f.index for f in faces if f.normal == pytest.approx((0, 0, -1)))
    end = next(f.index for f in faces if f.normal == pytest.approx((-1, 0, 0)))
    edge = next(e for e in describe_edges(path) if e.faces == tuple(sorted((bottom, end))))
    assert (edge.start, edge.end) == (pytest.approx((0, 0, 0)), pytest.approx((0, 10, 0)))
    assert edge.length == pytest.approx(10.0)


def test_a_step_cylinder_has_two_closed_circles_and_no_seam(tmp_path: Path) -> None:
    path = _write(tmp_path / "rod.step", lambda: gmsh.model.occ.addCylinder(0, 0, 0, 0, 0, 30, 5))
    edges = describe_edges(path)
    assert len(edges) == 2, "the seam on the curved face is not an edge"
    for e in edges:
        assert e.kind == "Circle" and e.curves == 1
        assert (e.start, e.end) == (None, None)
        assert e.length == pytest.approx(2 * 3.141592653589793 * 5)
        assert (e.point[0] ** 2 + e.point[1] ** 2) ** 0.5 == pytest.approx(5.0)
    assert sorted(round(e.point[2], 9) for e in edges) == [0.0, 30.0]


def test_an_stl_box_has_the_edges_of_the_step_box(tmp_path: Path) -> None:
    step = _write(tmp_path / "box.step", _box)
    stl = _write(tmp_path / "box.stl", _box, stl=True)
    # Face numbers differ between the formats; match faces by centroid.
    step_face = _by_centroid(describe_faces(step))
    stl_to_step = {i: step_face[c] for c, i in _by_centroid(describe_faces(stl)).items()}

    def key(e: EdgeInfo, faces: dict[int, int]) -> tuple[object, ...]:
        pair = tuple(sorted(faces.get(f, f) for f in e.faces))
        return pair, e.kind, e.curves, round(e.length, 6), e.start, e.end

    from_stl = sorted(key(e, stl_to_step) for e in describe_edges(stl))
    from_step = sorted(key(e, {}) for e in describe_edges(step))
    assert len(from_stl) == 12
    assert [k[:4] for k in from_stl] == [k[:4] for k in from_step]
    for a, b in zip(from_stl, from_step, strict=True):
        assert a[4] == pytest.approx(b[4]) and a[5] == pytest.approx(b[5])


def test_an_stl_cylinder_rim_is_a_closed_curved_edge(tmp_path: Path) -> None:
    def rod() -> object:
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 24)
        return gmsh.model.occ.addCylinder(0, 0, 0, 0, 0, 30, 5)

    edges = describe_edges(_write(tmp_path / "rod.stl", rod, stl=True))
    assert len(edges) == 2
    for e in edges:
        assert (e.kind, e.curves, e.start, e.end) == ("Curved", 1, None, None)
        assert e.length == pytest.approx(2 * 3.141592653589793 * 5, rel=0.02)


def test_faces_and_edges_are_read_together_once(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    reads: list[Path] = []
    parse = geometry.describe_step

    def counting(path: Path) -> tuple[list[FaceInfo], list[EdgeInfo]]:
        reads.append(path)
        return parse(path)

    monkeypatch.setattr(geometry, "describe_step", counting)
    path = _write(tmp_path / "box.step", _box)
    describe_faces(path)
    describe_edges(path).clear()  # a result served from the cache
    assert len(describe_edges(path)) == 12
    assert reads == [path]


def test_chains_tell_open_closed_and_separate_curves() -> None:
    assert chains([(1, 2), (2, 3)]) == [(1, 3)]
    assert chains([(1, 2), (2, 3), (3, 1)]) == [None]
    assert chains([(7, 7)]) == [None]
    assert chains([(1, 2), (5, 6)]) == [(1, 2), (5, 6)]
