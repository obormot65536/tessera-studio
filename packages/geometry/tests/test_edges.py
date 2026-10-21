"""describe_edges: where two faces meet, for STEP and STL alike (ADR-006)."""

import math
from collections import Counter
from collections.abc import Callable
from pathlib import Path

import gmsh
import pytest

import tessera.geometry as geometry
from tessera.geometry import describe_edges, describe_faces
from tessera.geometry.edges import Chain, chains, ordered_ends
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


def _slotted_rod() -> object:
    """A rod with a slot cut across it: the wall meets each slot side along two lines, and
    the slot floor and roof along two arcs."""
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 24)
    rod = gmsh.model.occ.addCylinder(0, 0, 0, 0, 0, 30, 5)
    slot = gmsh.model.occ.addBox(-10, -1, 10, 20, 2, 10)
    return gmsh.model.occ.cut([(3, rod)], [(3, slot)])


def _stl_to_step_faces(stl: Path, step: Path) -> dict[int, int]:
    """Face numbers differ between the formats; match each STL face to the STEP face with
    the nearest centroid (faceting moves curved faces' centroids a little)."""
    step_faces = describe_faces(step)
    mapping = {}
    for f in describe_faces(stl):
        nearest = min(step_faces, key=lambda g: math.dist(g.centroid, f.centroid))
        mapping[f.index] = nearest.index
    assert sorted(mapping.values()) == [f.index for f in step_faces], "one-to-one"
    return mapping


def _in_step_numbers(e: EdgeInfo, faces: dict[int, int]) -> tuple[int, int]:
    a, b = sorted(faces[f] for f in e.faces)
    return a, b


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
    stl_to_step = _stl_to_step_faces(stl, step)
    from_stl = {_in_step_numbers(e, stl_to_step): e for e in describe_edges(stl)}
    from_step = {e.faces: e for e in describe_edges(step)}
    assert len(from_stl) == 12 and from_stl.keys() == from_step.keys()
    for pair, a in from_stl.items():
        b = from_step[pair]
        assert (a.kind, a.curves) == (b.kind, b.curves) == ("Line", 1)
        assert a.length == pytest.approx(b.length)
        assert a.start == pytest.approx(b.start) and a.end == pytest.approx(b.end)


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


def test_faces_meeting_along_two_curves_are_one_edge(tmp_path: Path) -> None:
    step = _write(tmp_path / "slot.step", _slotted_rod)
    stl = _write(tmp_path / "slot.stl", _slotted_rod, stl=True)
    stl_to_step = _stl_to_step_faces(stl, step)
    from_step = {e.faces: e for e in describe_edges(step)}
    from_stl = {_in_step_numbers(e, stl_to_step): e for e in describe_edges(stl)}
    assert from_stl.keys() == from_step.keys() and len(from_step) == 10
    for pair, b in from_step.items():
        a = from_stl[pair]
        assert a.curves == b.curves
        assert a.length == pytest.approx(b.length, rel=0.01)
        if b.kind == "Line":
            assert a.kind == "Line", pair
    # The wall meets each slot side along two straight 10 mm lines, both held by one pair.
    sides = [e for e in from_step.values() if e.curves == 2 and e.kind == "Line"]
    assert len(sides) == 2
    for e in sides:
        assert e.length == pytest.approx(20.0) and (e.start, e.end) == (None, None)
    # The floor and roof meet the wall along two short arcs each.
    arcs = [e for e in from_step.values() if e.curves == 2 and e.kind == "Circle"]
    assert len(arcs) == 2


def test_end_order_ignores_rounding_noise() -> None:
    a, b = (0.0, 0.0, 5.0), (1e-13, 10.0, 0.0)
    assert ordered_ends(b, a, 10.0) == ordered_ends(a, b, 10.0) == (a, b)
    assert ordered_ends((2.0, 0, 0), (1.0, 0, 0), 10.0) == ((1.0, 0, 0), (2.0, 0, 0))


def test_chains_tell_open_closed_and_separate_curves() -> None:
    assert [c.ends for c in chains([(1, 2), (2, 3)])] == [(1, 3)]
    assert chains([(1, 2), (2, 3), (3, 1)]) == [Chain(None, [1, 2, 3])]
    assert [c.ends for c in chains([(7, 7)])] == [None]
    assert chains([(1, 2), (5, 6)]) == [Chain((1, 2), [1, 2]), Chain((5, 6), [5, 6])]
