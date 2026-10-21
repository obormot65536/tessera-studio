"""describe_faces reads a file once per content, and again as soon as the content changes."""

from collections.abc import Callable
from pathlib import Path

import gmsh
import pytest

import tessera.geometry as geometry
from tessera.model import EdgeInfo, FaceInfo, InputError, TriangleSurface

_Parsed = tuple[list[FaceInfo], list[EdgeInfo]]


def _write_box_step(path: Path, length: float) -> Path:
    gmsh.initialize(interruptible=False)
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.occ.addBox(0, 0, 0, length, 10, 10)
        gmsh.model.occ.synchronize()
        gmsh.write(str(path))
    finally:
        gmsh.finalize()
    return path


@pytest.fixture
def step_reads(monkeypatch: pytest.MonkeyPatch) -> list[Path]:
    """Records every STEP file geometry actually parses."""
    reads: list[Path] = []
    parse: Callable[[Path], _Parsed] = geometry.describe_step

    def counting(path: Path) -> _Parsed:
        reads.append(path)
        return parse(path)

    monkeypatch.setattr(geometry, "describe_step", counting)
    return reads


def test_an_unchanged_file_is_parsed_once(tmp_path: Path, step_reads: list[Path]) -> None:
    path = _write_box_step(tmp_path / "box.step", 40)
    first = geometry.describe_faces(path)
    second = geometry.describe_faces(path)
    assert first == second and len(first) == 6
    assert step_reads == [path]


def test_a_changed_file_is_parsed_again(tmp_path: Path, step_reads: list[Path]) -> None:
    path = _write_box_step(tmp_path / "box.step", 40)
    before = geometry.describe_faces(path)
    _write_box_step(path, 80)  # same name, new content
    after = geometry.describe_faces(path)
    assert len(step_reads) == 2
    assert max(f.area for f in after) > max(f.area for f in before)


def test_a_copy_under_another_name_reuses_the_faces(tmp_path: Path, step_reads: list[Path]) -> None:
    path = _write_box_step(tmp_path / "box.step", 40)
    copy = tmp_path / "copy.stp"
    copy.write_bytes(path.read_bytes())
    assert geometry.describe_faces(copy) == geometry.describe_faces(path)
    assert len(step_reads) == 1


def test_callers_cannot_change_the_cached_faces(tmp_path: Path, step_reads: list[Path]) -> None:
    path = _write_box_step(tmp_path / "box.step", 40)
    geometry.describe_faces(path)
    geometry.describe_faces(path).clear()  # a result served from the cache
    assert len(geometry.describe_faces(path)) == 6


def test_errors_are_not_cached(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, step_reads: list[Path]
) -> None:
    path = _write_box_step(tmp_path / "box.step", 40)
    parse = geometry.describe_step
    failures = [InputError("cannot read it this time")]

    def flaky(p: Path) -> _Parsed:
        if failures:
            raise failures.pop()
        return parse(p)

    monkeypatch.setattr(geometry, "describe_step", flaky)
    with pytest.raises(InputError, match="this time"):
        geometry.describe_faces(path)
    assert len(geometry.describe_faces(path)) == 6, "the same bytes are read again"
    assert len(step_reads) == 1


def test_faces_of_a_file_changed_during_the_read_are_not_kept(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, step_reads: list[Path]
) -> None:
    path = _write_box_step(tmp_path / "box.step", 40)
    original = path.read_bytes()
    parse = geometry.describe_step

    def edited_meanwhile(p: Path) -> _Parsed:
        found = parse(p)
        _write_box_step(p, 80)
        return found

    monkeypatch.setattr(geometry, "describe_step", edited_meanwhile)
    geometry.describe_faces(path)
    monkeypatch.setattr(geometry, "describe_step", parse)
    path.write_bytes(original)
    geometry.describe_faces(path)
    assert len(step_reads) == 2, "the first result was not remembered for the old content"


def test_the_least_recently_used_file_is_forgotten_first(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, step_reads: list[Path]
) -> None:
    monkeypatch.setattr(geometry, "FACE_CACHE_SIZE", 2)
    a, b, c = (_write_box_step(tmp_path / f"{n}.step", n) for n in (10, 20, 30))
    for path in (a, b, a, c):  # c pushes out b, the least recently used
        geometry.describe_faces(path)
    geometry.describe_faces(a)
    geometry.describe_faces(b)
    assert step_reads == [a, b, c, b]


def test_stl_faces_are_remembered_and_follow_the_bytes(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    reads: list[Path] = []
    read = geometry.read_stl

    def counting(path: Path) -> TriangleSurface:
        reads.append(path)
        return read(path)

    monkeypatch.setattr(geometry, "read_stl", counting)
    step = _write_box_step(tmp_path / "box.step", 40)
    stl = tmp_path / "box.stl"
    _step_to_stl(step, stl)
    first = geometry.describe_faces(stl)
    assert geometry.describe_faces(stl) == first and len(reads) == 1
    _step_to_stl(_write_box_step(step, 80), stl)
    assert geometry.describe_faces(stl) != first and len(reads) == 2


def test_an_unreadable_file_is_an_input_error(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    path = _write_box_step(tmp_path / "box.step", 40)

    def denied(*_: object, **__: object) -> None:
        raise PermissionError("permission denied")

    monkeypatch.setattr(Path, "open", denied)
    with pytest.raises(InputError, match="cannot read"):
        geometry.describe_faces(path)


def _step_to_stl(step: Path, stl: Path) -> None:
    gmsh.initialize(interruptible=False)
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.occ.importShapes(str(step))
        gmsh.model.occ.synchronize()
        gmsh.model.mesh.generate(2)
        gmsh.write(str(stl))
    finally:
        gmsh.finalize()
