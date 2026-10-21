"""Line supports through tessera.core: edges named by face pairs (ADR-006)."""

from pathlib import Path

import gmsh
import numpy as np
import pytest

import tessera.physics as physics
from tessera.core import describe_edges, describe_faces, edit, load_project, new_project, solve
from tessera.core.project_io import save_project
from tessera.mesh import mesh_step
from tessera.model import AnalysisKind, InputError, Project, SolveJob, SolveResult

L, W, H = 60.0, 10.0, 10.0


@pytest.fixture(scope="module")
def box_step(tmp_path_factory: pytest.TempPathFactory) -> Path:
    path = tmp_path_factory.mktemp("geom") / "beam.step"
    gmsh.initialize(interruptible=False)
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.occ.addBox(0, 0, 0, L, W, H)
        gmsh.model.occ.synchronize()
        gmsh.write(str(path))
    finally:
        gmsh.finalize()
    return path


@pytest.fixture
def project(box_step: Path) -> Project:
    p = new_project("beam")
    edit.add_part(p, box_step, material_id="steel-generic")
    return p


def _face(project: Project, normal: tuple[float, float, float]) -> int:
    return next(
        f.index for f in describe_faces(project, "beam") if f.normal == pytest.approx(normal)
    )


def _ends(project: Project) -> tuple[tuple[int, int], tuple[int, int], int]:
    """The bottom edges at both ends of the beam, and its top face."""
    bottom = _face(project, (0, 0, -1))
    left, right = _face(project, (-1, 0, 0)), _face(project, (1, 0, 0))
    return (left, bottom), (right, bottom), _face(project, (0, 0, 1))


def test_describe_edges_lists_the_box_edges(project: Project) -> None:
    edges = describe_edges(project, "beam")
    assert len(edges) == 12
    with pytest.raises(InputError, match="no part 'ghost'"):
        describe_edges(project, "ghost")


def test_add_edge_support_stores_sorted_pairs(project: Project) -> None:
    left, right, _ = _ends(project)
    pinned = edit.add_edge_support(project, "beam", [(left[1], left[0])])
    roller = edit.add_edge_support(project, "beam", [right], dofs=(3,), support_id="roller")
    assert (pinned.id, pinned.faces, pinned.dofs) == ("support-1", None, (1, 2, 3))
    assert pinned.edges is not None and pinned.edges.edges == (tuple(sorted(left)),)
    assert (roller.id, roller.part_id, roller.dofs) == ("roller", "beam", (3,))


@pytest.mark.parametrize(
    ("edges", "message"),
    [
        ([], "at least one edge"),
        ([(1, 1)], "two different faces"),
        ([(1, 99)], "no face 99"),
        ([(2, 1), (1, 2)], "listed twice"),
    ],
)
def test_bad_edges_are_refused_and_change_nothing(
    project: Project, edges: list[tuple[int, int]], message: str
) -> None:
    with pytest.raises(InputError, match=message):
        edit.add_edge_support(project, "beam", edges)
    assert project.study.supports == []


def test_faces_that_do_not_meet_are_refused_with_a_pointer(project: Project) -> None:
    left, right, top = _ends(project)
    with pytest.raises(InputError, match=r"do not meet at an edge.*tessera edges"):
        edit.add_edge_support(project, "beam", [left, (left[0], right[0])])
    with pytest.raises(InputError, match="no part 'ghost'"):
        edit.add_edge_support(project, "ghost", [left])
    with pytest.raises(InputError):
        edit.add_edge_support(project, "beam", [left], dofs=(4,))
    assert project.study.supports == []
    edit.add_edge_support(project, "beam", [(left[0], top)])  # left end meets the top


def test_a_part_held_on_a_line_cannot_be_removed(project: Project) -> None:
    left, _, _ = _ends(project)
    edit.add_edge_support(project, "beam", [left])
    with pytest.raises(InputError, match="support-1"):
        edit.remove(project, "beam")
    assert edit.remove(project, "support-1") == "support"


def test_line_supports_round_trip_through_a_project_file(project: Project, tmp_path: Path) -> None:
    import json

    left, right, _ = _ends(project)
    edit.add_edge_support(project, "beam", [left, right], dofs=(3,))
    edit.add_support(project, "beam", [left[0]])
    path = tmp_path / "beam.tessera.json"
    save_project(project, path)
    stored = json.loads(path.read_text(encoding="utf-8"))
    assert stored["schema_version"] == 5
    line, face = stored["study"]["supports"]
    assert line["edges"]["edges"] == sorted([sorted(left), sorted(right)])
    assert "faces" not in line and "edges" not in face
    assert load_project(path).study.supports == project.study.supports


def test_edge_nodes_on_a_box_mesh_lie_on_the_line(box_step: Path) -> None:
    mesh = mesh_step(box_step, part_id="beam", size=10.0)
    p = new_project("beam")
    edit.add_part(p, box_step)
    (left, bottom), _, _ = _ends(p)
    nodes = physics.edge_nodes(mesh, "beam", left, bottom)
    x = mesh.nodes[nodes]
    assert np.allclose(x[:, [0, 2]], 0.0), "on the line x = 0, z = 0"
    ys = np.sort(x[:, 1])
    assert ys[0] == pytest.approx(0.0) and ys[-1] == pytest.approx(W)
    # Corners and midsides alternate along the line: an odd count, evenly interleaved.
    assert len(ys) % 2 == 1 and len(ys) >= 3
    corners = mesh.nodes[np.unique(mesh.surfaces[f"face:{left}"][:, :3])]
    on_line = corners[np.isclose(corners[:, 2], 0.0)]
    assert len(ys) == 2 * len(on_line) - 1


class _FakeSolver:
    name = "fake"

    def __init__(self) -> None:
        self.job: SolveJob | None = None

    def capabilities(self) -> frozenset[AnalysisKind]:
        return frozenset({AnalysisKind.STATIC})

    def run(self, job: SolveJob, workdir: Path) -> SolveResult:
        self.job = job
        n = job.mesh.n_nodes
        return SolveResult(
            displacements=np.zeros((n, 3)), stresses=np.zeros((n, 6)), solver=self.name
        )


def test_solve_marks_line_supports_in_the_job(project: Project, tmp_path: Path) -> None:
    left, _, top = _ends(project)
    edit.add_edge_support(project, "beam", [left], support_id="pin")
    edit.add_support(project, "beam", [_face(project, (1, 0, 0))], support_id="clamp")
    edit.add_pressure_load(project, "beam", [top], 0.1)
    edit.set_mesh_size(project, 10.0)
    fake = _FakeSolver()
    solve(project, workdir=tmp_path, adapter=fake)
    assert fake.job is not None
    assert [(s.node_set, s.line) for s in fake.job.supports] == [
        ("support:pin", True),
        ("support:clamp", False),
    ]
    held = fake.job.mesh.nodes[fake.job.mesh.node_sets["support:pin"]]
    assert np.allclose(held[:, [0, 2]], 0.0)


@pytest.mark.solver
def test_a_simply_supported_box_sags_under_a_top_pressure(project: Project, tmp_path: Path) -> None:
    left, right, top = _ends(project)
    edit.add_edge_support(project, "beam", [left], support_id="pin")
    edit.add_edge_support(project, "beam", [right], dofs=(3,), support_id="roller")
    edit.add_pressure_load(project, "beam", [top], 1.0)
    edit.set_mesh_size(project, 5.0)
    result = solve(project, workdir=tmp_path)
    x, u = result.mesh.nodes, result.solution.displacements
    middle = np.isclose(x[:, 0], L / 2) & np.isclose(x[:, 2], 0.0)
    assert middle.any()
    assert (u[middle, 2] < 0).all(), "the middle of the span moves down"
    for node_set in ("support:pin", "support:roller"):
        held = result.mesh.node_sets[node_set]
        assert np.abs(u[held, 2]).max() < 1e-9
        assert np.abs(u[middle, 2]).min() > 100 * np.abs(u[held, 2]).max()
    pinned = result.mesh.node_sets["support:pin"]
    assert np.abs(u[pinned]).max() < 1e-9
    roller = result.mesh.node_sets["support:roller"]
    assert np.abs(u[roller, 0]).max() > 0, "the roller end slides along the beam"


def test_an_stl_part_gets_the_same_line_support(box_step: Path, tmp_path: Path) -> None:
    stl = tmp_path / "beam.stl"
    gmsh.initialize(interruptible=False)
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.occ.importShapes(str(box_step))
        gmsh.model.occ.synchronize()
        gmsh.model.mesh.generate(2)
        gmsh.write(str(stl))
    finally:
        gmsh.finalize()
    p = new_project("beam")
    edit.add_part(p, stl, material_id="steel-generic")
    left, _, top = _ends(p)
    edit.add_edge_support(p, "beam", [left], support_id="pin")
    edit.add_pressure_load(p, "beam", [top], 0.1)
    edit.set_mesh_size(p, 10.0)
    fake = _FakeSolver()
    solve(p, workdir=tmp_path / "run", adapter=fake)
    assert fake.job is not None and fake.job.supports[0].line
    held = fake.job.mesh.nodes[fake.job.mesh.node_sets["support:pin"]]
    assert np.allclose(held[:, [0, 2]], 0.0)
    assert sorted(held[:, 1])[0] == pytest.approx(0.0) and max(held[:, 1]) == pytest.approx(W)
