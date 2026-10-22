from pathlib import Path

import pytest
from pydantic import ValidationError

from tessera.model import Assembly, Connection, ConnectionKind, FaceSet, PartInstance, Project


def _part(pid: str) -> PartInstance:
    return PartInstance(id=pid, name=pid, geometry=Path(f"{pid}.step"))


def test_single_part_is_an_assembly() -> None:
    project = Project(name="bracket", assembly=Assembly(parts=[_part("p1")]))
    assert len(project.assembly.parts) == 1


def test_json_round_trip() -> None:
    assembly = Assembly(
        parts=[_part("p1"), _part("p2")],
        connections=[
            Connection(
                id="c1",
                kind=ConnectionKind.BOLT,
                a=FaceSet(part_id="p1", faces=(3,)),
                b=FaceSet(part_id="p2", faces=(7,)),
                preload_n=500.0,
            )
        ],
    )
    project = Project(name="hinge", assembly=assembly)
    assert Project.model_validate_json(project.model_dump_json()) == project


def test_connection_to_unknown_part_is_rejected() -> None:
    with pytest.raises(ValidationError, match="unknown part"):
        Assembly(
            parts=[_part("p1")],
            connections=[
                Connection(
                    id="c1",
                    kind=ConnectionKind.BONDED,
                    a=FaceSet(part_id="p1", faces=(0,)),
                    b=FaceSet(part_id="ghost", faces=(0,)),
                )
            ],
        )


def test_duplicate_part_ids_are_rejected() -> None:
    with pytest.raises(ValidationError, match="unique"):
        Assembly(parts=[_part("p1"), _part("p1")])


def test_face_set_needs_at_least_one_face() -> None:
    with pytest.raises(ValidationError, match="at least 1"):
        FaceSet(part_id="p1", faces=())


def test_support_fixes_all_directions_by_default_or_selected_ones() -> None:
    from tessera.model import Support

    faces = FaceSet(part_id="p1", faces=(1,))
    assert Support(id="s", faces=faces).dofs == (1, 2, 3)
    roller = Support(id="sym", faces=faces, dofs=(1,))
    assert Support.model_validate_json(roller.model_dump_json()) == roller


@pytest.mark.parametrize("dofs", [(), (1, 1), (4,), (0,)])
def test_support_rejects_invalid_dofs(dofs: tuple[int, ...]) -> None:
    from tessera.model import Support

    with pytest.raises(ValidationError):
        Support.model_validate({"id": "s", "faces": {"part_id": "p1", "faces": [1]}, "dofs": dofs})


def test_triangle_surface_validates_indices_and_face_ids() -> None:
    import numpy as np

    from tessera.model import TriangleSurface

    vertices = np.zeros((4, 3))
    tris = np.array([[0, 1, 2], [0, 2, 3]])
    TriangleSurface(vertices, tris, np.array([1, 2]))  # valid
    with pytest.raises(ValueError, match="vertex index"):
        TriangleSurface(vertices, np.array([[0, 1, 9], [0, 2, 3]]), np.array([1, 1]))
    with pytest.raises(ValueError, match=r"1\.\.n"):
        TriangleSurface(vertices, tris, np.array([1, 3]))
    surface = TriangleSurface(vertices, tris, np.array([1, 1]))
    assert surface == surface  # identity comparison, no ambiguous array truth value


def test_study_holds_force_and_pressure_loads_and_gravity() -> None:
    from tessera.model import ForceLoad, PressureLoad, Study

    faces = FaceSet(part_id="p1", faces=(1,))
    study = Study(
        loads=[
            ForceLoad(id="f", faces=faces, force=(0.0, 0.0, -10.0)),
            PressureLoad(id="p", faces=faces, pressure=2.5),
        ],
        gravity=(0.0, 0.0, -9810.0),
    )
    again = Study.model_validate_json(study.model_dump_json())
    assert again == study
    assert [type(load).__name__ for load in again.loads] == ["ForceLoad", "PressureLoad"]


def test_schema_2_force_loads_without_kind_still_load() -> None:
    from tessera.model import ForceLoad, Study

    study = Study.model_validate(
        {"loads": [{"id": "f", "faces": {"part_id": "p1", "faces": [1]}, "force": [1, 0, 0]}]}
    )
    assert isinstance(study.loads[0], ForceLoad)


def test_unknown_load_kind_is_rejected() -> None:
    from tessera.model import Study

    with pytest.raises(ValidationError):
        Study.model_validate(
            {"loads": [{"kind": "torque", "id": "m", "faces": {"part_id": "p1", "faces": [1]}}]}
        )


def test_moment_load_round_trips() -> None:
    from tessera.model import MomentLoad, Study

    load = MomentLoad(id="t", faces=FaceSet(part_id="p1", faces=(2,)), moment=(0, 0, 5e4))
    again = Study.model_validate_json(Study(loads=[load]).model_dump_json())
    assert again.loads == [load]


@pytest.mark.parametrize("bad", [float("nan"), float("inf"), -float("inf")])
def test_numbers_must_be_finite(bad: float) -> None:
    from tessera.model import ForceLoad, PressureLoad, Study

    faces = FaceSet(part_id="p1", faces=(1,))
    with pytest.raises(ValidationError):
        ForceLoad(id="f", faces=faces, force=(0.0, bad, 0.0))
    with pytest.raises(ValidationError):
        PressureLoad(id="p", faces=faces, pressure=bad)
    with pytest.raises(ValidationError):
        Study(gravity=(0.0, 0.0, bad))
    with pytest.raises(ValidationError):
        Study(mesh_size=bad)


def test_edge_set_stores_pairs_ascending() -> None:
    from tessera.model import EdgeSet

    edges = EdgeSet(part_id="p1", edges=((5, 3), (5, 6)))
    assert edges.edges == ((3, 5), (5, 6))
    assert EdgeSet.model_validate_json(edges.model_dump_json()) == edges
    assert edges.model_dump(mode="json") == {"part_id": "p1", "edges": [[3, 5], [5, 6]]}


@pytest.mark.parametrize(
    ("edges", "message"),
    [
        ((), "at least 1"),
        (((3, 3),), "two different faces"),
        (((3, 5), (5, 3)), "listed twice"),
        (((3, 5, 6),), "at most 2"),
    ],
)
def test_edge_set_rejects_bad_pairs(edges: object, message: str) -> None:
    from tessera.model import EdgeSet

    with pytest.raises(ValidationError, match=message):
        EdgeSet.model_validate({"part_id": "p1", "edges": edges})


def test_support_holds_faces_or_edges() -> None:
    from tessera.model import EdgeSet, Support

    line = Support(id="s", edges=EdgeSet(part_id="p2", edges=((1, 2),)), dofs=(3,))
    assert line.part_id == "p2" and line.faces is None
    assert Support.model_validate_json(line.model_dump_json()) == line
    assert "faces" not in line.model_dump()
    face = Support(id="f", faces=FaceSet(part_id="p1", faces=(1,)))
    assert face.part_id == "p1"
    assert "edges" not in face.model_dump(mode="json"), "face supports are written as before"
    with pytest.raises(ValidationError, match="exactly one"):
        Support(id="s")
    with pytest.raises(ValidationError, match="exactly one"):
        Support(
            id="s",
            faces=FaceSet(part_id="p1", faces=(1,)),
            edges=EdgeSet(part_id="p1", edges=((1, 2),)),
        )


def test_line_support_on_an_unknown_part_is_rejected() -> None:
    from tessera.model import EdgeSet, Study, Support

    study = Study(supports=[Support(id="s", edges=EdgeSet(part_id="ghost", edges=((1, 2),)))])
    with pytest.raises(ValidationError, match="unknown part 'ghost'"):
        Project(name="p", assembly=Assembly(parts=[_part("p1")]), study=study)


def test_schema_4_project_loads_unchanged() -> None:
    from tessera.model import SCHEMA_VERSION

    supports = [{"id": "s", "faces": {"part_id": "p1", "faces": [1]}, "dofs": [3]}]
    data = {
        "schema_version": 4,
        "name": "bracket",
        "assembly": {"parts": [{"id": "p1", "name": "p1", "geometry": "p1.step"}]},
        "study": {"supports": supports},
    }
    project = Project.model_validate(data)
    support = project.study.supports[0]
    assert (support.part_id, support.edges, support.dofs) == ("p1", None, (3,))
    assert SCHEMA_VERSION >= 5
    assert project.model_dump(mode="json")["study"]["supports"] == supports


# Catalog TO constants (ADR-009 test plan); Delta = 0.7427, positive definite.
_TO: dict[str, object] = {
    "e1": 3000.0,
    "e2": 2400.0,
    "e3": 1800.0,
    "nu12": 0.35,
    "nu13": 0.30,
    "nu23": 0.30,
    "g12": 1000.0,
    "g13": 700.0,
    "g23": 650.0,
}


def _orthotropic(mid: str = "pla-xyz", **changes: object) -> dict[str, object]:
    return {"kind": "orthotropic", "id": mid, **_TO, "density": 1.24e-9, **changes}


def test_project_materials_round_trip_through_json() -> None:
    from tessera.model import IsotropicDefinition, OrthotropicDefinition

    raw = {
        "name": "p",
        "materials": [
            _orthotropic(),
            {"kind": "isotropic", "id": "petg", "youngs_modulus": 2100.0, "poisson_ratio": 0.38},
        ],
    }
    project = Project.model_validate(raw)
    ortho, iso = project.materials
    assert isinstance(ortho, OrthotropicDefinition) and isinstance(iso, IsotropicDefinition)
    assert (ortho.e2, ortho.nu23, ortho.g23, ortho.density) == (2400.0, 0.30, 650.0, 1.24e-9)
    assert iso.density is None
    again = Project.model_validate_json(project.model_dump_json())
    assert again == project
    assert again.model_dump(mode="json")["materials"][0] == _orthotropic()


def test_projects_have_no_materials_by_default() -> None:
    assert Project(name="p").materials == []


@pytest.mark.parametrize(
    ("changes", "message"),
    [
        ({"nu12": 1.2, "e2": 3000.0}, "positive definite"),
        ({"g13": 0.0}, "g13 must be positive"),
        ({"e3": -1.0}, "e3 must be positive"),
    ],
)
def test_invalid_orthotropic_definition_is_rejected(
    changes: dict[str, object], message: str
) -> None:
    with pytest.raises(ValidationError, match=message) as info:
        Project.model_validate({"name": "p", "materials": [_orthotropic("bad", **changes)]})
    assert "material 'bad'" in str(info.value), "the message names the material"


@pytest.mark.parametrize(
    "changes",
    [{"youngs_modulus": 0.0}, {"poisson_ratio": 0.5}, {"poisson_ratio": -1.0}, {"density": 0.0}],
)
def test_invalid_isotropic_definition_is_rejected(changes: dict[str, object]) -> None:
    from tessera.model import IsotropicDefinition

    data = {"id": "m", "youngs_modulus": 2000.0, "poisson_ratio": 0.3, **changes}
    with pytest.raises(ValidationError):
        IsotropicDefinition.model_validate(data)


def test_unknown_material_kind_is_rejected() -> None:
    with pytest.raises(ValidationError, match="anisotropic"):
        Project.model_validate({"name": "p", "materials": [{"kind": "anisotropic", "id": "m"}]})


def test_duplicate_material_ids_are_rejected() -> None:
    iso = {"kind": "isotropic", "id": "pla", "youngs_modulus": 3000.0, "poisson_ratio": 0.35}
    with pytest.raises(ValidationError, match="material id 'pla' is used twice"):
        Project.model_validate({"name": "p", "materials": [iso, _orthotropic("pla")]})


def test_raster_angle_round_trips_and_defaults_to_zero() -> None:
    from tessera.model import PrintSetup

    assert PrintSetup(process="fdm").raster_angle == 0.0
    setup = PrintSetup(process="fdm", build_direction=(1.0, 0.0, 0.0), raster_angle=45.0)
    part = PartInstance(id="p1", name="p1", geometry=Path("p1.step"), print_setup=setup)
    project = Project(name="p", assembly=Assembly(parts=[part]))
    again = Project.model_validate_json(project.model_dump_json())
    assert again.assembly.parts[0].print_setup == setup
    assert (
        again.model_dump(mode="json")["assembly"]["parts"][0]["print_setup"]["raster_angle"] == 45
    )


@pytest.mark.parametrize("bad", [float("nan"), float("inf"), -float("inf")])
def test_raster_angle_must_be_finite(bad: float) -> None:
    from tessera.model import PrintSetup

    with pytest.raises(ValidationError, match="finite number"):
        PrintSetup(process="fdm", raster_angle=bad)


def test_schema_5_project_loads_unchanged() -> None:
    from tessera.model import SCHEMA_VERSION

    setup = {"process": "fdm", "slicer_profile": None, "build_direction": [0.0, 1.0, 0.0]}
    part = {"id": "p1", "name": "p1", "geometry": "p1.step", "material_id": "steel"}
    supports = [{"id": "s", "edges": {"part_id": "p1", "edges": [[1, 2]]}, "dofs": [3]}]
    data = {
        "schema_version": 5,
        "name": "bracket",
        "assembly": {"parts": [{**part, "print_setup": setup}]},
        "study": {"supports": supports},
    }
    project = Project.model_validate(data)
    assert SCHEMA_VERSION == 6
    assert project.materials == []
    loaded = project.assembly.parts[0]
    assert loaded.material_id == "steel" and loaded.print_setup is not None
    assert loaded.print_setup.build_direction == (0.0, 1.0, 0.0)
    assert loaded.print_setup.raster_angle == 0.0
    dumped = project.model_dump(mode="json")
    assert dumped["assembly"]["parts"][0]["print_setup"] == {**setup, "raster_angle": 0.0}
    assert dumped["study"]["supports"] == supports
    assert dumped["materials"] == []
