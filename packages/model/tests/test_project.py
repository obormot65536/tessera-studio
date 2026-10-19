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
