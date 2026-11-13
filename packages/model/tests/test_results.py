"""Result types for assemblies (ADR-011 PR-1)."""

from tessera.model.results import (
    FailureMode,
    FlagKind,
    StrengthSpot,
    StressFlag,
    StressSpot,
)


def test_bond_edge_flag_names_its_bond() -> None:
    flag = StressFlag(FlagKind.BOND_EDGE, bond="glue")
    assert str(flag) == "near the edge of bond 'glue'"
    assert flag.support is None
    assert FlagKind.BOND_EDGE.value == "bond-edge"


def test_existing_flags_read_as_before() -> None:
    assert str(StressFlag(FlagKind.SUPPORT, support="base")) == "near support 'base'"
    assert str(StressFlag(FlagKind.SHARP_INSIDE_EDGE)) == "near a sharp inside edge"
    assert StressFlag(FlagKind.SHARP_INSIDE_EDGE).bond is None


def test_spots_name_their_part() -> None:
    stress = StressSpot(value=1.0, node=0, location=(0.0, 0.0, 0.0), face=3, part="cap")
    assert stress.part == "cap"
    mode = next(iter(FailureMode))
    strength = StrengthSpot(
        ratio=2.0, node=0, element=0, location=(0.0, 0.0, 0.0), face=None, mode=mode, part="cap"
    )
    assert strength.part == "cap"


def test_spot_part_is_optional_until_post_fills_it() -> None:
    assert StressSpot(value=1.0, node=0, location=(0.0, 0.0, 0.0), face=3).part is None
