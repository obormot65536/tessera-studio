"""Strength data and strength-result types (ADR-010 PR-1)."""

import dataclasses
import math
from typing import Any

import numpy as np
import pytest
from pydantic import ValidationError

import tessera.model as model
from tessera.model import (
    FailureCriterion,
    FailureMode,
    IsotropicMaterial,
    IsotropicStrength,
    MeshCheck,
    OrthotropicMaterial,
    OrthotropicStrength,
    Project,
    StrengthResult,
    StrengthSpot,
    StudyResult,
)

TO: dict[str, Any] = {
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
TENSILE = {"t1": 40.0, "t2": 30.0, "t3": 15.0, "s12": 20.0, "s13": 12.0, "s23": 12.0}


def test_new_names_are_exported() -> None:
    names = (
        "OrthotropicStrength",
        "IsotropicStrength",
        "Strength",
        "FailureCriterion",
        "FailureMode",
        "StrengthSpot",
        "StrengthResult",
        "IsotropicStrengthDefinition",
        "OrthotropicStrengthDefinition",
    )
    for name in names:
        assert name in model.__all__
        assert hasattr(model, name)


# --- Orthotropic -------------------------------------------------------------------------------


def test_orthotropic_strength_fills_missing_compressive_values_and_lists_them() -> None:
    s = OrthotropicStrength(**TENSILE, c2=45.0)
    assert s.compressive == (40.0, 45.0, 15.0)
    assert s.assumed == ("c1", "c3")
    full = OrthotropicStrength(**TENSILE, c1=50.0, c2=45.0, c3=40.0)
    assert full.compressive == (50.0, 45.0, 40.0)
    assert full.assumed == ()


@pytest.mark.parametrize("name", [*TENSILE, "c1", "c2", "c3"])
@pytest.mark.parametrize("bad", [0.0, -1.0, math.nan, math.inf])
def test_orthotropic_strength_must_be_positive_and_finite(name: str, bad: float) -> None:
    with pytest.raises(ValueError, match=f"{name} must be a positive finite strength"):
        OrthotropicStrength(**{**TENSILE, name: bad})


def test_orthotropic_strength_is_frozen() -> None:
    s = OrthotropicStrength(**TENSILE)
    with pytest.raises(dataclasses.FrozenInstanceError):
        s.t1 = 1.0  # type: ignore[misc]


# --- Isotropic ---------------------------------------------------------------------------------


def test_ductile_isotropic_strength_is_the_yield_strength_alone() -> None:
    s = IsotropicStrength(250.0)
    assert (s.tensile, s.compressive, s.brittle, s.assumed) == (250.0, None, False, ())
    with pytest.raises(ValueError, match="compressive strength is used only for brittle"):
        IsotropicStrength(250.0, compressive=300.0)


def test_brittle_isotropic_strength_fills_a_missing_compressive_value() -> None:
    s = IsotropicStrength(50.0, brittle=True)
    assert s.compressive_or_tensile == 50.0
    assert s.assumed == ("compressive",)
    given = IsotropicStrength(50.0, compressive=80.0, brittle=True)
    assert given.compressive_or_tensile == 80.0
    assert given.assumed == ()


@pytest.mark.parametrize("bad", [0.0, -5.0, math.nan, math.inf])
def test_isotropic_strength_must_be_positive_and_finite(bad: float) -> None:
    with pytest.raises(ValueError, match="tensile must be a positive finite strength"):
        IsotropicStrength(bad)
    with pytest.raises(ValueError, match="compressive must be a positive finite strength"):
        IsotropicStrength(50.0, compressive=bad, brittle=True)


# --- On the materials --------------------------------------------------------------------------


def test_materials_have_no_strength_by_default() -> None:
    assert IsotropicMaterial("steel", 210000.0, 0.3).strength is None
    assert OrthotropicMaterial("to", **TO).strength is None


def test_materials_carry_a_strength_of_their_kind() -> None:
    iso = IsotropicMaterial("al", 68900.0, 0.33, strength=IsotropicStrength(276.0))
    ortho = OrthotropicMaterial("to", **TO, strength=OrthotropicStrength(**TENSILE))
    assert iso.strength == IsotropicStrength(276.0)
    assert ortho.strength == OrthotropicStrength(**TENSILE)


def test_a_strength_of_the_wrong_kind_is_rejected() -> None:
    with pytest.raises(ValueError, match="material 'al' is isotropic and needs an Isotropic"):
        IsotropicMaterial("al", 68900.0, 0.33, strength=OrthotropicStrength(**TENSILE))  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="material 'to' is orthotropic and needs an Orthotropic"):
        OrthotropicMaterial("to", **TO, strength=IsotropicStrength(50.0))  # type: ignore[arg-type]


# --- Project definitions (schema 7) ------------------------------------------------------------


def _ortho_definition(**strength: object) -> dict[str, object]:
    return {"kind": "orthotropic", "id": "my-pla", **TO, "strength": strength}


def test_project_materials_round_trip_with_strengths() -> None:
    data = {
        "name": "bracket",
        "materials": [
            _ortho_definition(**TENSILE, c3=40.0),
            {
                "kind": "isotropic",
                "id": "resin",
                "youngs_modulus": 2500.0,
                "poisson_ratio": 0.35,
                "strength": {"tensile": 50.0, "brittle": True},
            },
        ],
    }
    project = Project.model_validate(data)
    ortho, iso = project.materials
    assert isinstance(ortho.strength, model.OrthotropicStrengthDefinition)
    assert isinstance(iso.strength, model.IsotropicStrengthDefinition)
    assert ortho.strength.t3 == 15.0 and ortho.strength.c3 == 40.0 and ortho.strength.c1 is None
    assert iso.strength.brittle and iso.strength.compressive is None
    assert Project.model_validate(project.model_dump(mode="json")) == project


def test_project_material_strength_defaults_to_none() -> None:
    project = Project.model_validate(
        {
            "name": "x",
            "materials": [
                {"kind": "isotropic", "id": "a", "youngs_modulus": 1.0, "poisson_ratio": 0.3}
            ],
        }
    )
    assert project.materials[0].strength is None


@pytest.mark.parametrize(
    ("strength", "message"),
    [
        ({**TENSILE, "t3": 0.0}, "greater than 0"),
        ({**TENSILE, "c1": -1.0}, "greater than 0"),
        ({key: value for key, value in TENSILE.items() if key != "s13"}, "s13"),
        ({**TENSILE, "x": 1.0}, "Extra inputs are not permitted"),
    ],
)
def test_invalid_orthotropic_strength_is_rejected(strength: dict[str, float], message: str) -> None:
    with pytest.raises(ValidationError, match=message):
        Project.model_validate({"name": "x", "materials": [_ortho_definition(**strength)]})


def test_ductile_project_material_with_a_compressive_strength_is_rejected() -> None:
    material = {
        "kind": "isotropic",
        "id": "al",
        "youngs_modulus": 68900.0,
        "poisson_ratio": 0.33,
        "strength": {"tensile": 276.0, "compressive": 276.0},
    }
    with pytest.raises(ValidationError, match="material 'al': compressive strength is used only"):
        Project.model_validate({"name": "x", "materials": [material]})


def test_strength_definitions_build_the_solver_strengths() -> None:
    project = Project.model_validate(
        {"name": "x", "materials": [_ortho_definition(**TENSILE, c1=50.0)]}
    )
    definition = project.materials[0].strength
    assert definition is not None
    assert definition.to_strength() == OrthotropicStrength(**TENSILE, c1=50.0)
    iso = model.IsotropicStrengthDefinition(tensile=50.0, compressive=80.0, brittle=True)
    assert iso.to_strength() == IsotropicStrength(50.0, compressive=80.0, brittle=True)


def test_schema_6_project_loads_unchanged() -> None:
    from tessera.model import SCHEMA_VERSION

    material = {"kind": "orthotropic", "id": "my-pla", **TO, "density": None}
    data = {"schema_version": 6, "name": "bracket", "materials": [material]}
    project = Project.model_validate(data)
    assert SCHEMA_VERSION >= 7
    assert project.materials[0].strength is None
    dumped = project.model_dump(mode="json")
    assert dumped["materials"] == [{**material, "strength": None}]


# --- Results -----------------------------------------------------------------------------------


def _spot(ratio: float, **kw: object) -> StrengthSpot:
    return StrengthSpot(
        ratio=ratio,
        node=3,
        element=1,
        location=(1.0, 2.0, 3.0),
        face=4,
        mode=FailureMode.ACROSS_LAYERS_TENSION,
        **kw,  # type: ignore[arg-type]
    )


def test_strength_spot_knows_whether_it_is_flagged() -> None:
    from tessera.model import FlagKind, StressFlag

    assert not _spot(2.0).flagged
    assert _spot(2.0, flags=(StressFlag(FlagKind.SHARP_INSIDE_EDGE),)).flagged


def test_failure_mode_and_criterion_values() -> None:
    assert FailureCriterion.TSAI_WU == "tsai-wu"
    assert FailureCriterion.VON_MISES == "von-mises"
    assert FailureCriterion.MAX_PRINCIPAL == "max-principal"
    assert {m.value for m in FailureMode} == {
        "along-raster-tension",
        "along-raster-compression",
        "across-roads-tension",
        "across-roads-compression",
        "across-layers-tension",
        "across-layers-compression",
        "in-layer-shear",
        "interlayer-shear",
        "von-mises",
        "max-principal-tension",
        "max-principal-compression",
    }


def test_strength_result_defaults() -> None:
    result = StrengthResult(
        criteria={"to": FailureCriterion.TSAI_WU},
        lowest=_spot(1.5),
        clear_lowest=None,
        face_lowest=(_spot(1.5),),
    )
    assert result.assumed == {}
    assert result.unchecked == ()


def test_study_result_and_mesh_check_fields_default_to_none() -> None:
    fields = {f.name: f for f in dataclasses.fields(StudyResult)}
    assert fields["strength"].default is None and fields["strength"].kw_only
    check = {f.name: f for f in dataclasses.fields(MeshCheck)}
    assert check["lowest_strength_ratio"].default is None
    assert check["lowest_strength_ratio"].kw_only
    assert np.isinf(_spot(math.inf).ratio)
