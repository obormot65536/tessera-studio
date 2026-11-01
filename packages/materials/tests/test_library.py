from typing import Any

import pytest

from tessera.materials import get_isotropic, get_material, material_ids
from tessera.model import (
    InputError,
    IsotropicDefinition,
    IsotropicMaterial,
    OrthotropicDefinition,
    OrthotropicMaterial,
)


def test_steel() -> None:
    steel = get_isotropic("steel-generic")
    assert steel.youngs_modulus == 210_000.0
    assert steel.poisson_ratio == 0.30


def test_unknown_material_lists_known_ids() -> None:
    with pytest.raises(InputError, match="steel-generic"):
        get_isotropic("unobtainium")
    assert "aluminium-6061-t6" in material_ids()


# ADR-009 D1/D3: get_material resolves project definitions first, then the library.
_TO: dict[str, Any] = dict(
    e1=3000.0,
    e2=2400.0,
    e3=1800.0,
    nu12=0.35,
    nu13=0.30,
    nu23=0.30,
    g12=1000.0,
    g13=700.0,
    g23=650.0,
)


def test_get_material_finds_library_materials() -> None:
    assert get_material("steel-generic") == get_isotropic("steel-generic")


def test_project_definitions_win_over_the_library() -> None:
    mine = IsotropicDefinition(id="steel-generic", youngs_modulus=200_000.0, poisson_ratio=0.29)
    pla = OrthotropicDefinition.model_validate({"id": "pla-xyz", "density": 1.24e-9, **_TO})
    defined = [mine, pla]
    steel = get_material("steel-generic", defined)
    assert steel == IsotropicMaterial(
        name="steel-generic", youngs_modulus=200_000.0, poisson_ratio=0.29
    )
    assert get_material("pla-xyz", defined) == OrthotropicMaterial(
        name="pla-xyz", density=1.24e-9, **_TO
    )


def test_get_material_unknown_id_is_an_input_error_listing_both_sources() -> None:
    defined = [OrthotropicDefinition.model_validate({"id": "pla-xyz", **_TO})]
    with pytest.raises(InputError, match=r"unknown material 'unobtainium'.*pla-xyz.*steel"):
        get_material("unobtainium", defined)


def test_an_invalid_definition_edited_after_creation_is_an_input_error() -> None:
    bad = OrthotropicDefinition.model_validate({"id": "pla-xyz", **_TO})
    bad.nu12 = 1.2  # no validation on assignment
    with pytest.raises(InputError, match="pla-xyz"):
        get_material("pla-xyz", [bad])


# ADR-010 D3: strengths a project defines reach the material the solve and post see.
def test_project_strengths_reach_the_material() -> None:
    from tessera.model import IsotropicStrength, OrthotropicStrength

    pla = OrthotropicDefinition.model_validate(
        {"id": "pla-xyz", **_TO, "strength": {"t1": 40.0, "t2": 30.0, "t3": 15.0,
         "s12": 20.0, "s13": 12.0, "s23": 12.0, "c3": 40.0}}
    )  # fmt: skip
    resin = IsotropicDefinition.model_validate(
        {"id": "resin", "youngs_modulus": 2500.0, "poisson_ratio": 0.35,
         "strength": {"tensile": 50.0, "brittle": True}}
    )  # fmt: skip
    plain = IsotropicDefinition(id="plain", youngs_modulus=2500.0, poisson_ratio=0.35)
    defined = [pla, resin, plain]
    assert get_material("pla-xyz", defined).strength == OrthotropicStrength(
        t1=40.0, t2=30.0, t3=15.0, s12=20.0, s13=12.0, s23=12.0, c3=40.0
    )
    assert get_material("resin", defined).strength == IsotropicStrength(50.0, brittle=True)
    assert get_material("plain", defined).strength is None
