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


# ADR-010 action item 3 (#61): the first printed materials in the library. These pin the
# numbers the maintainer approved in the PR; changing one is a rule-6 change.
def test_printed_materials_are_listed() -> None:
    ids = material_ids()
    assert "pla-fff-gonabadi2020" in ids
    assert "petg-pellet-romeijn2022" in ids


def test_pla_entry_has_the_approved_constants_and_strengths() -> None:
    from tessera.model import OrthotropicStrength
    from tessera.model.elasticity import check_orthotropic

    pla = get_material("pla-fff-gonabadi2020")
    assert pla == OrthotropicMaterial(
        name="pla-fff-gonabadi2020",
        e1=3500.0, e2=4040.0, e3=2000.0,
        nu12=0.32, nu13=0.35, nu23=0.37,
        g12=1270.0, g13=1210.0, g23=950.0,
        density=1.22e-9,
        strength=OrthotropicStrength(
            t1=55.0, t2=45.16, t3=5.0, s12=30.0, s13=17.5, s23=17.5, c1=95.79, c2=94.87
        ),
    )  # fmt: skip
    assert isinstance(pla, OrthotropicMaterial)
    check_orthotropic(
        pla.e1, pla.e2, pla.e3, pla.nu12, pla.nu13, pla.nu23, pla.g12, pla.g13, pla.g23
    )


def test_petg_entry_is_elastic_only() -> None:
    from tessera.model.elasticity import check_orthotropic

    petg = get_material("petg-pellet-romeijn2022")
    assert petg == OrthotropicMaterial(
        name="petg-pellet-romeijn2022",
        e1=1898.7, e2=1272.4, e3=1894.8,
        nu12=0.401, nu13=0.462, nu23=0.373,
        g12=748.8, g13=611.7, g23=685.5,
        density=1.2742e-9,
    )  # fmt: skip
    assert isinstance(petg, OrthotropicMaterial)
    assert petg.strength is None  # no interlayer strengths published (#61): not checked
    check_orthotropic(
        petg.e1, petg.e2, petg.e3, petg.nu12, petg.nu13, petg.nu23, petg.g12, petg.g13,
        petg.g23,
    )  # fmt: skip


def test_a_project_definition_wins_over_a_printed_library_material() -> None:
    mine = OrthotropicDefinition.model_validate({"id": "pla-fff-gonabadi2020", **_TO})
    assert get_material("pla-fff-gonabadi2020", [mine]) == OrthotropicMaterial(
        name="pla-fff-gonabadi2020", **_TO
    )


def test_printed_materials_are_not_isotropic() -> None:
    with pytest.raises(InputError, match="unknown material"):
        get_isotropic("pla-fff-gonabadi2020")
