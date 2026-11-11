"""Built-in material library, and lookup of materials defined in a project.

Isotropic entries are textbook metals used by benchmarks and non-printed parts. Orthotropic
entries are printed materials curated with citations (ADR-010 action item 3, #61). Their ids
name the source, because measured values differ widely between printers and settings. A project
can still define its own (ADR-009 D3), and ``get_material`` finds those first.

Axes: 1 along the raster, 2 across the road in the layer, 3 the build direction. The full
evidence for every number is in ``docs/materials/pla-petg.md``.
"""

from collections.abc import Iterable

from tessera.model import (
    InputError,
    IsotropicDefinition,
    IsotropicMaterial,
    Material,
    MaterialDefinition,
    OrthotropicMaterial,
    OrthotropicStrength,
)

_ISOTROPIC: dict[str, IsotropicMaterial] = {
    # Typical handbook values for structural steel and 6061-T6 aluminium (densities approved
    # by the maintainer in iteration 5, when gravity loads first used them).
    "steel-generic": IsotropicMaterial(
        name="steel-generic", youngs_modulus=210_000.0, poisson_ratio=0.30, density=7.85e-9
    ),
    "aluminium-6061-t6": IsotropicMaterial(
        name="aluminium-6061-t6", youngs_modulus=68_900.0, poisson_ratio=0.33, density=2.70e-9
    ),
}


_ORTHOTROPIC: dict[str, OrthotropicMaterial] = {
    # PLA, desktop FFF. Gonabadi, Yadav, Bull, Int J Adv Manuf Technol 111 (2020) 695-709,
    # doi 10.1007/s00170-020-06138-4: Ender-3, 0.15 mm layers, 220 C, 100% infill; ASTM D638
    # and D5379 (Iosipescu) with DIC. Gaps filled from Song et al., Materials & Design 123
    # (2017) 154-164, doi 10.1016/j.matdes.2017.03.051 (0.2 mm layers, 220 C): e2, nu23, t2,
    # c1, c2. Assumed: the interlayer shear test does not say whether its plane is 13 or 23,
    # so g23 and s23 take it and s13 takes it too (lower than the 23 MPa across-road value).
    # Strengths are yield values where the source gives one; t1 and t3 are peak (t3 brittle).
    "pla-fff-gonabadi2020": OrthotropicMaterial(
        name="pla-fff-gonabadi2020",
        e1=3500.0, e2=4040.0, e3=2000.0,  # e2 Song Table 6
        nu12=0.32, nu13=0.35, nu23=0.37,  # nu23 Song Table 6
        g12=1270.0, g13=1210.0, g23=950.0,  # Gonabadi Table 7; g23 assumed (see above)
        density=1.22e-9,  # filament, measured (Gonabadi p698)
        strength=OrthotropicStrength(
            t1=55.0, t2=45.16, t3=5.0,  # t2 Song Table 4 (0.2% yield)
            s12=30.0, s13=17.5, s23=17.5,  # Gonabadi Table 7 yield; s13, s23 assumed
            c1=95.79, c2=94.87,  # Song Table 4 (0.2% yield); c3 not published
        ),
    ),
    # PETG, pellet-fed large-format extrusion, NOT desktop filament. Romeijn, Behrens, Paul,
    # Wei, Additive Manufacturing (2022) 103145, doi 10.1016/j.addma.2022.103145, Tables 5-8:
    # 5 mm roads, 2 mm layers, 220 C, heat-treated 10 min at 100 C before ASTM D3039 tests with
    # DIC. G from 45-degree off-axis tension. nu13 derived from the measured nu31 (0.461) as
    # nu31 * E1 / E3 (Table 6 labels it nu13; Eq 6 defines it as nu31). No strength: no source
    # gives interlayer tension or shear for PETG, so a PETG part gets the orthotropic note
    # instead of a strength ratio (decision on #61).
    "petg-pellet-romeijn2022": OrthotropicMaterial(
        name="petg-pellet-romeijn2022",
        e1=1898.7, e2=1272.4, e3=1894.8,
        nu12=0.401, nu13=0.462, nu23=0.373,
        g12=748.8, g13=611.7, g23=685.5,
        density=1.2742e-9,
    ),
}  # fmt: skip


def get_isotropic(material_id: str) -> IsotropicMaterial:
    try:
        return _ISOTROPIC[material_id]
    except KeyError:
        known = ", ".join(sorted(_ISOTROPIC))
        raise InputError(f"unknown material {material_id!r}; known: {known}") from None


def material_ids() -> list[str]:
    """Ids of the built-in materials."""
    return sorted([*_ISOTROPIC, *_ORTHOTROPIC])


def _from_definition(definition: MaterialDefinition) -> Material:
    """The material a project definition describes, with its strengths (ADR-010 D3)."""
    d = definition
    try:
        if isinstance(d, IsotropicDefinition):
            return IsotropicMaterial(
                name=d.id,
                youngs_modulus=d.youngs_modulus,
                poisson_ratio=d.poisson_ratio,
                density=d.density,
                strength=None if d.strength is None else d.strength.to_strength(),
            )
        return OrthotropicMaterial(
            name=d.id, e1=d.e1, e2=d.e2, e3=d.e3, nu12=d.nu12, nu13=d.nu13, nu23=d.nu23,
            g12=d.g12, g13=d.g13, g23=d.g23, density=d.density,
            strength=None if d.strength is None else d.strength.to_strength(),
        )  # fmt: skip
    except ValueError as exc:  # a definition changed after pydantic validated it
        raise InputError(f"material {d.id!r}: {exc}") from None


def get_material(material_id: str, defined: Iterable[MaterialDefinition] = ()) -> Material:
    """The material ``material_id``: a project definition with that id, else the library.

    Project definitions win, so a file keeps its meaning if a later library adds the same id
    (ADR-009 D3). Raises ``InputError`` for an unknown id or an invalid definition.
    """
    defined = list(defined)
    for definition in defined:
        if definition.id == material_id:
            return _from_definition(definition)
    if material_id in _ISOTROPIC:
        return _ISOTROPIC[material_id]
    if material_id in _ORTHOTROPIC:
        return _ORTHOTROPIC[material_id]
    known = ", ".join([*(d.id for d in defined), *material_ids()])
    raise InputError(f"unknown material {material_id!r}; known: {known}")
