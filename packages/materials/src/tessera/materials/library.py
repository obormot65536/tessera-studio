"""Built-in material library, and lookup of materials defined in a project.

v0.1 holds only textbook isotropic metals used by benchmarks and non-printed parts. Printed
materials (orthotropic, with citations) arrive in v0.5 through the materials-curator workflow.
Until then a project defines its own (ADR-009 D3), and ``get_material`` finds those first.
"""

from collections.abc import Iterable

from tessera.model import (
    InputError,
    IsotropicDefinition,
    IsotropicMaterial,
    Material,
    MaterialDefinition,
    OrthotropicMaterial,
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


def get_isotropic(material_id: str) -> IsotropicMaterial:
    try:
        return _ISOTROPIC[material_id]
    except KeyError:
        known = ", ".join(sorted(_ISOTROPIC))
        raise InputError(f"unknown material {material_id!r}; known: {known}") from None


def material_ids() -> list[str]:
    """Ids of the built-in materials."""
    return sorted(_ISOTROPIC)


def _from_definition(definition: MaterialDefinition) -> Material:
    if isinstance(definition, IsotropicDefinition):
        return IsotropicMaterial(
            name=definition.id,
            youngs_modulus=definition.youngs_modulus,
            poisson_ratio=definition.poisson_ratio,
            density=definition.density,
        )
    d = definition
    try:
        return OrthotropicMaterial(
            name=d.id, e1=d.e1, e2=d.e2, e3=d.e3, nu12=d.nu12, nu13=d.nu13, nu23=d.nu23,
            g12=d.g12, g13=d.g13, g23=d.g23, density=d.density,
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
    known = ", ".join([*(d.id for d in defined), *sorted(_ISOTROPIC)])
    raise InputError(f"unknown material {material_id!r}; known: {known}")
