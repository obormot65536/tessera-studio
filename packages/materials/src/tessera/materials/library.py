"""Built-in material library.

v0.1 holds only textbook isotropic metals used by benchmarks and non-printed parts. Printed
materials (orthotropic, with citations) arrive in v0.5 through the materials-curator workflow.
"""

from tessera.model import InputError, IsotropicMaterial

_ISOTROPIC: dict[str, IsotropicMaterial] = {
    # Typical handbook values for structural steel and 6061-T6 aluminium.
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
    return sorted(_ISOTROPIC)
