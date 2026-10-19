import pytest

from tessera.materials import get_isotropic, material_ids
from tessera.model import InputError


def test_steel() -> None:
    steel = get_isotropic("steel-generic")
    assert steel.youngs_modulus == 210_000.0
    assert steel.poisson_ratio == 0.30


def test_unknown_material_lists_known_ids() -> None:
    with pytest.raises(InputError, match="steel-generic"):
        get_isotropic("unobtainium")
    assert "aluminium-6061-t6" in material_ids()
