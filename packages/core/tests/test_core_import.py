import importlib


def test_package_imports() -> None:
    assert importlib.import_module("tessera.core") is not None


def test_core_carries_the_thin_wall_minimum_for_clients() -> None:
    # ADR-007 Amendment 2: clients may import only core and model, and need N to advise.
    from tessera.core import WALL_ELEMENTS_MIN
    from tessera.post import WALL_ELEMENTS_MIN as POST_MINIMUM

    assert WALL_ELEMENTS_MIN == POST_MINIMUM == 5.0
