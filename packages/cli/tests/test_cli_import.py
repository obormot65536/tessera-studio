import importlib


def test_package_imports() -> None:
    assert importlib.import_module("tessera.cli") is not None
