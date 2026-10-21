from collections.abc import Iterator

import pytest

import tessera.geometry as geometry


@pytest.fixture(autouse=True)
def _fresh_face_cache() -> Iterator[None]:
    """Every test starts and ends with no remembered faces, so tests cannot see each other's."""
    geometry.clear_face_cache()
    yield
    geometry.clear_face_cache()
