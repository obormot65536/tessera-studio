import numpy as np
import pytest

from tessera.post import von_mises


def test_uniaxial_and_pure_shear() -> None:
    s = np.array([[100.0, 0, 0, 0, 0, 0], [0, 0, 0, 50.0, 0, 0], [7.0, 7.0, 7.0, 0, 0, 0]])
    np.testing.assert_allclose(von_mises(s), [100.0, 50.0 * np.sqrt(3), 0.0], atol=1e-12)


def test_rejects_wrong_shape() -> None:
    with pytest.raises(ValueError):
        von_mises(np.zeros((3, 5)))
