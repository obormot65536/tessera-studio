import numpy as np
import pytest

from tessera.post import von_mises


def test_uniaxial_and_pure_shear() -> None:
    s = np.array([[100.0, 0, 0, 0, 0, 0], [0, 0, 0, 50.0, 0, 0], [7.0, 7.0, 7.0, 0, 0, 0]])
    np.testing.assert_allclose(von_mises(s), [100.0, 50.0 * np.sqrt(3), 0.0], atol=1e-12)


def test_keeps_the_leading_axes() -> None:
    # (M, Q, 6) integration-point stresses give (M, Q), not a transposed (Q, M) (ADR-008).
    s = np.random.default_rng(1).normal(size=(3, 4, 6))
    vm = von_mises(s)
    assert vm.shape == (3, 4)
    np.testing.assert_allclose(vm, von_mises(s.reshape(12, 6)).reshape(3, 4))
    assert von_mises(s[1, 2]).shape == ()
    assert von_mises(s[1, 2]) == pytest.approx(vm[1, 2])


def test_rejects_wrong_shape() -> None:
    with pytest.raises(ValueError):
        von_mises(np.zeros((3, 5)))
