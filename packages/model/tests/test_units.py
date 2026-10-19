import pytest
from hypothesis import given
from hypothesis import strategies as st

from tessera.model.units import from_internal, to_internal


def test_known_conversions() -> None:
    assert to_internal(1.0, "length", "in") == pytest.approx(25.4)
    assert to_internal(2.0, "force", "kN") == pytest.approx(2000.0)
    assert to_internal(70.0, "stress", "GPa") == pytest.approx(70_000.0)


@given(
    st.floats(min_value=-1e9, max_value=1e9, allow_nan=False), st.sampled_from(["mm", "m", "in"])
)
def test_round_trip(value: float, unit: str) -> None:
    assert from_internal(to_internal(value, "length", unit), "length", unit) == pytest.approx(value)


def test_unknown_unit_is_rejected() -> None:
    with pytest.raises(ValueError, match="unknown length unit"):
        to_internal(1.0, "length", "furlong")
