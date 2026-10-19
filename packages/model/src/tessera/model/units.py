"""Internal unit system (ADR-002).

All quantities inside Tessera use the consistent mm-N-MPa system:

| Quantity | Unit |
|---|---|
| length | mm |
| force | N |
| stress, pressure, modulus | MPa (N/mm^2) |
| mass | tonne (t) |
| density | t/mm^3 |
| time | s |

Conversion happens only at the edges (file import, UI, reports) through
``to_internal`` and ``from_internal``.
"""

from typing import Final

LENGTH_UNIT: Final = "mm"
FORCE_UNIT: Final = "N"
STRESS_UNIT: Final = "MPa"

_LENGTH_TO_MM: Final[dict[str, float]] = {"mm": 1.0, "cm": 10.0, "m": 1000.0, "in": 25.4}
_FORCE_TO_N: Final[dict[str, float]] = {"N": 1.0, "kN": 1000.0, "lbf": 4.4482216152605}
_STRESS_TO_MPA: Final[dict[str, float]] = {"MPa": 1.0, "Pa": 1e-6, "kPa": 1e-3, "GPa": 1000.0}

_FACTORS: Final[dict[str, dict[str, float]]] = {
    "length": _LENGTH_TO_MM,
    "force": _FORCE_TO_N,
    "stress": _STRESS_TO_MPA,
}


def _factor(kind: str, unit: str) -> float:
    try:
        return _FACTORS[kind][unit]
    except KeyError:
        known = ", ".join(_FACTORS.get(kind, {}))
        raise ValueError(f"unknown {kind} unit {unit!r}; known: {known}") from None


def to_internal(value: float, kind: str, unit: str) -> float:
    """Convert ``value`` given in ``unit`` to the internal unit for ``kind``."""
    return value * _factor(kind, unit)


def from_internal(value: float, kind: str, unit: str) -> float:
    """Convert an internal ``value`` of ``kind`` to ``unit``."""
    return value / _factor(kind, unit)
