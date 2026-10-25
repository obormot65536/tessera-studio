"""Material strengths for the failure criteria (ADR-010 D3).

Strengths are in MPa, compressive ones as positive magnitudes. They describe first failure
(onset), as measured on coupons; the criteria that use them live in ``tessera.post``.
"""

import math
from dataclasses import dataclass, fields


def _check_positive(owner: object) -> None:
    for f in fields(owner):  # type: ignore[arg-type]
        value = getattr(owner, f.name)
        if value is None or isinstance(value, bool):
            continue
        if not (math.isfinite(value) and value > 0.0):
            raise ValueError(f"{f.name} must be a positive finite strength in MPa, got {value}")


@dataclass(frozen=True)
class OrthotropicStrength:
    """Strengths on the material axes (1 = raster, 2 = across the road in the layer, 3 = build).

    ``t*`` are tensile and ``c*`` compressive strengths; ``s12`` is in-layer shear and ``s13``,
    ``s23`` interlayer shear. A missing compressive strength is taken equal to the tensile one
    and listed in ``assumed`` (ADR-010 D3); shear strengths are required.
    """

    t1: float
    t2: float
    t3: float
    s12: float
    s13: float
    s23: float
    c1: float | None = None
    c2: float | None = None
    c3: float | None = None

    def __post_init__(self) -> None:
        _check_positive(self)

    @property
    def compressive(self) -> tuple[float, float, float]:
        """(C1, C2, C3), each the tensile value where it was not given."""
        return (
            self.t1 if self.c1 is None else self.c1,
            self.t2 if self.c2 is None else self.c2,
            self.t3 if self.c3 is None else self.c3,
        )

    @property
    def assumed(self) -> tuple[str, ...]:
        """Names of the compressive strengths taken equal to the tensile ones."""
        return tuple(name for name in ("c1", "c2", "c3") if getattr(self, name) is None)


@dataclass(frozen=True)
class IsotropicStrength:
    """Strength of an isotropic material.

    For a ductile material (``brittle=False``, metals) ``tensile`` is the yield strength and the
    criterion is von Mises; a compressive strength is refused there because nothing would use
    it. For a brittle one (resin) ``tensile`` is the tensile strength and the criterion is the
    maximum principal stress; a missing ``compressive`` is taken equal to ``tensile`` and listed
    in ``assumed`` (ADR-010 D2, D3).
    """

    tensile: float
    compressive: float | None = None
    brittle: bool = False

    def __post_init__(self) -> None:
        _check_positive(self)
        if self.compressive is not None and not self.brittle:
            raise ValueError(
                "compressive strength is used only for brittle materials; a ductile material "
                "is checked with von Mises against its yield strength (tensile)"
            )

    @property
    def compressive_or_tensile(self) -> float:
        return self.tensile if self.compressive is None else self.compressive

    @property
    def assumed(self) -> tuple[str, ...]:
        return ("compressive",) if self.brittle and self.compressive is None else ()


Strength = IsotropicStrength | OrthotropicStrength
