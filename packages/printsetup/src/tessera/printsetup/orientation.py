"""Material frame from a print setup (ADR-009 D2).

Material axes: 1 = raster/extrusion, 2 = across the road in the layer, 3 = build direction.
"""

import math

from tessera.model import MaterialFrame
from tessera.model.project import Vec3


def material_frame(build_direction: Vec3, raster_angle: float = 0.0) -> MaterialFrame:
    """Material frame of a part printed with ``build_direction`` up and rasters at ``raster_angle``.

    ``build_direction`` is the part-coordinate vector pointing up on the printer (any non-zero
    length). Axis 3 is that vector, normalized. Axis 1 is printer X carried by the smallest
    rotation taking printer +Z onto axis 3 (Rodrigues; a part upside down, axis 3 = -Z, uses
    180 degrees about printer X), then turned about axis 3 by ``raster_angle`` degrees,
    counter-clockwise looking down axis 3 (right-handed). Axis 2 = axis 3 x axis 1.

    Raises ``ValueError`` for a zero or non-finite build direction or a non-finite raster angle.
    """
    x, y, z = (float(c) for c in build_direction)
    if not all(math.isfinite(c) for c in (x, y, z)):
        raise ValueError(f"build direction must be finite, got {build_direction!r}")
    if not math.isfinite(raster_angle):
        raise ValueError(f"raster angle must be finite, got {raster_angle!r}")
    # Scale first so hypot neither overflows nor underflows for extreme magnitudes.
    scale = max(abs(x), abs(y), abs(z))
    if scale == 0.0:
        raise ValueError(f"build direction must be non-zero, got {build_direction!r}")
    x, y, z = x / scale, y / scale, z / scale
    norm = math.hypot(x, y, z)
    e3 = (x / norm, y / norm, z / norm)

    # Rodrigues for +Z -> e3 about k = Z x e3 = (-e3y, e3x, 0), applied to (1, 0, 0):
    #   e1 = (1 - cx^2 (1 - e3z), -cx cy (1 - e3z), -e3x),  (cx, cy) = (e3x, e3y) / rho.
    # Writing (1 - cos) / sin^2 as 1 / (1 + e3z) would blow up near -Z; the unit tilt direction
    # (cx, cy) with (1 - e3z) stays exact. At rho = 0, (cx, cy) = (0, 1) gives the identity for
    # +Z and the 180 degree turn about X for -Z (the ADR convention).
    rho = math.hypot(e3[0], e3[1])
    cx, cy = (e3[0] / rho, e3[1] / rho) if rho > 0.0 else (0.0, 1.0)
    k = 1.0 - e3[2]
    t1 = (1.0 - cx * cx * k, -cx * cy * k, -e3[0])
    t2 = _cross(e3, t1)

    a = math.radians(raster_angle)
    ca, sa = math.cos(a), math.sin(a)
    e1 = tuple(ca * p + sa * q for p, q in zip(t1, t2, strict=True))
    e2 = tuple(ca * q - sa * p for p, q in zip(t1, t2, strict=True))
    return MaterialFrame((e1[0], e1[1], e1[2]), (e2[0], e2[1], e2[2]))


def _cross(a: Vec3, b: Vec3) -> Vec3:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )
