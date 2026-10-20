"""Analysis definitions (static now; modal, thermal, CFD later)."""

from tessera.physics.loads import (
    body_force_load,
    orient_outward,
    pressure_face_load,
    tet10_volume_weights,
    tri6_load_weights,
    tri6_normal_weights,
    uniform_face_load,
)
from tessera.physics.static import build_static_job

__all__ = [
    "body_force_load",
    "build_static_job",
    "orient_outward",
    "pressure_face_load",
    "tet10_volume_weights",
    "tri6_load_weights",
    "tri6_normal_weights",
    "uniform_face_load",
]
