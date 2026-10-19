"""Analysis definitions (static now; modal, thermal, CFD later)."""

from tessera.physics.loads import tri6_load_weights, uniform_face_load
from tessera.physics.static import build_static_job

__all__ = ["build_static_job", "tri6_load_weights", "uniform_face_load"]
