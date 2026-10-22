"""Results post-processing."""

from tessera.post.flags import (
    FLAG_DISTANCE_FACTOR,
    INSIDE_EDGE_ANGLE,
    flag_nodes,
    flag_points,
    median_edge_length,
    sharp_inside_edge_nodes,
    support_is_singular,
)
from tessera.post.mesh_check import (
    MESH_CHECK_REFINEMENT_LIMIT,
    MESH_CHECK_THRESHOLD,
    compare_meshes,
    local_edge_length,
    spot_check,
)
from tessera.post.points import element_nodal_stresses, point_von_mises
from tessera.post.summary import summarize, von_mises

__all__ = [
    "FLAG_DISTANCE_FACTOR",
    "INSIDE_EDGE_ANGLE",
    "MESH_CHECK_REFINEMENT_LIMIT",
    "MESH_CHECK_THRESHOLD",
    "compare_meshes",
    "element_nodal_stresses",
    "flag_nodes",
    "flag_points",
    "local_edge_length",
    "median_edge_length",
    "point_von_mises",
    "sharp_inside_edge_nodes",
    "spot_check",
    "summarize",
    "support_is_singular",
    "von_mises",
]
