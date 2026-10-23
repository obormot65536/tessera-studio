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
    WALL_ELEMENTS_MIN,
    compare_meshes,
    local_edge_length,
    node_edge_length,
    spot_check,
)
from tessera.post.points import (
    DegenerateElementsWarning,
    element_nodal_stresses,
    point_von_mises,
)
from tessera.post.summary import summarize, von_mises
from tessera.post.wall import (
    FACE_CLUSTER_ANGLE,
    OuterSurface,
    outer_surface,
    samples_on,
    wall_span,
)

__all__ = [
    "FACE_CLUSTER_ANGLE",
    "FLAG_DISTANCE_FACTOR",
    "INSIDE_EDGE_ANGLE",
    "MESH_CHECK_REFINEMENT_LIMIT",
    "MESH_CHECK_THRESHOLD",
    "WALL_ELEMENTS_MIN",
    "DegenerateElementsWarning",
    "OuterSurface",
    "compare_meshes",
    "element_nodal_stresses",
    "flag_nodes",
    "flag_points",
    "local_edge_length",
    "median_edge_length",
    "node_edge_length",
    "outer_surface",
    "point_von_mises",
    "samples_on",
    "sharp_inside_edge_nodes",
    "spot_check",
    "summarize",
    "support_is_singular",
    "von_mises",
    "wall_span",
]
