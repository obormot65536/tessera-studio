"""Solver-neutral finite element mesh (ADR-004).

Bulk numeric data uses numpy arrays in frozen dataclasses rather than pydantic models.
"""

from dataclasses import dataclass, field

import numpy as np
from numpy.typing import NDArray

# Quadratic tetrahedron in CalculiX/Abaqus C3D10 order (0-based): corners 0-3, then midside
# nodes on edges (0,1), (1,2), (2,0), (0,3), (1,3), (2,3).
TET10_EDGES: tuple[tuple[int, int], ...] = ((0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3))

# Quadratic triangle: corners 0-2, then midside nodes on edges (0,1), (1,2), (2,0).
TRI6_EDGES: tuple[tuple[int, int], ...] = ((0, 1), (1, 2), (2, 0))


@dataclass(frozen=True)
class Mesh:
    nodes: NDArray[np.float64]
    """(N, 3) coordinates in mm."""
    elements: NDArray[np.int64]
    """(M, 10) quadratic tetrahedra, 0-based node indices, C3D10 order."""
    element_sets: dict[str, NDArray[np.int64]] = field(default_factory=dict)
    """Named sets of element indices (e.g. one per part or print region)."""
    node_sets: dict[str, NDArray[np.int64]] = field(default_factory=dict)
    """Named sets of node indices."""
    surfaces: dict[str, NDArray[np.int64]] = field(default_factory=dict)
    """Named (K, 6) quadratic triangle faces in TRI6 order, e.g. one per geometry face."""

    def __post_init__(self) -> None:
        if self.nodes.ndim != 2 or self.nodes.shape[1] != 3:
            raise ValueError(f"nodes must have shape (N, 3), got {self.nodes.shape}")
        if self.elements.ndim != 2 or self.elements.shape[1] != 10:
            raise ValueError(f"elements must have shape (M, 10), got {self.elements.shape}")

    @property
    def n_nodes(self) -> int:
        return int(self.nodes.shape[0])

    @property
    def n_elements(self) -> int:
        return int(self.elements.shape[0])
