"""Node sets on geometry edges, for line supports (ADR-006)."""

import numpy as np
from numpy.typing import NDArray

from tessera.model import TRI6_EDGES, InputError, Mesh

# TRI6 midside node of each element edge, in TRI6_EDGES order.
_MIDSIDE = (3, 4, 5)


def _element_edges(tris: NDArray[np.int64]) -> NDArray[np.int64]:
    """(3K, 3) element edges of tri6 faces: two corners (ascending) and the midside node."""
    return np.concatenate(
        [
            np.column_stack(
                [
                    np.minimum(tris[:, i], tris[:, j]),
                    np.maximum(tris[:, i], tris[:, j]),
                    tris[:, m],
                ]
            )
            for (i, j), m in zip(TRI6_EDGES, _MIDSIDE, strict=True)
        ]
    )


def edge_nodes(mesh: Mesh, part_id: str, a: int, b: int) -> NDArray[np.int64]:
    """Sorted nodes on the edge between faces ``a`` and ``b``: corners and midsides.

    These are the nodes of the tri6 element edges that occur in both ``face:<a>`` and
    ``face:<b>``. Faces that only touch at a corner share no element edge, so that is refused
    like faces that do not meet at all.
    """
    surfaces = []
    for face in (a, b):
        try:
            surfaces.append(mesh.surfaces[f"face:{face}"])
        except KeyError:
            raise InputError(f"part {part_id!r} has no face {face}") from None
    first, second = (_element_edges(s) for s in surfaces)
    # One integer per element edge, from its two corner nodes.
    key_first, key_second = (e[:, 0] * mesh.n_nodes + e[:, 1] for e in (first, second))
    on_edge = first[np.isin(key_first, key_second)]
    if on_edge.size == 0:
        raise InputError(
            f"faces {a} and {b} of part {part_id!r} do not meet along an edge "
            "(list the edges with `tessera edges`)"
        )
    return np.unique(on_edge)
