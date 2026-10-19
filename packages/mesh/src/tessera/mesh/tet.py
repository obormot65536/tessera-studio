"""Quadratic tetrahedral meshing with gmsh."""

from pathlib import Path

import gmsh
import numpy as np

from tessera.mesh._gmsh import gmsh_model
from tessera.model import Mesh

_GMSH_TET10 = 11
_GMSH_TRI6 = 9
# gmsh orders the last two tet10 midside nodes (2,3), (1,3); C3D10 wants (1,3), (2,3).
_GMSH_TO_C3D10 = np.array([0, 1, 2, 3, 4, 5, 6, 7, 9, 8])


def mesh_step(path: Path, *, part_id: str, size: float | None = None) -> Mesh:
    """Mesh a single-solid STEP file with C3D10 tetrahedra.

    The result has element set ``part:<part_id>`` and, for every geometry face ``i``, a
    surface ``face:<i>`` (tri6) and a node set ``face:<i>``.
    """
    if not path.is_file():
        raise FileNotFoundError(path)
    with gmsh_model():
        gmsh.model.occ.importShapes(str(path))
        gmsh.model.occ.synchronize()
        if len(gmsh.model.getEntities(3)) != 1:
            raise ValueError(f"{path} must contain exactly one solid in v0.1")
        if size is not None:
            gmsh.option.setNumber("Mesh.MeshSizeMax", size)
            gmsh.option.setNumber("Mesh.MeshSizeMin", size / 4)
        gmsh.option.setNumber("Mesh.ElementOrder", 2)
        gmsh.model.mesh.generate(3)

        raw_tags, raw_coords, _ = gmsh.model.mesh.getNodes()
        tags, coords = np.asarray(raw_tags, dtype=np.int64), np.asarray(raw_coords)
        index = np.full(int(tags.max()) + 1, -1, dtype=np.int64)
        index[tags] = np.arange(tags.size)
        nodes = coords.reshape(-1, 3).astype(np.float64)

        _, tet_nodes = gmsh.model.mesh.getElementsByType(_GMSH_TET10)
        elements = index[np.asarray(tet_nodes, dtype=np.int64)].reshape(-1, 10)[:, _GMSH_TO_C3D10]

        surfaces: dict[str, np.ndarray] = {}
        node_sets: dict[str, np.ndarray] = {}
        for _, face in gmsh.model.getEntities(2):
            _, tri_nodes = gmsh.model.mesh.getElementsByType(_GMSH_TRI6, face)
            tris = index[np.asarray(tri_nodes, dtype=np.int64)].reshape(-1, 6)
            surfaces[f"face:{face}"] = tris
            node_sets[f"face:{face}"] = np.unique(tris)

    return Mesh(
        nodes=nodes,
        elements=elements,
        element_sets={f"part:{part_id}": np.arange(elements.shape[0], dtype=np.int64)},
        node_sets=node_sets,
        surfaces=surfaces,
    )
