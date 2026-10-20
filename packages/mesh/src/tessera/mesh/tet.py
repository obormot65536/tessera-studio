"""Quadratic tetrahedral meshing with gmsh."""

import math
from pathlib import Path

import gmsh
import numpy as np

from tessera.mesh._gmsh import gmsh_model
from tessera.model import InputError, Mesh, TriangleSurface

_GMSH_TRI3 = 2
_GMSH_TET10 = 11
_GMSH_TRI6 = 9
# gmsh orders the last two tet10 midside nodes (2,3), (1,3); C3D10 wants (1,3), (2,3).
_GMSH_TO_C3D10 = np.array([0, 1, 2, 3, 4, 5, 6, 7, 9, 8])


def _mesh_options(size: float | None) -> dict[str, float | str]:
    options: dict[str, float | str] = {"Mesh.ElementOrder": 2}
    if size is not None:
        options |= {"Mesh.MeshSizeMax": size, "Mesh.MeshSizeMin": size / 4}
    return options


def mesh_step(path: Path, *, part_id: str, size: float | None = None) -> Mesh:
    """Mesh a single-solid STEP file with C3D10 tetrahedra.

    The result has element set ``part:<part_id>`` and, for every geometry face ``i``, a
    surface ``face:<i>`` (tri6) and a node set ``face:<i>``.
    """
    if not path.is_file():
        raise FileNotFoundError(path)
    with gmsh_model(_mesh_options(size)):
        try:
            gmsh.model.occ.importShapes(str(path))
        except Exception as exc:  # gmsh raises plain Exception
            raise InputError(f"cannot read {path} as STEP: {exc}") from exc
        gmsh.model.occ.synchronize()
        if len(gmsh.model.getEntities(3)) != 1:
            raise InputError(f"{path} must contain exactly one solid in v0.1")
        gmsh.model.mesh.generate(3)
        return _extract(part_id, {face: face for _, face in gmsh.model.getEntities(2)})


def mesh_surface(surface: TriangleSurface, *, part_id: str, size: float | None = None) -> Mesh:
    """Mesh the solid bounded by a closed triangulated surface (ADR-005).

    gmsh rebuilds and remeshes the surface in patches; each patch is mapped back to the
    Tessera face its original triangles belong to, so ``face:<i>`` matches
    ``describe_surface`` for every face id ``i``.
    """
    # STL geometry is flat facets: midside nodes belong on straight edges, not on gmsh's
    # smooth reparametrized patches (which would bend flat faces).
    with gmsh_model(_mesh_options(size) | {"Mesh.SecondOrderLinear": 1}):
        entity = gmsh.model.addDiscreteEntity(2)
        node_tags = np.arange(1, surface.vertices.shape[0] + 1)
        gmsh.model.mesh.addNodes(2, entity, node_tags.tolist(), surface.vertices.ravel().tolist())
        gmsh.model.mesh.addElementsByType(
            entity, _GMSH_TRI3, [], (surface.triangles + 1).ravel().tolist()
        )
        # Split at every Tessera face boundary (and possibly more): patches never straddle faces.
        angle = math.radians(_feature_angle(surface)) * 0.999
        gmsh.model.mesh.classifySurfaces(angle, True, True, math.pi)
        gmsh.model.mesh.createGeometry()

        face_of_triangle = {
            tuple(sorted(tri)): int(face)
            for tri, face in zip((surface.triangles + 1).tolist(), surface.face_ids, strict=True)
        }
        patches = [tag for _, tag in gmsh.model.getEntities(2)]
        face_of_entity = {patch: _face_of_patch(patch, face_of_triangle) for patch in patches}

        loop = gmsh.model.geo.addSurfaceLoop(patches)
        gmsh.model.geo.addVolume([loop])
        gmsh.model.geo.synchronize()
        gmsh.model.mesh.generate(3)
        return _extract(part_id, face_of_entity)


def _feature_angle(surface: TriangleSurface) -> float:
    """Smallest dihedral angle (degrees) across a Tessera face boundary."""
    x = surface.vertices[surface.triangles]
    normals = np.cross(x[:, 1] - x[:, 0], x[:, 2] - x[:, 0])
    normals /= np.linalg.norm(normals, axis=1)[:, None]
    edges = np.sort(surface.triangles[:, [0, 1, 1, 2, 2, 0]].reshape(-1, 2), axis=1)
    owner = np.repeat(np.arange(surface.triangles.shape[0]), 3)
    order = np.lexsort((edges[:, 1], edges[:, 0]))
    edges, owner = edges[order], owner[order]
    pair = np.all(edges[:-1] == edges[1:], axis=1)
    a, b = owner[:-1][pair], owner[1:][pair]
    boundary = surface.face_ids[a] != surface.face_ids[b]
    if not boundary.any():
        return 180.0
    cos = np.clip(np.einsum("ij,ij->i", normals[a[boundary]], normals[b[boundary]]), -1, 1)
    return float(np.degrees(np.arccos(cos)).min())


def _face_of_patch(patch: int, face_of_triangle: dict[tuple[int, ...], int]) -> int:
    _, tri_nodes = gmsh.model.mesh.getElementsByType(_GMSH_TRI3, patch)
    faces = {
        face_of_triangle[tuple(sorted(t))] for t in np.asarray(tri_nodes).reshape(-1, 3).tolist()
    }
    if len(faces) != 1:
        raise RuntimeError(f"gmsh patch {patch} spans Tessera faces {sorted(faces)}")
    return faces.pop()


def _extract(part_id: str, face_of_entity: dict[int, int]) -> Mesh:
    """Read the generated gmsh mesh; group surface entities into Tessera faces."""
    raw_tags, raw_coords, _ = gmsh.model.mesh.getNodes()
    tags, coords = np.asarray(raw_tags, dtype=np.int64), np.asarray(raw_coords)
    index = np.full(int(tags.max()) + 1, -1, dtype=np.int64)
    index[tags] = np.arange(tags.size)
    nodes = coords.reshape(-1, 3).astype(np.float64)

    _, tet_nodes = gmsh.model.mesh.getElementsByType(_GMSH_TET10)
    elements = index[np.asarray(tet_nodes, dtype=np.int64)].reshape(-1, 10)[:, _GMSH_TO_C3D10]

    grouped: dict[int, list[np.ndarray]] = {}
    for entity, face in face_of_entity.items():
        _, tri_nodes = gmsh.model.mesh.getElementsByType(_GMSH_TRI6, entity)
        grouped.setdefault(face, []).append(
            index[np.asarray(tri_nodes, dtype=np.int64)].reshape(-1, 6)
        )
    surfaces = {f"face:{face}": np.concatenate(parts) for face, parts in sorted(grouped.items())}
    node_sets = {name: np.unique(tris) for name, tris in surfaces.items()}

    return Mesh(
        nodes=nodes,
        elements=elements,
        element_sets={f"part:{part_id}": np.arange(elements.shape[0], dtype=np.int64)},
        node_sets=node_sets,
        surfaces=surfaces,
    )
