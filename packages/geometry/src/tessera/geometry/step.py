"""STEP import and face queries (via gmsh's OpenCASCADE kernel)."""

from pathlib import Path

import gmsh
import numpy as np

from tessera.geometry._gmsh import gmsh_model
from tessera.model import FaceInfo


def _import_single_solid(path: Path) -> int:
    if not path.is_file():
        raise FileNotFoundError(path)
    gmsh.model.occ.importShapes(str(path))
    gmsh.model.occ.synchronize()
    volumes = gmsh.model.getEntities(3)
    if len(volumes) != 1:
        raise ValueError(f"{path} contains {len(volumes)} solids; v0.1 supports exactly one")
    return volumes[0][1]


def _outward_normal(face: int, volume: int, centroid: np.ndarray, size: float) -> np.ndarray:
    uv = gmsh.model.getParametrization(2, face, centroid.tolist())
    normal = np.asarray(gmsh.model.getNormal(face, uv), dtype=float)
    normal /= np.linalg.norm(normal)
    probe = centroid + normal * 1e-4 * size
    if gmsh.model.isInside(3, volume, probe.tolist()):
        normal = -normal
    return normal


def describe_faces(path: Path) -> list[FaceInfo]:
    """List the faces of a single-solid STEP file, ordered by face index."""
    with gmsh_model():
        volume = _import_single_solid(path)
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(3, volume)
        size = float(np.linalg.norm([xmax - xmin, ymax - ymin, zmax - zmin]))
        faces = []
        for _, tag in sorted(gmsh.model.getEntities(2)):
            kind = gmsh.model.getType(2, tag)
            centroid = np.asarray(gmsh.model.occ.getCenterOfMass(2, tag))
            normal = None
            if kind == "Plane":
                n = _outward_normal(tag, volume, centroid, size)
                normal = (float(n[0]), float(n[1]), float(n[2]))
            faces.append(
                FaceInfo(
                    index=tag,
                    kind=kind,
                    area=float(gmsh.model.occ.getMass(2, tag)),
                    centroid=(float(centroid[0]), float(centroid[1]), float(centroid[2])),
                    normal=normal,
                )
            )
        return faces
