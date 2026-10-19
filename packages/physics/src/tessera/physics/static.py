"""Build a discrete static ``SolveJob`` from a study on a meshed part."""

import dataclasses

import numpy as np

from tessera.model import (
    FaceSet,
    FixedSupport,
    InputError,
    IsotropicMaterial,
    Mesh,
    NodalLoads,
    PartInstance,
    Section,
    SolveJob,
    Study,
)
from tessera.physics.loads import uniform_face_load


def _surface(mesh: Mesh, faces: FaceSet) -> np.ndarray:
    try:
        return np.concatenate([mesh.surfaces[f"face:{f}"] for f in faces.faces])
    except KeyError as missing:
        face = str(missing.args[0]).removeprefix("face:")
        raise InputError(f"part {faces.part_id!r} has no face {face}") from None


def build_static_job(
    study: Study, part: PartInstance, mesh: Mesh, material: IsotropicMaterial
) -> SolveJob:
    if not study.supports:
        raise InputError("the study needs at least one support")
    node_sets = dict(mesh.node_sets)
    supports = []
    for support in study.supports:
        name = f"support:{support.id}"
        node_sets[name] = np.unique(_surface(mesh, support.faces))
        supports.append(FixedSupport(node_set=name, dofs=support.dofs))

    loads = []
    for load in study.loads:
        node_ids, forces = uniform_face_load(
            mesh.nodes, _surface(mesh, load.faces), np.asarray(load.force, dtype=float)
        )
        loads.append(NodalLoads(nodes=node_ids, forces=forces))

    return SolveJob(
        mesh=dataclasses.replace(mesh, node_sets=node_sets),
        materials=(material,),
        sections=(Section(element_set=f"part:{part.id}", material=material.name),),
        supports=tuple(supports),
        loads=tuple(loads),
    )
