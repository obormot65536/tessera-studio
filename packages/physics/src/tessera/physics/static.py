"""Build a discrete static ``SolveJob`` from a study on a meshed part."""

import dataclasses
from typing import assert_never

import numpy as np

from tessera.model import (
    FaceSet,
    FixedSupport,
    ForceLoad,
    InputError,
    IsotropicMaterial,
    Mesh,
    MomentLoad,
    NodalLoads,
    PartInstance,
    PressureLoad,
    Section,
    SolveJob,
    Study,
)
from tessera.physics.loads import (
    body_force_load,
    moment_face_load,
    orient_outward,
    pressure_face_load,
    uniform_face_load,
)


def _surface(mesh: Mesh, faces: FaceSet, part: PartInstance) -> np.ndarray:
    if faces.part_id != part.id:
        raise InputError(f"faces of part {faces.part_id!r} used in a study of part {part.id!r}")
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
        node_sets[name] = np.unique(_surface(mesh, support.faces, part))
        supports.append(FixedSupport(node_set=name, dofs=support.dofs))

    elements = mesh.elements[mesh.element_sets[f"part:{part.id}"]]
    loads = []
    for load in study.loads:
        tris = _surface(mesh, load.faces, part)
        match load:
            case ForceLoad():
                force = np.asarray(load.force, dtype=float)
                node_ids, forces = uniform_face_load(mesh.nodes, tris, force)
            case PressureLoad():
                tris = orient_outward(mesh.nodes, elements, tris)
                node_ids, forces = pressure_face_load(mesh.nodes, tris, load.pressure)
            case MomentLoad():
                moment = np.asarray(load.moment, dtype=float)
                node_ids, forces = moment_face_load(mesh.nodes, tris, moment)
            case _:
                assert_never(load)
        loads.append(NodalLoads(nodes=node_ids, forces=forces))
    if study.gravity is not None:
        if material.density is None:
            raise InputError(f"material {material.name!r} has no density, needed for gravity")
        node_ids, forces = body_force_load(mesh.nodes, elements, material.density, study.gravity)
        loads.append(NodalLoads(nodes=node_ids, forces=forces))

    return SolveJob(
        mesh=dataclasses.replace(mesh, node_sets=node_sets),
        materials=(material,),
        sections=(Section(element_set=f"part:{part.id}", material=material.name),),
        supports=tuple(supports),
        loads=tuple(loads),
    )
