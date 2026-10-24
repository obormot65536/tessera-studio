"""build_static_job on a hand-built one-element mesh."""

import dataclasses
from pathlib import Path

import numpy as np
import pytest

from tessera.model import (
    TET10_EDGES,
    EdgeSet,
    FaceSet,
    ForceLoad,
    InputError,
    IsotropicMaterial,
    MaterialFrame,
    Mesh,
    OrthotropicMaterial,
    PartInstance,
    Section,
    Study,
    Support,
)
from tessera.physics import build_static_job, edge_nodes

_CORNERS = np.array([[0.0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
_NODES = np.vstack([_CORNERS, [0.5 * (_CORNERS[i] + _CORNERS[j]) for i, j in TET10_EDGES]])
_MESH = Mesh(
    nodes=_NODES,
    elements=np.arange(10)[None, :],
    element_sets={"part:p": np.array([0])},
    node_sets={"face:1": np.array([0, 1, 2, 4, 5, 6])},
    surfaces={"face:1": np.array([[0, 2, 1, 6, 5, 4]])},
)
_PART = PartInstance(id="p", name="p", geometry=Path("p.step"))
_STEEL = IsotropicMaterial(name="steel", youngs_modulus=210_000.0, poisson_ratio=0.3)


def _study(support_part: str, load_part: str) -> Study:
    face = (1,)
    return Study(
        supports=[Support(id="s", faces=FaceSet(part_id=support_part, faces=face))],
        loads=[ForceLoad(id="f", faces=FaceSet(part_id=load_part, faces=face), force=(0, 0, 1))],
    )


def test_faces_of_the_part_build_a_job() -> None:
    job = build_static_job(_study("p", "p"), _PART, _MESH, _STEEL)
    assert job.loads[0].forces.sum(axis=0) == pytest.approx([0.0, 0.0, 1.0])


@pytest.mark.parametrize(("support_part", "load_part"), [("other", "p"), ("p", "other")])
def test_faces_of_another_part_are_refused(support_part: str, load_part: str) -> None:
    with pytest.raises(InputError, match="part 'other'"):
        build_static_job(_study(support_part, load_part), _PART, _MESH, _STEEL)


# Face 2 is the y = 0 face (corners 0, 1, 3); it meets face 1 (z = 0) along corners 0-1, whose
# midside node is 4. Face 3 is a stray triangle that touches face 1 only at corner 0.
_EXTRA = np.array([[-1.0, -1, -1], [-2, -1, -1], [-1, -2, -1], [-1.5, -1, -1], [-1.5, -1.5, -1]])
_EDGE_MESH = Mesh(
    nodes=np.vstack([_NODES, _EXTRA]),
    elements=np.arange(10)[None, :],
    element_sets={"part:p": np.array([0])},
    surfaces={
        "face:1": np.array([[0, 2, 1, 6, 5, 4]]),
        "face:2": np.array([[0, 1, 3, 4, 8, 7]]),
        "face:3": np.array([[0, 11, 12, 10, 13, 14]]),
    },
)


def test_edge_nodes_are_the_corners_and_midside_of_shared_element_edges() -> None:
    assert edge_nodes(_EDGE_MESH, "p", 1, 2).tolist() == [0, 1, 4]
    assert edge_nodes(_EDGE_MESH, "p", 2, 1).tolist() == [0, 1, 4]


def test_faces_touching_at_a_corner_have_no_edge() -> None:
    with pytest.raises(InputError, match="do not meet along an edge"):
        edge_nodes(_EDGE_MESH, "p", 1, 3)
    with pytest.raises(InputError, match="no face 7"):
        edge_nodes(_EDGE_MESH, "p", 1, 7)


def test_a_line_support_holds_the_edge_nodes_and_is_marked() -> None:
    line = Support(id="line", edges=EdgeSet(part_id="p", edges=((1, 2),)), dofs=(3,))
    face = Support(id="face", faces=FaceSet(part_id="p", faces=(1,)))
    job = build_static_job(Study(supports=[line, face]), _PART, _EDGE_MESH, _STEEL)
    assert job.mesh.node_sets["support:line"].tolist() == [0, 1, 4]
    assert [(s.node_set, s.dofs, s.line) for s in job.supports] == [
        ("support:line", (3,), True),
        ("support:face", (1, 2, 3), False),
    ]


def test_a_line_support_on_another_part_is_refused() -> None:
    line = Support(id="line", edges=EdgeSet(part_id="other", edges=((1, 2),)))
    with pytest.raises(InputError, match="part 'other'"):
        build_static_job(Study(supports=[line]), _PART, _EDGE_MESH, _STEEL)


# ADR-009 PR-6: an orthotropic material with a frame, and its density under gravity.
_TO = OrthotropicMaterial(
    name="to",
    e1=3000.0,
    e2=2400.0,
    e3=1800.0,
    nu12=0.35,
    nu13=0.30,
    nu23=0.30,
    g12=1000.0,
    g13=700.0,
    g23=650.0,
    density=1.2e-9,
)
_SIDE = MaterialFrame((0.0, 0.0, -1.0), (0.0, 1.0, 0.0))  # printed with build direction +x


def test_an_orthotropic_material_gets_its_frame_on_the_section() -> None:
    job = build_static_job(_study("p", "p"), _PART, _MESH, _TO, frame=_SIDE)
    assert job.materials == (_TO,)
    assert job.sections == (Section(element_set="part:p", material="to", orientation=_SIDE),)


def test_without_a_frame_the_section_uses_the_part_axes() -> None:
    job = build_static_job(_study("p", "p"), _PART, _MESH, _STEEL)
    assert job.sections[0].orientation is None


def test_gravity_uses_the_orthotropic_density() -> None:
    study = _study("p", "p")
    study.loads.clear()
    study.gravity = (0.0, 0.0, -9810.0)
    job = build_static_job(study, _PART, _MESH, _TO, frame=_SIDE)
    weight = 1.2e-9 * (1.0 / 6.0) * -9810.0  # density x volume of the unit tet x g
    assert job.loads[0].forces.sum(axis=0) == pytest.approx([0.0, 0.0, weight])
    no_density = dataclasses.replace(_TO, density=None)
    with pytest.raises(InputError, match="'to' has no density"):
        build_static_job(study, _PART, _MESH, no_density, frame=_SIDE)
