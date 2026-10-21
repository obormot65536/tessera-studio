"""build_static_job on a hand-built one-element mesh."""

from pathlib import Path

import numpy as np
import pytest

from tessera.model import (
    TET10_EDGES,
    FaceSet,
    ForceLoad,
    InputError,
    IsotropicMaterial,
    Mesh,
    PartInstance,
    Study,
    Support,
)
from tessera.physics import build_static_job

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
