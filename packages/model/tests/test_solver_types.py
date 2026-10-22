"""Solve-job types for orthotropic materials and orientation (ADR-009 PR-1)."""

import dataclasses
import math

import numpy as np
import pytest

import tessera.model as model
from tessera.model.mesh import Mesh
from tessera.model.solver import (
    IsotropicMaterial,
    MaterialFrame,
    OrthotropicMaterial,
    Section,
    SolveJob,
)

TO = {
    "e1": 3000.0,
    "e2": 2400.0,
    "e3": 1800.0,
    "nu12": 0.35,
    "nu13": 0.30,
    "nu23": 0.30,
    "g12": 1000.0,
    "g13": 700.0,
    "g23": 650.0,
}


def test_new_names_are_exported() -> None:
    for name in ("OrthotropicMaterial", "Material", "MaterialFrame"):
        assert name in model.__all__
        assert hasattr(model, name)


def test_orthotropic_material_holds_constants_and_is_frozen() -> None:
    mat = OrthotropicMaterial("to", **TO, density=1.24e-9)
    assert (mat.e1, mat.nu23, mat.g23, mat.density) == (3000.0, 0.30, 650.0, 1.24e-9)
    assert OrthotropicMaterial("to", **TO).density is None
    with pytest.raises(dataclasses.FrozenInstanceError):
        mat.e1 = 1.0  # type: ignore[misc]


def test_orthotropic_material_documents_major_poisson_convention() -> None:
    doc = OrthotropicMaterial.__doc__ or ""
    assert "-eps_j / eps_i" in doc
    assert "nu_ji = nu_ij E_j / E_i" in doc


@pytest.mark.parametrize(
    "changes",
    [{"e1": 0.0}, {"g13": -5.0}, {"g23": math.nan}, {"nu12": 1.2, "e2": 3000.0}],
)
def test_orthotropic_material_is_validated(changes: dict[str, float]) -> None:
    with pytest.raises(ValueError):
        OrthotropicMaterial("bad", **{**TO, **changes})


def test_frame_axes_and_rotation() -> None:
    frame = MaterialFrame((0.0, 0.0, -1.0), (0.0, 1.0, 0.0))
    assert frame.axis_3 == pytest.approx((1.0, 0.0, 0.0))
    np.testing.assert_array_equal(frame.rotation, [[0, 0, -1], [0, 1, 0], [1, 0, 0]])
    assert np.linalg.det(frame.rotation) == pytest.approx(1.0)


def test_identity_frame() -> None:
    frame = MaterialFrame((1.0, 0.0, 0.0), (0.0, 1.0, 0.0))
    assert frame.axis_3 == (0.0, 0.0, 1.0)
    np.testing.assert_array_equal(frame.rotation, np.eye(3))


@pytest.mark.parametrize(
    ("axis_1", "axis_2", "match"),
    [
        ((1.0, 0.0, 0.0), (0.0, 1.0 + 1e-8, 0.0), "unit"),
        ((2.0, 0.0, 0.0), (0.0, 1.0, 0.0), "unit"),
        ((0.0, 0.0, 0.0), (0.0, 1.0, 0.0), "unit"),
        ((1.0, 0.0, 0.0), (math.sqrt(0.5), math.sqrt(0.5), 0.0), "orthogonal"),
        ((1.0, 0.0, 0.0), (1e-8, 1.0, 0.0), "orthogonal"),
        ((math.nan, 0.0, 0.0), (0.0, 1.0, 0.0), "finite"),
    ],
)
def test_frame_rejects_non_orthonormal_axes(
    axis_1: tuple[float, float, float], axis_2: tuple[float, float, float], match: str
) -> None:
    with pytest.raises(ValueError, match=match):
        MaterialFrame(axis_1, axis_2)


def test_frame_accepts_rounding_within_tolerance() -> None:
    s = math.sqrt(0.5)
    MaterialFrame((s, s, 0.0), (-s, s, 1e-12))


def test_from_axes_gram_schmidt_and_right_handed() -> None:
    a, b = (2.0, 0.0, 0.0), (5.0, 3.0, 0.0)
    frame = MaterialFrame.from_axes(a, b)
    assert frame.axis_1 == pytest.approx((1.0, 0.0, 0.0))
    assert frame.axis_2 == pytest.approx((0.0, 1.0, 0.0))
    assert frame.axis_3 == pytest.approx((0.0, 0.0, 1.0))

    tilted = MaterialFrame.from_axes((1.0, 2.0, -0.5), (0.3, -1.0, 2.0))
    r = tilted.rotation
    np.testing.assert_allclose(r @ r.T, np.eye(3), atol=1e-12)
    assert np.linalg.det(r) == pytest.approx(1.0)
    np.testing.assert_allclose(np.cross(r[0], r[1]), r[2], atol=1e-12)
    assert all(isinstance(x, float) for x in tilted.axis_1 + tilted.axis_2)


@pytest.mark.parametrize(
    ("a", "b"),
    [((0.0, 0.0, 0.0), (0.0, 1.0, 0.0)), ((1.0, 0.0, 0.0), (-3.0, 0.0, 0.0))],
)
def test_from_axes_rejects_zero_or_parallel(
    a: tuple[float, float, float], b: tuple[float, float, float]
) -> None:
    with pytest.raises(ValueError):
        MaterialFrame.from_axes(a, b)


def test_section_orientation_defaults_to_global_axes() -> None:
    assert Section("part", "pla").orientation is None
    frame = MaterialFrame((1.0, 0.0, 0.0), (0.0, 1.0, 0.0))
    assert Section("part", "pla", orientation=frame).orientation is frame


def test_solve_job_takes_mixed_materials() -> None:
    mesh = Mesh(
        nodes=np.zeros((0, 3)),
        elements=np.zeros((0, 10), dtype=np.int64),
        element_sets={},
        node_sets={},
    )
    iso = IsotropicMaterial("steel", 210000.0, 0.3)
    ortho = OrthotropicMaterial("to", **TO)
    job = SolveJob(mesh=mesh, materials=(iso, ortho), sections=())
    assert job.materials == (iso, ortho)
