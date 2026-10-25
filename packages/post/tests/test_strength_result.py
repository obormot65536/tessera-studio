"""Strength ratios in a study summary and the two-mesh check (ADR-010 D4, D5).

Small structured tet10 grids with straight edges and synthetic integration-point stresses, so
no mesher or solver runs.
"""

import itertools
from collections.abc import Callable

import numpy as np
import pytest
from numpy.typing import NDArray

from tessera.model import (
    TET10_EDGES,
    FailureCriterion,
    FailureMode,
    FixedSupport,
    IsotropicMaterial,
    IsotropicStrength,
    Material,
    MaterialFrame,
    Mesh,
    OrthotropicMaterial,
    OrthotropicStrength,
    PointStresses,
    Section,
    SolveJob,
    SolveResult,
    Verdict,
)
from tessera.post import compare_meshes, summarize

Field = Callable[[NDArray[np.float64]], NDArray[np.float64]]

STRENGTH = OrthotropicStrength(
    t1=40.0, t2=30.0, t3=15.0, s12=20.0, s13=12.0, s23=12.0, c1=50.0, c2=45.0, c3=40.0
)
STIFFNESS = dict(
    e1=3000.0, e2=2400.0, e3=1800.0, nu12=0.35, nu13=0.30, nu23=0.30, g12=1000.0, g13=700.0,
    g23=650.0,
)  # fmt: skip
TO = OrthotropicMaterial("to", **STIFFNESS, strength=STRENGTH)
SIDE = MaterialFrame((0.0, 0.0, -1.0), (0.0, 1.0, 0.0))  # printed on its side: axis 3 = +x
A, B = 0.5854101966249685, 0.1381966011250105
RULE = np.array([[A, B, B, B], [B, A, B, B], [B, B, A, B], [B, B, B, A]])
"""Barycentric coordinates of the four tet4 integration points."""


def grid(cells: tuple[int, int, int], size: float = 1.0) -> Mesh:
    """Straight-edged quadratic tets on a grid of cubes of side ``size``, six per cube."""
    nx, ny, nz = cells

    def index(i: int, j: int, k: int) -> int:
        return (i * (ny + 1) + j) * (nz + 1) + k

    corners = []
    for i, j, k in itertools.product(range(nx), range(ny), range(nz)):
        for perm in itertools.permutations(range(3)):
            v = [i, j, k]
            path = [index(*v)]
            for axis in perm:
                v[axis] += 1
                path.append(index(*v))
            corners.append(path)
    corners = np.array(corners, dtype=np.int64)
    xyz = size * np.stack(
        np.meshgrid(*(np.arange(n + 1) for n in cells), indexing="ij"), axis=-1
    ).reshape(-1, 3).astype(np.float64)
    a, b, c, d = (xyz[corners[:, n]] for n in range(4))
    negative = np.einsum("ij,ij->i", np.cross(b - a, c - a), d - a) < 0
    corners[negative] = corners[negative][:, [0, 2, 1, 3]]
    pairs = np.sort(corners[:, np.array(TET10_EDGES)], axis=2)
    edges, inverse = np.unique(pairs.reshape(-1, 2), axis=0, return_inverse=True)
    nodes = np.vstack([xyz, xyz[edges].mean(axis=1)])
    elements = np.hstack([corners, len(xyz) + inverse.reshape(-1, 6)])
    lo, hi = nodes.min(axis=0), nodes.max(axis=0)
    node_sets = {}
    for axis in range(3):
        node_sets[f"face:{2 * axis + 1}"] = np.flatnonzero(np.isclose(nodes[:, axis], lo[axis]))
        node_sets[f"face:{2 * axis + 2}"] = np.flatnonzero(np.isclose(nodes[:, axis], hi[axis]))
    return Mesh(
        nodes=nodes,
        elements=elements,
        element_sets={"part:p": np.arange(len(elements))},
        node_sets=node_sets,
    )


def split(mesh: Mesh, left: NDArray[np.bool_]) -> Mesh:
    """``mesh`` with two element sets, ``left`` and ``right``, instead of one."""
    return Mesh(
        nodes=mesh.nodes,
        elements=mesh.elements,
        element_sets={"left": np.flatnonzero(left), "right": np.flatnonzero(~left)},
        node_sets=mesh.node_sets,
    )


def uniform(sxx: float) -> Field:
    def field(x: NDArray[np.float64]) -> NDArray[np.float64]:
        out = np.zeros((*x.shape[:-1], 6))
        out[..., 0] = sxx
        return out

    return field


def solution(mesh: Mesh, field: Field | None) -> SolveResult:
    """Zero nodal results and, unless ``field`` is None, its values at the integration points."""
    points = None
    if field is not None:
        corners = mesh.nodes[mesh.elements[:, :4]]  # (M, 4, 3)
        at = np.einsum("qk,mkd->mqd", RULE, corners)
        points = PointStresses(points=at, stresses=field(at))
    n = mesh.n_nodes
    return SolveResult(
        displacements=np.zeros((n, 3)),
        stresses=np.zeros((n, 6)),
        solver="synthetic",
        point_stresses=points,
    )


def job(
    mesh: Mesh,
    materials: tuple[Material, ...] = (TO,),
    sections: tuple[Section, ...] | None = None,
    supports: tuple[FixedSupport, ...] = (),
) -> SolveJob:
    return SolveJob(
        mesh=mesh,
        materials=materials,
        sections=sections or (Section("part:p", materials[0].name),),
        supports=supports,
    )


def test_uniform_stress_gives_the_closed_form_ratio_everywhere() -> None:
    mesh = grid((2, 2, 2))
    result = summarize(job(mesh), solution(mesh, uniform(10.0)))
    strength = result.strength
    assert strength is not None
    assert strength.criteria == {"to": FailureCriterion.TSAI_WU}
    assert strength.lowest.ratio == pytest.approx(4.0, rel=1e-12)
    assert strength.lowest.mode is FailureMode.ALONG_RASTER_TENSION
    assert strength.clear_lowest is not None
    assert strength.clear_lowest.ratio == pytest.approx(4.0, rel=1e-12)
    assert sorted(s.face or 0 for s in strength.face_lowest) == [1, 2, 3, 4, 5, 6]
    assert all(s.ratio == pytest.approx(4.0, rel=1e-12) for s in strength.face_lowest)
    assert strength.assumed == {} and strength.unchecked == ()


def test_no_point_stresses_means_no_strength_result() -> None:
    mesh = grid((1, 1, 1))
    assert summarize(job(mesh), solution(mesh, None)).strength is None


def test_no_strength_data_means_no_strength_result() -> None:
    mesh = grid((1, 1, 1))
    plain = IsotropicMaterial("steel", 210000.0, 0.3)
    result = summarize(job(mesh, (plain,)), solution(mesh, uniform(10.0)))
    assert result.strength is None


def test_the_frame_of_each_section_is_used_and_the_lowest_element_wins_at_a_node() -> None:
    """Left half printed upright (R = 4 along the raster), right half on its side (R = 1.5
    across the layers). Nodes on the interface belong to both and take the lower ratio."""
    base = grid((2, 1, 1))
    centroid_x = base.nodes[base.elements[:, :4]].mean(axis=1)[:, 0]
    mesh = split(base, centroid_x < 1.0)
    # The interface x = 1 as a pretend face 7, to read the ratio its nodes get.
    interface = np.flatnonzero(np.isclose(mesh.nodes[:, 0], 1.0))
    mesh = Mesh(
        nodes=mesh.nodes,
        elements=mesh.elements,
        element_sets=mesh.element_sets,
        node_sets={**mesh.node_sets, "face:7": interface},
    )
    sections = (Section("left", "to"), Section("right", "to", SIDE))
    result = summarize(job(mesh, sections=sections), solution(mesh, uniform(10.0)))
    strength = result.strength
    assert strength is not None
    assert strength.lowest.ratio == pytest.approx(1.5, rel=1e-12)
    assert strength.lowest.mode is FailureMode.ACROSS_LAYERS_TENSION
    assert strength.lowest.element in set(mesh.element_sets["right"].tolist())
    by_face = {s.face: s for s in strength.face_lowest}
    assert by_face[1].ratio == pytest.approx(4.0, rel=1e-12)  # x = 0: upright only
    assert by_face[2].ratio == pytest.approx(1.5, rel=1e-12)  # x = 2: on its side
    assert by_face[7].ratio == pytest.approx(1.5, rel=1e-12)  # the interface takes the lower
    assert [s.ratio for s in strength.face_lowest] == sorted(s.ratio for s in strength.face_lowest)


def test_materials_without_strength_are_listed_and_skipped() -> None:
    base = grid((2, 1, 1))
    centroid_x = base.nodes[base.elements[:, :4]].mean(axis=1)[:, 0]
    mesh = split(base, centroid_x < 1.0)
    weak = IsotropicMaterial("plain", 1000.0, 0.3)
    sections = (Section("left", "plain"), Section("right", "to"))
    result = summarize(job(mesh, (weak, TO), sections), solution(mesh, uniform(10.0)))
    strength = result.strength
    assert strength is not None
    assert strength.unchecked == ("plain",)
    assert strength.criteria == {"to": FailureCriterion.TSAI_WU}
    left_only = np.flatnonzero(mesh.nodes[:, 0] < 1.0 - 1e-9)
    assert strength.lowest.node not in set(left_only.tolist())


def test_assumed_strengths_are_listed_by_material() -> None:
    mesh = grid((1, 1, 1))
    tension_only = OrthotropicStrength(t1=40.0, t2=30.0, t3=15.0, s12=20.0, s13=12.0, s23=12.0)
    material = OrthotropicMaterial("pla", **STIFFNESS, strength=tension_only)
    result = summarize(job(mesh, (material,)), solution(mesh, uniform(-10.0)))
    assert result.strength is not None
    assert result.strength.assumed == {"pla": ("c1", "c2", "c3")}
    assert result.strength.lowest.ratio == pytest.approx(4.0, rel=1e-12)  # C1 = T1 = 40
    assert result.strength.lowest.mode is FailureMode.ALONG_RASTER_COMPRESSION


def test_isotropic_materials_use_their_criteria() -> None:
    mesh = grid((1, 1, 1))
    steel = IsotropicMaterial("steel", 210000.0, 0.3, strength=IsotropicStrength(250.0))
    result = summarize(job(mesh, (steel,)), solution(mesh, uniform(-100.0)))
    assert result.strength is not None
    assert result.strength.criteria == {"steel": FailureCriterion.VON_MISES}
    assert result.strength.lowest.ratio == pytest.approx(2.5, rel=1e-12)
    assert result.strength.lowest.mode is FailureMode.VON_MISES


def test_a_flagged_lowest_ratio_never_becomes_the_clear_one() -> None:
    """Stress rises towards a clamped face at x = 0, so the lowest ratio sits next to it."""
    mesh = grid((6, 2, 2))
    clamp = FixedSupport("support:root")
    mesh = Mesh(
        nodes=mesh.nodes,
        elements=mesh.elements,
        element_sets=mesh.element_sets,
        node_sets={**mesh.node_sets, "support:root": mesh.node_sets["face:1"]},
    )

    def rising(x: NDArray[np.float64]) -> NDArray[np.float64]:
        out = np.zeros((*x.shape[:-1], 6))
        out[..., 0] = 20.0 - 2.0 * x[..., 0]
        return out

    result = summarize(job(mesh, supports=(clamp,)), solution(mesh, rising))
    strength = result.strength
    assert strength is not None
    assert strength.lowest.flagged
    assert strength.lowest.ratio == pytest.approx(2.0, rel=1e-9)  # 40 / 20 at x = 0
    assert strength.clear_lowest is not None and not strength.clear_lowest.flagged
    assert strength.clear_lowest.ratio > strength.lowest.ratio


def test_mesh_check_judges_the_lowest_ratio() -> None:
    coarse_mesh, fine_mesh = grid((2, 2, 2), 1.5), grid((3, 3, 3), 1.0)
    coarse = summarize(job(coarse_mesh), solution(coarse_mesh, uniform(10.0)))
    fine = summarize(job(fine_mesh), solution(fine_mesh, uniform(10.0)))
    check = compare_meshes(coarse, fine).lowest_strength_ratio
    assert check is not None
    assert check.coarse == pytest.approx(4.0) and check.fine == pytest.approx(4.0)
    assert check.verdict is Verdict.SETTLED


def test_mesh_check_has_no_strength_verdict_without_strengths_or_stress() -> None:
    coarse_mesh, fine_mesh = grid((2, 2, 2), 1.5), grid((3, 3, 3), 1.0)
    plain = (IsotropicMaterial("m", 1000.0, 0.3),)
    coarse = summarize(job(coarse_mesh, plain), solution(coarse_mesh, uniform(10.0)))
    fine = summarize(job(fine_mesh, plain), solution(fine_mesh, uniform(10.0)))
    assert compare_meshes(coarse, fine).lowest_strength_ratio is None
    coarse = summarize(job(coarse_mesh), solution(coarse_mesh, uniform(0.0)))
    fine = summarize(job(fine_mesh), solution(fine_mesh, uniform(0.0)))
    assert fine.strength is not None and np.isinf(fine.strength.lowest.ratio)
    assert compare_meshes(coarse, fine).lowest_strength_ratio is None
