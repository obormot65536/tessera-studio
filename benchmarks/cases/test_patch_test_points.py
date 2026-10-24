"""Integration-point stresses end to end (ADR-008): TS-03 and pure bending through real ccx.

Problems (mm, N, MPa per ADR-002)
---------------------------------
1. **TS-03 at integration points.** The TS-03 cube (``test_patch_test.py``): steel, 10 mm,
   rollers on x = 0, y = 0, z = 0, a uniform 10 MPa traction on x = 10, unstructured gmsh
   C3D10 mesh, solved through ``tessera.core.solve``. The solver adapter is wrapped so the
   ``SolveJob`` asks for point stresses (``point_stresses=True``); core is unchanged.
   Exact: sigma_xx = 10 MPa, every other component 0, everywhere.
2. **Pure bending.** A 10 x 4 x 2 mm block, its corner at (1000, -500, 250) mm so the elements
   sit far from the origin, loaded on all six faces by the consistent nodal forces of
   sigma_xx = k (z - z_mid), k = 2 MPa/mm (peak 2 MPa), with a 3-2-1 support. Quadratic
   tetrahedra reproduce the quadratic displacement field exactly, so the stress, linear in z,
   is exact at every integration point and so is the affine field through them.

Measured quantities
-------------------
At every integration point (``SolveResult.point_stresses``) and at every element's ten nodes
(``tessera.post.element_nodal_stresses``, unaveraged): the stress against the exact field,
relative to the peak stress. For comparison, the nodal ``.frd`` stress of the same run.

Precision
---------
``job.dat`` prints stresses as ``E13.6`` (seven significant digits): each value is off by up to
5e-7 of itself. The point coordinates are computed from the mesh in double precision (not
read from ccx's COORD), so the fit adds only its own weights: at most 1.927 + 3 x 0.309 = 2.854
at a corner node, so element-nodal values are within 2.854 x 5e-7 = 1.43e-6 of the element's
largest stress (``E13_6_FIT_BOUND``). TS-03's 10.000000 prints exactly, so there the floor is
solver round-off. The nodal ``.frd`` values carry their -1e-5 extrapolation bias (TS-03).

Measured (CalculiX 2.21, SPOOLES):

    case                     elements   points    element-nodal   .frd nodal   (max |error| / peak)
    TS-03, h = 3.0               1110   4.1e-13     7.0e-13         1.0e-05
    bending, 1000 mm away          96   2.4e-07     4.4e-07         8.7e-06

With the points read from ccx's COORD instead (7 digits of a coordinate near 1000 mm), the
bending block measured 1.2e-5 at the points and 2.7e-5 at the element nodes, worse than the
.frd values: the reason the adapter computes the points from the mesh.

The gate compares against ``references/patch_test_points.toml`` (proposed in the PR, created
by the maintainer); until it exists the gate is skipped, but every run still checks the
output-precision bound above.
"""

from __future__ import annotations

import dataclasses
import itertools
import shutil
from collections.abc import Callable
from pathlib import Path
from typing import Any, Literal

import numpy as np
import pytest
from numpy.typing import NDArray

pytestmark = [
    pytest.mark.benchmark,
    pytest.mark.solver,
    pytest.mark.skipif(shutil.which("ccx") is None, reason="CalculiX (ccx) not on PATH"),
]

REFERENCE = "patch_test_points"

PRINT_HALF_ULP = 5e-7
"""Half a unit in the seventh significant digit: the rounding of an E13.6 value, relative."""
E13_6_FIT_BOUND = 2.854 * PRINT_HALF_ULP
"""Largest error of an element-nodal stress from E13.6 point stresses, relative to the
element's largest stress (the tet4 extrapolation weights sum to at most 2.854 in magnitude)."""

# TS-03 (as in test_patch_test.py)
SIDE = 10.0
TRACTION = 10.0
GATE_MESH_SIZE = 3.0

# Pure bending
CORNER = np.array([1000.0, -500.0, 250.0])
BLOCK = (10.0, 4.0, 2.0)
CELLS = (4, 2, 2)
CURVATURE_STRESS = 2.0  # k, MPa/mm
YOUNGS_MODULUS, POISSON_RATIO = 1000.0, 0.25

_TET10_EDGES = ((0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3))
_TRI6_MASS = (
    np.array(
        [
            [6, -1, -1, 0, -4, 0],
            [-1, 6, -1, 0, 0, -4],
            [-1, -1, 6, -4, 0, 0],
            [0, 0, -4, 32, 16, 16],
            [-4, 0, 0, 16, 32, 16],
            [0, -4, 0, 16, 16, 32],
        ],
        dtype=np.float64,
    )
    / 180.0
)
"""Consistent TRI6 matrix / area: f = area * M t is exact for a traction interpolated at the
six nodes (corners, then midsides of edges 01, 12, 20)."""


@dataclasses.dataclass(frozen=True)
class Errors:
    """Max |stress - exact| / peak: at the points, at the element nodes, at the .frd nodes."""

    points: float
    element_nodal: float
    frd_nodal: float
    n_elements: int


class _AskForPoints:
    """CalculiX, with every job asking for integration-point stresses."""

    name = "calculix"

    def __init__(self) -> None:
        from tessera.solvers import CalculixAdapter

        self._inner = CalculixAdapter(timeout=600)

    def capabilities(self) -> frozenset[Any]:
        return self._inner.capabilities()

    def run(self, job: Any, workdir: Path) -> Any:
        return self._inner.run(dataclasses.replace(job, point_stresses=True), workdir)


def _errors(
    mesh: Any, solution: Any, exact: Callable[[NDArray[np.float64]], NDArray[np.float64]]
) -> Errors:
    from tessera.post import element_nodal_stresses

    ps = solution.point_stresses
    assert ps is not None and solution.metadata["point_rule"] == "tet4"
    peak = float(np.abs(exact(ps.points)).max())
    nodal = element_nodal_stresses(mesh, ps)
    return Errors(
        points=float(np.abs(ps.stresses - exact(ps.points)).max() / peak),
        element_nodal=float(np.abs(nodal - exact(mesh.nodes[mesh.elements])).max() / peak),
        frd_nodal=float(np.abs(solution.stresses - exact(mesh.nodes)).max() / peak),
        n_elements=mesh.n_elements,
    )


# --- TS-03 -------------------------------------------------------------------------------


def _face_with_normal(faces: list[Any], normal: tuple[float, float, float]) -> int:
    matches = [
        f.index for f in faces if f.normal is not None and np.allclose(f.normal, normal, atol=1e-6)
    ]
    assert len(matches) == 1, f"expected one face with normal {normal}, got {matches}"
    return matches[0]


def uniaxial(points: NDArray[np.float64]) -> NDArray[np.float64]:
    out = np.zeros((*points.shape[:-1], 6))
    out[..., 0] = TRACTION
    return out


def _solve_ts03(step_path: Path, workdir: Path) -> Errors:
    from tessera.core import describe_faces, new_project, solve
    from tessera.model import FaceSet, ForceLoad, PartInstance, Support

    p = new_project("patch-test-points")
    p.assembly.parts.append(
        PartInstance(id="cube", name="cube", geometry=step_path, material_id="steel-generic")
    )
    faces = describe_faces(p, "cube")
    loaded = _face_with_normal(faces, (1.0, 0.0, 0.0))
    rollers: tuple[tuple[str, tuple[float, float, float], Literal[1, 2, 3]], ...] = (
        ("roller-x", (-1.0, 0.0, 0.0), 1),
        ("roller-y", (0.0, -1.0, 0.0), 2),
        ("roller-z", (0.0, 0.0, -1.0), 3),
    )
    for name, normal, dof in rollers:
        face = _face_with_normal(faces, normal)
        p.study.supports.append(
            Support(id=name, faces=FaceSet(part_id="cube", faces=(face,)), dofs=(dof,))
        )
    p.study.loads.append(
        ForceLoad(
            id="traction",
            faces=FaceSet(part_id="cube", faces=(loaded,)),
            force=(TRACTION * SIDE * SIDE, 0.0, 0.0),
        )
    )
    p.study.mesh_size = GATE_MESH_SIZE
    workdir.mkdir(parents=True, exist_ok=True)
    result = solve(p, workdir=workdir, adapter=_AskForPoints())
    return _errors(result.mesh, result.solution, uniaxial)


# --- Pure bending ------------------------------------------------------------------------


def bending(points: NDArray[np.float64]) -> NDArray[np.float64]:
    out = np.zeros((*points.shape[:-1], 6))
    out[..., 0] = CURVATURE_STRESS * (points[..., 2] - (CORNER[2] + BLOCK[2] / 2))
    return out


def _block_tet10() -> tuple[NDArray[np.float64], NDArray[np.int64]]:
    """Tet10 mesh of the block: each hex cell split into six Kuhn tetrahedra."""
    nx, ny, nz = CELLS
    axes = [np.linspace(0, BLOCK[k], CELLS[k] + 1) for k in range(3)]
    nodes: list[NDArray[np.float64]] = [
        CORNER + np.array([x, y, z]) for x in axes[0] for y in axes[1] for z in axes[2]
    ]

    def vid(i: int, j: int, k: int) -> int:
        return (i * (ny + 1) + j) * (nz + 1) + k

    midside: dict[tuple[int, int], int] = {}
    elements: list[list[int]] = []
    for i, j, k in itertools.product(range(nx), range(ny), range(nz)):
        for perm in itertools.permutations(range(3)):
            v = [i, j, k]
            corners = [vid(*v)]
            for axis in perm:
                v[axis] += 1
                corners.append(vid(*v))
            a, b, c, d = (nodes[n] for n in corners)
            if np.dot(np.cross(b - a, c - a), d - a) < 0:
                corners[1], corners[2] = corners[2], corners[1]
            conn = list(corners)
            for e0, e1 in _TET10_EDGES:
                key = (min(corners[e0], corners[e1]), max(corners[e0], corners[e1]))
                if key not in midside:
                    midside[key] = len(nodes)
                    nodes.append((nodes[key[0]] + nodes[key[1]]) / 2)
                conn.append(midside[key])
            elements.append(conn)
    return np.array(nodes), np.array(elements, dtype=np.int64)


def _bending_job() -> Any:
    from tessera.model import (
        FixedSupport,
        IsotropicMaterial,
        Mesh,
        NodalLoads,
        Section,
        SolveJob,
    )

    nodes, elements = _block_tet10()
    low, high = CORNER, CORNER + np.array(BLOCK)
    tol = 1e-9 * max(BLOCK)
    edge_index = {frozenset(e): 4 + n for n, e in enumerate(_TET10_EDGES)}
    load_nodes: list[int] = []
    load_forces: list[NDArray[np.float64]] = []
    for conn in elements:
        for face in ((0, 1, 2), (0, 1, 3), (1, 2, 3), (0, 2, 3)):
            corners = nodes[conn[list(face)]]
            for axis, side in itertools.product(range(3), (-1.0, 1.0)):
                plane = low[axis] if side < 0 else high[axis]
                if not np.all(np.abs(corners[:, axis] - plane) < tol):
                    continue
                ring = ((face[0], face[1]), (face[1], face[2]), (face[2], face[0]))
                tri6 = [int(conn[c]) for c in face] + [
                    int(conn[edge_index[frozenset(e)]]) for e in ring
                ]
                p0, p1, p2 = corners
                area = 0.5 * float(np.linalg.norm(np.cross(p1 - p0, p2 - p0)))
                normal = np.zeros(3)
                normal[axis] = side
                s = bending(nodes[tri6])  # (6, 6) xx yy zz xy yz zx
                tensor = np.zeros((6, 3, 3))
                tensor[:, 0, 0] = s[:, 0]
                load_nodes.extend(tri6)
                load_forces.append(area * _TRI6_MASS @ (tensor @ normal))
    forces = np.vstack(load_forces)

    def nearest(point: NDArray[np.float64]) -> NDArray[np.int64]:
        return np.array([np.argmin(np.linalg.norm(nodes - point, axis=1))])

    mesh = Mesh(
        nodes=nodes,
        elements=elements,
        element_sets={"block": np.arange(len(elements), dtype=np.int64)},
        node_sets={
            "a": nearest(low),
            "b": nearest(low + np.array([BLOCK[0], 0, 0])),
            "c": nearest(low + np.array([0, BLOCK[1], 0])),
        },
    )
    return SolveJob(
        mesh=mesh,
        materials=(IsotropicMaterial("m", YOUNGS_MODULUS, POISSON_RATIO),),
        sections=(Section("block", "m"),),
        supports=(FixedSupport("a"), FixedSupport("b", (2, 3)), FixedSupport("c", (3,))),
        loads=(NodalLoads(np.array(load_nodes), forces),),
        point_stresses=True,
    )


def _solve_bending(workdir: Path) -> Errors:
    from tessera.solvers import CalculixAdapter

    job = _bending_job()
    solution = CalculixAdapter(timeout=600).run(job, workdir)
    return _errors(job.mesh, solution, bending)


# --- Tests -------------------------------------------------------------------------------


@pytest.fixture(scope="module")
def ts03(box_step: Callable[..., Path], tmp_path_factory: pytest.TempPathFactory) -> Errors:
    step = box_step("patch_cube", (0.0, 0.0, 0.0), (SIDE, SIDE, SIDE))
    return _solve_ts03(step, tmp_path_factory.mktemp("patch-test-points"))


@pytest.fixture(scope="module")
def pure_bending(tmp_path_factory: pytest.TempPathFactory) -> Errors:
    return _solve_bending(tmp_path_factory.mktemp("bending-points"))


def test_bending_loads_are_self_balanced() -> None:
    job = _bending_job()
    forces = job.loads[0].forces
    nodes = job.mesh.nodes[job.loads[0].nodes]
    assert np.abs(forces.sum(axis=0)).max() < 1e-9 * np.abs(forces).max()
    moment = np.cross(nodes - CORNER, forces).sum(axis=0)
    assert np.abs(moment).max() < 1e-9 * np.abs(forces).max() * max(BLOCK)


def _report(name: str, errors: Errors) -> None:
    print(
        f"\n{name}: {errors.n_elements} elements, max |error| / peak: points "
        f"{errors.points:.2e}, element-nodal {errors.element_nodal:.2e}, "
        f".frd nodal {errors.frd_nodal:.2e}"
    )


@pytest.mark.parametrize("case", ["ts03", "pure_bending"])
def test_point_stresses_at_output_precision(case: str, request: pytest.FixtureRequest) -> None:
    errors: Errors = request.getfixturevalue(case)
    _report(case, errors)
    # Every point value is the exact stress, printed to seven digits.
    assert errors.points <= PRINT_HALF_ULP
    # The affine fit through them adds only its weights (see the module docstring).
    assert errors.element_nodal <= E13_6_FIT_BOUND
    # And both beat the extrapolated, averaged .frd nodal values (TS-03's -1e-5 bias).
    assert errors.element_nodal < errors.frd_nodal


@pytest.mark.parametrize("case", ["ts03", "pure_bending"])
def test_point_stresses_gate(
    case: str, request: pytest.FixtureRequest, load_reference: Callable[[str], dict[str, Any]]
) -> None:
    # Solve first: fail (not skip) if the solve breaks, before the reference check can skip.
    errors: Errors = request.getfixturevalue(case)
    ref = load_reference(REFERENCE)
    assert errors.points <= ref["rel_tolerance_point_stress"], f"{case}: {errors}"
    assert errors.element_nodal <= ref["rel_tolerance_element_nodal_stress"], f"{case}: {errors}"
