"""Orthotropic CalculiX runs (need ``ccx``): ADR-009 PR-2, the early answer before TO-01.

A structured C3D10 cube under uniaxial tension along material axis 1, with rollers on the three
coordinate planes through the origin. The exact field is linear (u = eps x), which quadratic
tetrahedra reproduce, so every node must match the closed-form compliance up to ``.frd`` output
precision. Run with (a) the identity frame and (b) the frame for a part printed on its side,
build direction (1, 0, 0): e1 = (0, 0, -1), e2 = (0, 1, 0), e3 = (1, 0, 0). Case (b) loads along
global z, so it fails if ccx swaps axes, ignores the orientation, or returns stresses in
material axes instead of global ones.
"""

from itertools import permutations
from pathlib import Path

import numpy as np
import pytest
from numpy.typing import NDArray

from tessera.model.elasticity import orthotropic_compliance, rotate_stress
from tessera.model.mesh import TET10_EDGES, Mesh
from tessera.model.solver import (
    FixedSupport,
    MaterialFrame,
    NodalLoads,
    OrthotropicMaterial,
    Section,
    SolveJob,
)
from tessera.solvers import CalculixAdapter

# ADR-009 catalog TO constants (made-up test values, not material data).
CONSTANTS = (3000.0, 2400.0, 1800.0, 0.35, 0.30, 0.30, 1000.0, 700.0, 650.0)
L, N, SIGMA = 10.0, 2, 10.0

# frd output precision, as measured in test_calculix_patch: E12.5 rounding (5e-6 relative) and
# the -1e-5 bias of ccx's C3D10 nodal stress extrapolation.
DISP_REL_TOL = 1e-5
STRESS_REL_TOL = 2e-5

IDENTITY = MaterialFrame((1.0, 0.0, 0.0), (0.0, 1.0, 0.0))
SIDE = MaterialFrame((0.0, 0.0, -1.0), (0.0, 1.0, 0.0))  # build direction (1, 0, 0)


def cube_tet10(length: float, n: int) -> tuple[NDArray[np.float64], NDArray[np.int64]]:
    """Tet10 mesh of a cube: each hex split into 6 Kuhn tetrahedra (as in test_calculix_patch)."""
    xs = np.linspace(0.0, length, n + 1)

    def vid(i: int, j: int, k: int) -> int:
        return (i * (n + 1) + j) * (n + 1) + k

    nodes: list[list[float]] = [[x, y, z] for x in xs for y in xs for z in xs]
    unit = np.eye(3, dtype=np.int64)
    midside: dict[tuple[int, int], int] = {}
    elements: list[list[int]] = []
    for i in range(n):
        for j in range(n):
            for k in range(n):
                for perm in permutations(range(3)):
                    p = np.array([i, j, k])
                    corners = [vid(*p)]
                    for axis in perm:
                        p = p + unit[axis]
                        corners.append(vid(*p))
                    a, b, c, d = (np.array(nodes[v]) for v in corners)
                    if np.dot(np.cross(b - a, c - a), d - a) < 0:
                        corners[1], corners[2] = corners[2], corners[1]
                    conn = list(corners)
                    for e0, e1 in TET10_EDGES:
                        u, v = corners[e0], corners[e1]
                        key = (min(u, v), max(u, v))
                        if key not in midside:
                            midside[key] = len(nodes)
                            nodes.append(list((np.array(nodes[u]) + np.array(nodes[v])) / 2))
                        conn.append(midside[key])
                    elements.append(conn)
    return np.array(nodes, dtype=np.float64), np.array(elements, dtype=np.int64)


def tension_job(frame: MaterialFrame) -> tuple[SolveJob, int]:
    """Uniaxial SIGMA along material axis 1, applied on the cube face it points to (+/-)."""
    nodes, elements = cube_tet10(L, N)
    axis = int(np.argmax(np.abs(frame.axis_1)))
    tol = 1e-9 * L
    coord = nodes[:, axis]
    # Uniform traction on a tri6: zero at the corners, A/3 at each midside node.
    faces = ((0, 1, 2), (0, 1, 3), (1, 2, 3), (0, 2, 3))
    edge_index = {frozenset(e): 4 + k for k, e in enumerate(TET10_EDGES)}
    load_nodes: list[int] = []
    load_forces: list[NDArray[np.float64]] = []
    area_total = 0.0
    for conn in elements:
        for face in faces:
            if all(abs(coord[conn[c]] - L) < tol for c in face):
                p0, p1, p2 = (nodes[conn[c]] for c in face)
                area = 0.5 * float(np.linalg.norm(np.cross(p1 - p0, p2 - p0)))
                area_total += area
                for a, b in ((face[0], face[1]), (face[1], face[2]), (face[2], face[0])):
                    load_nodes.append(int(conn[edge_index[frozenset((a, b))]]))
                    load_forces.append(SIGMA * area / 3.0 * np.eye(3)[axis])
    assert abs(area_total - L * L) < 1e-9 * L * L

    x, y, z = nodes.T
    mesh = Mesh(
        nodes=nodes,
        elements=elements,
        element_sets={"cube": np.arange(len(elements), dtype=np.int64)},
        node_sets={
            "x0": np.flatnonzero(np.abs(x) < tol),
            "y0": np.flatnonzero(np.abs(y) < tol),
            "z0": np.flatnonzero(np.abs(z) < tol),
        },
    )
    job = SolveJob(
        mesh=mesh,
        materials=(OrthotropicMaterial("to", *CONSTANTS),),
        sections=(Section("cube", "to", frame),),
        supports=(FixedSupport("x0", (1,)), FixedSupport("y0", (2,)), FixedSupport("z0", (3,))),
        loads=(NodalLoads(np.array(load_nodes), np.array(load_forces)),),
    )
    return job, axis


def exact_strain(frame: MaterialFrame, axis: int) -> NDArray[np.float64]:
    """Global strain tensor under SIGMA along global ``axis``, from the closed-form compliance."""
    sigma = np.zeros((3, 3))
    sigma[axis, axis] = SIGMA
    r = frame.rotation
    s_local = r @ sigma @ r.T
    voigt = np.array(
        [s_local[0, 0], s_local[1, 1], s_local[2, 2], s_local[0, 1], s_local[1, 2], s_local[2, 0]]
    )
    e = orthotropic_compliance(*CONSTANTS) @ voigt  # engineering shear strains
    eps_local = np.array(
        [[e[0], e[3] / 2, e[5] / 2], [e[3] / 2, e[1], e[4] / 2], [e[5] / 2, e[4] / 2, e[2]]]
    )
    return r.T @ eps_local @ r


@pytest.mark.solver
@pytest.mark.parametrize(
    ("frame", "expected_diag"),
    [
        # eps (x, y, z): sigma/E1, -nu12 sigma/E1, -nu13 sigma/E1
        (IDENTITY, (10 / 3000, -3.5 / 3000, -3.0 / 3000)),
        # e1 = -z, e2 = y, e3 = x: eps_x = -nu13 sigma/E1, eps_y = -nu12 sigma/E1, eps_z = sigma/E1
        (SIDE, (-3.0 / 3000, -3.5 / 3000, 10 / 3000)),
    ],
    ids=["identity", "build_x"],
)
def test_uniaxial_along_material_axis_1(
    tmp_path: Path, frame: MaterialFrame, expected_diag: tuple[float, float, float]
) -> None:
    job, axis = tension_job(frame)
    eps = exact_strain(frame, axis)
    np.testing.assert_allclose(eps, np.diag(expected_diag), rtol=0, atol=1e-15)

    result = CalculixAdapter(timeout=120).run(job, tmp_path)
    assert result.metadata["ccx_name.orientation.cube"] == "O_cube"

    nodes = job.mesh.nodes
    # Rollers on the planes through the origin and a diagonal strain: u = eps x exactly.
    exact_u = nodes @ eps.T
    u_max = SIGMA / CONSTANTS[0] * L
    np.testing.assert_allclose(result.displacements, exact_u, rtol=0, atol=DISP_REL_TOL * u_max)
    far = np.abs(nodes[:, axis] - L) < 1e-9
    np.testing.assert_allclose(result.displacements[far, axis], u_max, rtol=DISP_REL_TOL)

    # Strain = least-squares gradient of the nodal displacements.
    grad, *_ = np.linalg.lstsq(np.column_stack([nodes, np.ones(len(nodes))]), result.displacements)
    measured = grad[:3].T
    np.testing.assert_allclose(
        (measured + measured.T) / 2, eps, rtol=0, atol=DISP_REL_TOL * SIGMA / CONSTANTS[0]
    )

    # Stresses come back in GLOBAL axes: sigma along the loaded global axis, the rest zero.
    exact_s = np.zeros((job.mesh.n_nodes, 6))
    exact_s[:, axis] = SIGMA
    np.testing.assert_allclose(result.stresses, exact_s, rtol=0, atol=STRESS_REL_TOL * SIGMA)
    # ... and rotated into the material frame they are sigma_11 = SIGMA.
    local = rotate_stress(result.stresses, frame)
    exact_local = np.zeros_like(local)
    exact_local[:, 0] = SIGMA
    np.testing.assert_allclose(local, exact_local, rtol=0, atol=STRESS_REL_TOL * SIGMA)
