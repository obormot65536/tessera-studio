"""CalculiX adapter end-to-end tests (need ``ccx``).

Patch test: a rectangular block of C3D10 elements under uniaxial tension. The exact solution is
a linear displacement field, which quadratic tetrahedra reproduce exactly, so every node must
match it up to solver round-off and ``.frd`` output precision.
"""

from itertools import permutations
from pathlib import Path

import numpy as np
import pytest
from numpy.typing import NDArray

from tessera.model.mesh import TET10_EDGES, Mesh
from tessera.model.solver import (
    FixedSupport,
    IsotropicMaterial,
    NodalLoads,
    Section,
    SolveJob,
    SolverError,
)
from tessera.solvers import CalculixAdapter

# Values chosen so every exact nodal result has a short decimal expansion: the ASCII .frd file
# stores six significant digits (E12.5), so a generic value would carry up to 5e-6 relative
# rounding error and could not be checked to 1e-6. With these numbers frd rounding is exact.
E, NU, SIGMA = 1000.0, 0.25, 2.0
L, W, H = 10.0, 4.0, 2.0
NX, NY, NZ = 4, 2, 2

REL_TOL = 1e-6
# Nodal stresses in the .frd are extrapolated by ccx from the integration points. The
# integration-point stresses are exact here (2.000000 in the .dat), but the frd nodal values come
# out as 1.99998 at every node: a systematic -1e-5 relative bias, consistent with ccx's C3D10
# extrapolation coefficients being stored to five decimals (1.92705 - 3 * 0.30902 = 0.99999).
STRESS_REL_TOL = 2e-5


def block_tet10(
    lx: float, ly: float, lz: float, nx: int, ny: int, nz: int
) -> tuple[NDArray[np.float64], NDArray[np.int64]]:
    """Tet10 mesh of a box: each hex split into 6 Kuhn tetrahedra (conforming across hexes)."""
    xs, ys, zs = np.linspace(0, lx, nx + 1), np.linspace(0, ly, ny + 1), np.linspace(0, lz, nz + 1)

    def vid(i: int, j: int, k: int) -> int:
        return (i * (ny + 1) + j) * (nz + 1) + k

    nodes: list[list[float]] = [[x, y, z] for x in xs for y in ys for z in zs]
    unit = np.eye(3, dtype=np.int64)
    tets: list[list[int]] = []
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                for perm in permutations(range(3)):
                    p = np.array([i, j, k])
                    corners = [vid(*p)]
                    for axis in perm:
                        p = p + unit[axis]
                        corners.append(vid(*p))
                    a, b, c, d = (np.array(nodes[v]) for v in corners)
                    if np.dot(np.cross(b - a, c - a), d - a) < 0:
                        corners[1], corners[2] = corners[2], corners[1]
                    tets.append(corners)

    midside: dict[tuple[int, int], int] = {}
    elements: list[list[int]] = []
    for corners in tets:
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


def tension_job() -> SolveJob:
    nodes, elements = block_tet10(L, W, H, NX, NY, NZ)
    tol = 1e-9 * L
    x, y, z = nodes.T

    # Faces on x = L: tet faces whose three corners all lie on x = L. Under uniform traction a
    # tri6 gets zero at its corners and A/3 * traction at each midside node.
    faces = ((0, 1, 2), (0, 1, 3), (1, 2, 3), (0, 2, 3))
    edge_index = {frozenset(e): 4 + n for n, e in enumerate(TET10_EDGES)}
    load_nodes: list[int] = []
    load_forces: list[list[float]] = []
    area_total = 0.0
    for conn in elements:
        for face in faces:
            if all(abs(x[conn[c]] - L) < tol for c in face):
                p0, p1, p2 = (nodes[conn[c]] for c in face)
                area = 0.5 * float(np.linalg.norm(np.cross(p1 - p0, p2 - p0)))
                area_total += area
                for a, b in ((face[0], face[1]), (face[1], face[2]), (face[2], face[0])):
                    load_nodes.append(int(conn[edge_index[frozenset((a, b))]]))
                    load_forces.append([SIGMA * area / 3.0, 0.0, 0.0])
    assert abs(area_total - W * H) < 1e-9 * W * H

    mesh = Mesh(
        nodes=nodes,
        elements=elements,
        element_sets={"block": np.arange(len(elements), dtype=np.int64)},
        node_sets={
            "x0": np.flatnonzero(np.abs(x) < tol),
            "y0": np.flatnonzero(np.abs(y) < tol),
            "z0": np.flatnonzero(np.abs(z) < tol),
        },
    )
    return SolveJob(
        mesh=mesh,
        materials=(IsotropicMaterial("pla", E, NU),),
        sections=(Section("block", "pla"),),
        supports=(FixedSupport("x0", (1,)), FixedSupport("y0", (2,)), FixedSupport("z0", (3,))),
        loads=(NodalLoads(np.array(load_nodes), np.array(load_forces)),),
    )


def test_block_mesh_is_valid() -> None:
    nodes, elements = block_tet10(L, W, H, NX, NY, NZ)
    a, b, c, d = (nodes[elements[:, n]] for n in range(4))
    vol = np.einsum("ij,ij->i", np.cross(b - a, c - a), d - a) / 6.0
    assert (vol > 0).all()
    assert vol.sum() == pytest.approx(L * W * H)
    for n, (e0, e1) in enumerate(TET10_EDGES):
        mid = (nodes[elements[:, e0]] + nodes[elements[:, e1]]) / 2
        np.testing.assert_allclose(nodes[elements[:, 4 + n]], mid)
    assert len(np.unique(nodes.round(9), axis=0)) == len(nodes)  # midsides are shared


@pytest.mark.solver
def test_uniaxial_tension_patch(tmp_path: Path) -> None:
    job = tension_job()
    result = CalculixAdapter(timeout=120).run(job, tmp_path)

    assert result.solver == "calculix"
    assert result.metadata.get("ccx_version", "").startswith("Version 2.")
    assert result.metadata["ccx_name.elset.block"] == "E_block"
    n = job.mesh.n_nodes
    assert result.displacements.shape == (n, 3)
    assert result.stresses.shape == (n, 6)

    x, y, z = job.mesh.nodes.T
    exact_u = np.column_stack([SIGMA * x / E, -NU * SIGMA * y / E, -NU * SIGMA * z / E])
    np.testing.assert_allclose(result.displacements, exact_u, rtol=0, atol=REL_TOL * SIGMA * L / E)

    on_end = np.abs(x - L) < 1e-9
    np.testing.assert_allclose(result.displacements[on_end, 0], SIGMA * L / E, rtol=REL_TOL)
    on_side = np.abs(y - W) < 1e-9
    np.testing.assert_allclose(result.displacements[on_side, 1], -NU * SIGMA * W / E, rtol=REL_TOL)
    on_top = np.abs(z - H) < 1e-9
    np.testing.assert_allclose(result.displacements[on_top, 2], -NU * SIGMA * H / E, rtol=REL_TOL)

    exact_s = np.zeros((n, 6))
    exact_s[:, 0] = SIGMA
    np.testing.assert_allclose(result.stresses, exact_s, rtol=0, atol=STRESS_REL_TOL * SIGMA)


@pytest.mark.solver
def test_inverted_element_raises(tmp_path: Path) -> None:
    job = tension_job()
    elements = job.mesh.elements.copy()
    # Swap corners 1 and 2 (and the matching midsides) of the first element: negative volume.
    elements[0, [1, 2]] = elements[0, [2, 1]]
    elements[0, [4, 6, 8, 9]] = elements[0, [6, 4, 9, 8]]
    mesh = Mesh(job.mesh.nodes, elements, job.mesh.element_sets, job.mesh.node_sets)
    broken = SolveJob(mesh, job.materials, job.sections, job.supports, job.loads)
    with pytest.raises(SolverError, match="nonpositive jacobian") as info:
        CalculixAdapter(timeout=120).run(broken, tmp_path)
    assert info.value.log  # the solver log travels with the error
