"""CalculiX adapter end-to-end tests (need ``ccx``).

Patch test: a rectangular block of C3D10 elements under uniaxial tension. The exact solution is
a linear displacement field, which quadratic tetrahedra reproduce exactly, so every node must
match it up to solver round-off and ``.frd`` output precision.

Integration-point stresses (ADR-008) are checked on fields the elements reproduce exactly: the
patch test, a rotated uniaxial stress with all six components distinct, and pure bending (a
linear stress field). Where every point of an element is exact, so is the affine field through
them that ``tessera.post`` evaluates at the element's nodes.
"""

import dataclasses
import shutil
import sys
from collections.abc import Callable
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
from tessera.solvers import CalculixAdapter, calculix
from tessera.solvers.calculix import CCX_ENV_VAR, parse_dat

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


# --------------------------------------------------------------------------------------------
# Integration-point stresses (ADR-008)
# --------------------------------------------------------------------------------------------

POINT_REL_TOL = 1e-6
"""Integration-point stresses are exact for these fields; job.dat prints 7 significant digits."""

# C3D10 4-point rule: point k has barycentric coordinate B at corner k and A at the other three.
A = (5.0 - np.sqrt(5.0)) / 20.0
B = (5.0 + 3.0 * np.sqrt(5.0)) / 20.0
RULE_BARYCENTRIC = np.full((4, 4), A) + np.eye(4) * (B - A)

# Consistent TRI6 "mass" matrix / area: f = area * TRI6_MASS @ t integrates N_i t exactly for a
# traction t interpolated at the six nodes (corners, then midsides of edges 01, 12, 20).
TRI6_MASS = (
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

StressField = Callable[[NDArray[np.float64]], NDArray[np.float64]]
"""Points (P, 3) -> stress tensors (P, 3, 3) in MPa."""


def test_point_rule_constants() -> None:
    assert pytest.approx(0.1381966, abs=1e-7) == A
    assert pytest.approx(0.5854102, abs=1e-7) == B
    assert pytest.approx(1 - 3 * A) == B
    assert TRI6_MASS.sum() == pytest.approx(1.0)


def tet10_map(nodes: NDArray[np.float64], xi: NDArray[np.float64]) -> NDArray[np.float64]:
    """Isoparametric tet10 map of natural coordinates ``xi`` (P, 3) for one element's nodes."""
    lam = np.column_stack([1.0 - xi.sum(axis=1), xi])  # (P, 4)
    shape = [lam[:, k] * (2 * lam[:, k] - 1) for k in range(4)]
    shape += [4 * lam[:, i] * lam[:, j] for i, j in TET10_EDGES]
    return np.column_stack(shape) @ nodes


def traction_job(
    field: StressField, size: tuple[float, float, float], cells: tuple[int, int, int]
) -> SolveJob:
    """Box of C3D10 elements loaded on all six faces by the traction of a linear stress field.

    The nodal forces are consistent (``TRI6_MASS`` applied to t = sigma . n at the face's six
    nodes), which is exact for a linear field. A 3-2-1 support removes the rigid-body motion
    without restraining the self-balanced solution.
    """
    nodes, elements = block_tet10(*size, *cells)
    extent = np.array(size)
    tol = 1e-9 * extent.max()
    edge_index = {frozenset(e): 4 + k for k, e in enumerate(TET10_EDGES)}
    faces = ((0, 1, 2), (0, 1, 3), (1, 2, 3), (0, 2, 3))
    load_nodes: list[int] = []
    load_forces: list[NDArray[np.float64]] = []
    for conn in elements:
        for face in faces:
            corners = nodes[conn[list(face)]]
            for axis in range(3):
                for side, value in ((-1.0, 0.0), (1.0, extent[axis])):
                    if not np.all(np.abs(corners[:, axis] - value) < tol):
                        continue
                    ring = ((face[0], face[1]), (face[1], face[2]), (face[2], face[0]))
                    tri6 = [int(conn[c]) for c in face]
                    tri6 += [int(conn[edge_index[frozenset(e)]]) for e in ring]
                    p0, p1, p2 = corners
                    area = 0.5 * float(np.linalg.norm(np.cross(p1 - p0, p2 - p0)))
                    normal = np.zeros(3)
                    normal[axis] = side
                    traction = field(nodes[tri6]) @ normal  # (6, 3)
                    load_nodes.extend(tri6)
                    load_forces.append(area * TRI6_MASS @ traction)
    forces = np.vstack(load_forces)
    assert np.abs(forces.sum(axis=0)).max() < 1e-9 * np.abs(forces).max()  # self-balanced

    def nearest(point: tuple[float, float, float]) -> NDArray[np.int64]:
        return np.array([np.argmin(np.linalg.norm(nodes - np.array(point), axis=1))])

    mesh = Mesh(
        nodes=nodes,
        elements=elements,
        element_sets={"block": np.arange(len(elements), dtype=np.int64)},
        node_sets={
            "a": nearest((0, 0, 0)),
            "b": nearest((size[0], 0, 0)),
            "c": nearest((0, size[1], 0)),
        },
    )
    return SolveJob(
        mesh=mesh,
        materials=(IsotropicMaterial("pla", E, NU),),
        sections=(Section("block", "pla"),),
        supports=(FixedSupport("a"), FixedSupport("b", (2, 3)), FixedSupport("c", (3,))),
        loads=(NodalLoads(np.array(load_nodes), forces),),
        point_stresses=True,
    )


def straight_rule_points(mesh: Mesh) -> NDArray[np.float64]:
    """(M, 4, 3) rule points of straight-sided elements: barycentric weights of the corners."""
    return np.einsum("qk,mkd->mqd", RULE_BARYCENTRIC, mesh.nodes[mesh.elements[:, :4]])


@pytest.mark.solver
def test_patch_points_are_exact_without_the_frd_bias(tmp_path: Path) -> None:
    job = dataclasses.replace(tension_job(), point_stresses=True)
    result = CalculixAdapter(timeout=120).run(job, tmp_path)
    ps = result.point_stresses
    assert ps is not None
    m = job.mesh.n_elements
    assert ps.stresses.shape == (m, 4, 6)
    assert ps.points.shape == (m, 4, 3)
    assert result.metadata["point_rule"] == "tet4"

    exact = np.zeros((m, 4, 6))
    exact[..., 0] = SIGMA
    np.testing.assert_allclose(ps.stresses, exact, rtol=0, atol=POINT_REL_TOL * SIGMA)
    # The nodal .frd values of the same run carry ccx's extrapolation bias (about -1e-5, ADR-004
    # item 7); the integration points do not.
    frd_bias = float(np.abs(result.stresses[:, 0] - SIGMA).max())
    assert frd_bias > 5 * POINT_REL_TOL * SIGMA
    np.testing.assert_allclose(ps.points, straight_rule_points(job.mesh), rtol=0, atol=1e-6 * L)


@pytest.mark.solver
def test_rotated_uniaxial_tension_keeps_all_six_components_apart(tmp_path: Path) -> None:
    # sigma = 14 n (x) n with n = (1, 2, 3) / sqrt 14: xx yy zz xy xz yz = 1 4 9 2 3 6 MPa, so a
    # swap of any two components shows, in particular of xz and yz, which job.dat and job.frd
    # print in different orders.
    tensor = np.outer([1.0, 2.0, 3.0], [1.0, 2.0, 3.0])

    def field(points: NDArray[np.float64]) -> NDArray[np.float64]:
        return np.broadcast_to(tensor, (len(points), 3, 3))

    job = traction_job(field, (4.0, 3.0, 2.0), (2, 2, 1))
    result = CalculixAdapter(timeout=120).run(job, tmp_path)
    ps = result.point_stresses
    assert ps is not None
    neutral = np.array([1.0, 4.0, 9.0, 2.0, 6.0, 3.0])  # xx yy zz xy yz zx
    np.testing.assert_allclose(
        ps.stresses, np.broadcast_to(neutral, ps.stresses.shape), rtol=0, atol=POINT_REL_TOL * 9
    )
    # The nodal .frd reader agrees on the order, up to its bias.
    np.testing.assert_allclose(
        result.stresses,
        np.broadcast_to(neutral, result.stresses.shape),
        rtol=0,
        atol=STRESS_REL_TOL * 9,
    )


@pytest.mark.solver
def test_pure_bending_points_are_exact_and_below_the_surface_peak(tmp_path: Path) -> None:
    # sigma_xx = k (z - H / 2): a linear stress field, which quadratic tetrahedra reproduce.
    k = 2.0
    peak = k * H / 2

    def field(points: NDArray[np.float64]) -> NDArray[np.float64]:
        out = np.zeros((len(points), 3, 3))
        out[:, 0, 0] = k * (points[:, 2] - H / 2)
        return out

    job = traction_job(field, (L, W, H), (NX, NY, NZ))
    result = CalculixAdapter(timeout=120).run(job, tmp_path)
    ps = result.point_stresses
    assert ps is not None
    exact = np.zeros(ps.stresses.shape)
    exact[..., 0] = k * (ps.points[..., 2] - H / 2)
    np.testing.assert_allclose(ps.stresses, exact, rtol=0, atol=POINT_REL_TOL * peak)

    # The points lie inside the elements, so their maximum stays below the surface peak that
    # the nodal values reach: why reported peaks stay nodal (ADR-008, decision 4).
    nodal_peak = float(np.abs(result.stresses[:, 0]).max())
    assert nodal_peak == pytest.approx(peak, rel=STRESS_REL_TOL)
    point_peak = float(np.abs(ps.stresses[..., 0]).max())
    assert point_peak < 0.9 * nodal_peak


def _edit_deck(monkeypatch: pytest.MonkeyPatch, edit: Callable[[str], str]) -> None:
    write = calculix.write_deck
    monkeypatch.setattr(calculix, "write_deck", lambda job: edit(write(job)))


def _also_print_coord(deck: str) -> str:
    """Ask ccx for COORD too: the oracle for ``rule_points`` (production decks do not)."""
    return deck.replace("*END STEP\n", "*EL PRINT, ELSET=EALL\nCOORD\n*END STEP\n")


@pytest.mark.solver
def test_point_coordinates_follow_the_rule_on_a_curved_tet(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    corners = np.array([[0, 0, 0], [1.2, 0.1, 0], [0.2, 1.1, 0.1], [0.1, 0.2, 0.9]], dtype=float)
    mid = np.array([(corners[i] + corners[j]) / 2 for i, j in TET10_EDGES])
    bulge = np.array(
        [
            [0.0, -0.06, 0.04],
            [0.05, 0.05, 0.0],
            [-0.06, 0.0, 0.03],
            [-0.04, 0.03, 0.0],
            [0.06, 0.0, 0.05],
            [0.0, 0.06, 0.04],
        ]
    )
    nodes = np.vstack([corners, mid + bulge])
    mesh = Mesh(
        nodes=nodes,
        elements=np.arange(10, dtype=np.int64).reshape(1, 10),
        element_sets={"part": np.array([0])},
        node_sets={"base": np.array([0, 1, 2, 4, 5, 6])},
    )
    job = SolveJob(
        mesh=mesh,
        materials=(IsotropicMaterial("pla", E, NU),),
        sections=(Section("part", "pla"),),
        supports=(FixedSupport("base"),),
        loads=(NodalLoads(np.array([3]), np.array([[0.5, -0.25, 1.0]])),),
        point_stresses=True,
    )
    _edit_deck(monkeypatch, _also_print_coord)
    result = CalculixAdapter(timeout=120).run(job, tmp_path)
    ps = result.point_stresses
    assert ps is not None
    # ccx's own COORD output pins the rule and its point order (ADR-008, finding 4) ...
    coord = parse_dat((tmp_path / "job.dat").read_text()).blocks["global coordinates"]
    np.testing.assert_allclose(ps.points[0], coord.values[0], rtol=0, atol=1e-6)
    # ... and both are the isoparametric map of (a, a, a), (b, a, a), (a, b, a), (a, a, b).
    curved = tet10_map(nodes, RULE_BARYCENTRIC[:, 1:])
    np.testing.assert_allclose(ps.points[0], curved, rtol=0, atol=1e-12)
    straight = RULE_BARYCENTRIC @ corners
    assert np.abs(ps.points[0] - straight).max() > 1e-2  # the curved midsides are visible
    assert np.isfinite(ps.stresses).all()


PRINT_CARDS = "*EL PRINT, ELSET=EALL\nS\n"


def _run_with_edited_deck(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, edit: Callable[[str], str]
) -> None:
    """Run the patch job (point stresses on, plus a set ``half``) with ``edit`` on its deck."""
    job = tension_job()
    mesh = job.mesh
    half = {"half": np.arange(mesh.n_elements // 2, dtype=np.int64)}
    mesh = Mesh(mesh.nodes, mesh.elements, {**mesh.element_sets, **half}, mesh.node_sets)
    _edit_deck(monkeypatch, edit)
    job = dataclasses.replace(job, mesh=mesh, point_stresses=True)
    CalculixAdapter(timeout=120).run(job, tmp_path)


@pytest.mark.parametrize(
    ("cards", "message"),
    [
        ("", r"job\.dat is missing or empty"),  # no *EL PRINT: an empty job.dat (finding 7)
        ("*EL PRINT\nS\n", r"job\.dat is missing or empty"),  # no set: a warning, empty .dat
        ("*EL PRINT, ELSET=EALL\nCOORD\n", "'stresses' block"),  # S block missing
        ("*EL PRINT, ELSET=E_half\nS\n", r"no integration-point result for \d+ element"),
    ],
)
@pytest.mark.solver
def test_missing_or_incomplete_dat_is_a_solver_error(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, cards: str, message: str
) -> None:
    with pytest.raises(SolverError, match=message) as info:
        _run_with_edited_deck(tmp_path, monkeypatch, lambda deck: deck.replace(PRINT_CARDS, cards))
    assert info.value.log


@pytest.mark.solver
def test_duplicate_dat_rows_are_a_solver_error(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # A second definition of EALL makes ccx append to it and print element 1 twice, with no
    # warning (finding 8). Tessera never writes one; the parser is the second guard.
    def define_again(deck: str) -> str:
        return deck.replace("*STEP\n", "*ELSET, ELSET=EALL\n1\n*STEP\n")

    with pytest.raises(SolverError, match="duplicate") as info:
        _run_with_edited_deck(tmp_path, monkeypatch, define_again)
    assert info.value.log


# --- Off by default -----------------------------------------------------------------------


def _ccx_that_then_writes(tmp_path: Path, dat: str | None) -> str:
    """A wrapper around the real ccx that replaces job.dat afterwards (None: deletes it)."""
    ccx = shutil.which("ccx")
    assert ccx is not None
    after = "rm -f job.dat" if dat is None else f"printf '%s' '{dat}' > job.dat"
    script = tmp_path / "ccx-wrapper.sh"
    script.write_text(f'#!/bin/sh\n"{ccx}" "$@"\nstatus=$?\n{after}\nexit $status\n')
    script.chmod(0o755)
    return str(script)


@pytest.mark.parametrize("dat", [None, "", "garbage 1 2 3\n stresses (elem, integ.pnt.,x"])
@pytest.mark.solver
@pytest.mark.skipif(sys.platform == "win32", reason="POSIX shell wrapper")
def test_default_job_never_reads_job_dat(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, dat: str | None
) -> None:
    workdir = tmp_path / "run"
    monkeypatch.setenv(CCX_ENV_VAR, _ccx_that_then_writes(tmp_path, dat))
    result = CalculixAdapter(timeout=120).run(tension_job(), workdir)
    assert result.point_stresses is None
    assert "point_rule" not in result.metadata
    deck = (workdir / "job.inp").read_text()
    assert "*EL PRINT" not in deck and "EALL" not in deck
    assert (workdir / "job.dat").exists() == (dat is not None)  # the wrapper really ran
