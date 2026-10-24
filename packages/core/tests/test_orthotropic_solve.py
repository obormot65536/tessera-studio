"""An orthotropic part solved end to end through tessera.core (needs ``ccx``), ADR-009 PR-6.

A smoke test of the wiring, not the TO-01 benchmark: a 10 mm cube with the catalog's TO
constants, rollers on the three coordinate planes through the origin and 10 MPa of tension on
the x = L face. The exact field is linear, u = eps x, so the quadratic tetrahedra reproduce it
up to ``.frd`` output precision. Printed upright (build Z), x is material axis 1; printed on its
side (build (1, 0, 0)), x is material axis 3 and y stays axis 2.
"""

from pathlib import Path

import gmsh
import numpy as np
import pytest

from tessera.core import describe_faces, edit, new_project, solve
from tessera.model import OrthotropicDefinition

L, SIGMA = 10.0, 10.0
E1, E2, E3 = 3000.0, 2400.0, 1800.0
NU12, NU13, NU23 = 0.35, 0.30, 0.30
TO = OrthotropicDefinition(
    id="to-catalog", e1=E1, e2=E2, e3=E3, nu12=NU12, nu13=NU13, nu23=NU23,
    g12=1000.0, g13=700.0, g23=650.0,
)  # fmt: skip

# Strains (x, y, z) under SIGMA along x, from the major ratios and nu_ji = nu_ij E_j / E_i.
UPRIGHT = (SIGMA / E1, -NU12 * SIGMA / E1, -NU13 * SIGMA / E1)  # x = 1, y = 2, z = 3
ON_SIDE = (SIGMA / E3, -NU23 * SIGMA / E2, -NU13 * SIGMA / E1)  # x = 3, y = 2, z = 1


@pytest.fixture(scope="module")
def cube_step(tmp_path_factory: pytest.TempPathFactory) -> Path:
    path = tmp_path_factory.mktemp("geom") / "cube.step"
    gmsh.initialize(interruptible=False)
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.occ.addBox(0, 0, 0, L, L, L)
        gmsh.model.occ.synchronize()
        gmsh.write(str(path))
    finally:
        gmsh.finalize()
    return path


@pytest.mark.solver
@pytest.mark.parametrize(
    ("build", "strain"),
    [((0.0, 0.0, 1.0), UPRIGHT), ((1.0, 0.0, 0.0), ON_SIDE)],
    ids=["build-z", "build-x"],
)
def test_uniaxial_cube_matches_closed_form(
    cube_step: Path, tmp_path: Path, build: tuple[float, float, float], strain: tuple[float, ...]
) -> None:
    project = new_project("cube")
    edit.add_material(project, TO)
    edit.add_part(project, cube_step, part_id="cube", material_id=TO.id)
    edit.set_print_setup(project, "cube", "fdm", build, 0.0)
    faces = describe_faces(project, "cube")

    def face(normal: tuple[float, float, float]) -> int:
        return next(f.index for f in faces if f.normal == pytest.approx(normal))

    # Rollers: u_x = 0 on x = 0, u_y = 0 on y = 0, u_z = 0 on z = 0.
    for dof, normal in enumerate([(-1.0, 0.0, 0.0), (0.0, -1.0, 0.0), (0.0, 0.0, -1.0)], 1):
        edit.add_support(project, "cube", [face(normal)], dofs=(dof,))
    edit.add_force_load(project, "cube", [face((1.0, 0.0, 0.0))], (SIGMA * L * L, 0.0, 0.0))
    edit.set_mesh_size(project, 5.0)

    result = solve(project, workdir=tmp_path)

    assert result.orthotropic
    nodes, disp = result.mesh.nodes, result.solution.displacements
    far = np.isclose(nodes[:, 0], L)
    assert far.any()
    np.testing.assert_allclose(disp[far, 0], strain[0] * L, rtol=1e-5)
    np.testing.assert_allclose(disp, nodes * np.asarray(strain), atol=1e-6)
