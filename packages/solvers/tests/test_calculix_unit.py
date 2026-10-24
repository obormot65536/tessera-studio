"""CalculiX adapter unit tests that do not need the ccx executable."""

import dataclasses
from pathlib import Path

import numpy as np
import pytest
from hypothesis import given
from hypothesis import strategies as st

from tessera.model.mesh import Mesh
from tessera.model.solver import (
    FixedSupport,
    IsotropicMaterial,
    MaterialFrame,
    NodalLoads,
    OrthotropicMaterial,
    PointStresses,
    Section,
    SolveJob,
    SolverAdapter,
    SolverError,
)
from tessera.solvers import CalculixAdapter, calculix
from tessera.solvers.calculix import (
    CCX_ENV_VAR,
    DatBlock,
    DatResults,
    deck_names,
    parse_dat,
    parse_frd,
    point_stresses,
    read_point_stresses,
    rule_points,
    write_deck,
)

# Unit right-angle tet10: corners, then midsides on (0,1) (1,2) (2,0) (0,3) (1,3) (2,3).
TET_NODES = np.array(
    [
        [0, 0, 0],
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1],
        [0.5, 0, 0],
        [0.5, 0.5, 0],
        [0, 0.5, 0],
        [0, 0, 0.5],
        [0.5, 0, 0.5],
        [0, 0.5, 0.5],
    ],
    dtype=np.float64,
)


def tiny_job() -> SolveJob:
    mesh = Mesh(
        nodes=TET_NODES,
        elements=np.arange(10, dtype=np.int64).reshape(1, 10),
        element_sets={"part": np.array([0], dtype=np.int64)},
        node_sets={"fixed": np.array([0, 2, 3, 6, 7, 9], dtype=np.int64)},
    )
    return SolveJob(
        mesh=mesh,
        materials=(IsotropicMaterial("resin", 2500.0, 0.35),),
        sections=(Section("part", "resin"),),
        supports=(FixedSupport("fixed", (1, 3)),),
        loads=(
            NodalLoads(np.array([1, 1, 4]), np.array([[1.5, 0, 0], [0.5, 0, 0], [0, 0, -2.0]])),
            NodalLoads(np.array([4]), np.array([[0.0, 0.0, 0.0]])),
        ),
    )


def test_adapter_satisfies_protocol() -> None:
    adapter: SolverAdapter = CalculixAdapter()
    assert adapter.name == "calculix"
    assert adapter.capabilities() == frozenset({"static"})


def test_deck_contents() -> None:
    deck = write_deck(tiny_job())
    lines = deck.splitlines()
    for kw in ("*NODE", "*ELEMENT, TYPE=C3D10", "*STEP", "*STATIC", "*END STEP"):
        assert kw in lines
    n = calculix._num
    assert f"1, {n(0)}, {n(0)}, {n(0)}" in lines
    assert f"10, {n(0)}, {n(0.5)}, {n(0.5)}" in lines
    assert "1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10" in lines  # 0-based -> 1-based, order kept
    assert lines[lines.index("*ELSET, ELSET=E_part") + 1] == "1"
    assert lines[lines.index("*NSET, NSET=N_fixed") + 1] == "1, 3, 4, 7, 8, 10"
    assert lines[lines.index("*MATERIAL, NAME=M_resin") + 1 :][:2] == [
        "*ELASTIC",
        f"{n(2500)}, {n(0.35)}",
    ]
    assert "*SOLID SECTION, ELSET=E_part, MATERIAL=M_resin" in lines
    b = lines.index("*BOUNDARY")
    assert lines[b + 1 : b + 3] == ["N_fixed, 1, 1", "N_fixed, 3, 3"]
    c = lines.index("*CLOAD")
    # duplicates summed (node 2: 1.5 + 0.5), zero components skipped
    assert lines[c + 1 : c + 3] == [f"2, 1, {n(2)}", f"5, 3, {n(-2)}"]
    assert lines[c + 3] == "*NODE FILE"
    assert lines[lines.index("*NODE FILE") + 1] == "U"
    assert lines[lines.index("*EL FILE") + 1] == "S"
    # Point stresses are opt-in (ADR-008): by default nothing goes to job.dat.
    assert lines[lines.index("*EL FILE") :] == ["*EL FILE", "S", "*END STEP"]
    assert "*EL PRINT" not in deck
    assert "EALL" not in deck


def test_point_stresses_add_only_the_print_request_to_the_deck() -> None:
    # With point_stresses the deck defines EALL on *ELEMENT and prints S for it; nothing else
    # changes, and COORD is not asked for (the points are computed from the mesh).
    job = tiny_job()
    default = write_deck(job)
    deck = write_deck(dataclasses.replace(job, point_stresses=True))
    lines = deck.splitlines()
    assert "*ELEMENT, TYPE=C3D10, ELSET=EALL" in lines
    assert lines[lines.index("*EL FILE") :] == [
        "*EL FILE",
        "S",
        "*EL PRINT, ELSET=EALL",
        "S",
        "*END STEP",
    ]
    assert "COORD" not in deck
    undone = deck.replace(", ELSET=EALL\n", "\n", 1).replace("*EL PRINT, ELSET=EALL\nS\n", "")
    assert undone == default


def test_neutral_names_are_mapped_to_safe_unique_names() -> None:
    job = tiny_job()
    one = np.array([0], dtype=np.int64)
    long_name = "face:" + "x" * 100
    mesh = Mesh(
        job.mesh.nodes,
        job.mesh.elements,
        element_sets={"part:beam": one, "Part:Beam": one, "part-beam": one, "3d": one},
        node_sets={
            "support:root": np.array([0, 2]),
            long_name: np.array([1]),
            "face_" + "x" * 100: one,
        },
    )
    mats = (IsotropicMaterial("PLA generic", 3500.0, 0.36),)
    sections = (Section("part:beam", "PLA generic"), Section("Part:Beam", "PLA generic"))
    job = SolveJob(mesh, mats, sections, (FixedSupport("support:root"),))
    names = deck_names(job)
    assert names.element_sets == {
        "part:beam": "E_part_beam",
        "Part:Beam": "E_Part_Beam_2",
        "part-beam": "E_part_beam_3",
        "3d": "E_3d",
    }
    assert names.node_sets["support:root"] == "N_support_root"
    assert len(names.node_sets[long_name]) == 80
    assert names.node_sets["face_" + "x" * 100] == names.node_sets[long_name][:78] + "_2"
    assert names.materials == {"PLA generic": "M_PLA_generic"}
    for table in (names.element_sets, names.node_sets):
        upper = [v.upper() for v in table.values()]
        assert len(set(upper)) == len(upper)
    assert deck_names(job) == names  # deterministic

    lines = write_deck(job).splitlines()
    assert "** ELSET 'part:beam' -> E_part_beam" in lines
    assert "*SOLID SECTION, ELSET=E_Part_Beam_2, MATERIAL=M_PLA_generic" in lines
    assert "N_support_root, 2, 2" in lines
    assert not any(":" in line for line in lines if not line.startswith("**"))
    assert {"ccx_name.elset.part:beam": "E_part_beam"}.items() <= names.as_metadata().items()


def test_long_sets_are_wrapped() -> None:
    job = tiny_job()
    nodes = np.arange(40, dtype=np.int64) % 10
    mesh = Mesh(job.mesh.nodes, job.mesh.elements, job.mesh.element_sets, {"many": nodes})
    deck = write_deck(SolveJob(mesh, job.materials, job.sections))
    lines = deck.splitlines()
    k = lines.index("*NSET, NSET=N_many")
    assert [len(line.split(",")) for line in lines[k + 1 : k + 4]] == [16, 16, 8]


@pytest.mark.parametrize(
    "change",
    [
        {"sections": (Section("nope", "resin"),)},
        {"sections": (Section("part", "steel"),)},
        {"sections": ()},
        {"supports": (FixedSupport("nope"),)},
        {"supports": (FixedSupport("fixed", (4,)),)},
        {"materials": (IsotropicMaterial("resin", 1.0, 0.3), IsotropicMaterial("resin", 2.0, 0.3))},
    ],
)
def test_deck_rejects_invalid_jobs(change: dict[str, object]) -> None:
    job = tiny_job()
    kwargs: dict[str, object] = {
        "mesh": job.mesh,
        "materials": job.materials,
        "sections": job.sections,
        "supports": job.supports,
        "loads": job.loads,
    }
    kwargs.update(change)
    with pytest.raises(SolverError):
        write_deck(SolveJob(**kwargs))  # type: ignore[arg-type]


# ADR-009 catalog TO constants (made-up test values, not material data).
TO = OrthotropicMaterial("to", 3000.0, 2400.0, 1800.0, 0.35, 0.3, 0.3, 1000.0, 700.0, 650.0)
SIDE = MaterialFrame((0.0, 0.0, -1.0), (0.0, 1.0, 0.0))  # build direction (1, 0, 0)


def orthotropic_job(*sections: Section) -> SolveJob:
    job = tiny_job()
    return SolveJob(job.mesh, (TO,), sections or (Section("part", "to", SIDE),), job.supports)


def test_orthotropic_elastic_block_is_engineering_constants() -> None:
    lines = write_deck(orthotropic_job()).splitlines()
    k = lines.index("*MATERIAL, NAME=M_to")
    n = calculix._num
    assert lines[k + 1 : k + 4] == [
        "*ELASTIC, TYPE=ENGINEERING CONSTANTS",
        ", ".join(n(v) for v in (3000, 2400, 1800, 0.35, 0.3, 0.3, 1000, 700)),
        n(650),
    ]
    assert lines[k + 4].startswith("*")  # no density given, no *DENSITY
    with_density = OrthotropicMaterial("to", 3.0, 2.4, 1.8, 0.35, 0.3, 0.3, 1.0, 0.7, 0.65, 1.2e-9)
    job = orthotropic_job()
    lines = write_deck(SolveJob(job.mesh, (with_density,), job.sections)).splitlines()
    k = lines.index("*DENSITY")
    assert lines[k - 1] == n(0.65)
    assert lines[k + 1] == n(1.2e-9)


def test_orientation_card_and_section_reference() -> None:
    lines = write_deck(orthotropic_job()).splitlines()
    n = calculix._num
    k = lines.index("*ORIENTATION, NAME=O_part, SYSTEM=RECTANGULAR")
    assert lines[k + 1] == ", ".join(n(v) for v in (0, 0, -1, 0, 1, 0))
    assert lines[k + 2] == "*SOLID SECTION, ELSET=E_part, MATERIAL=M_to, ORIENTATION=O_part"
    assert lines.index("*MATERIAL, NAME=M_to") < k  # card order verified on ccx 2.21
    assert "** ORIENTATION 'part' -> O_part" in lines


def test_no_orientation_card_without_orientation() -> None:
    deck = write_deck(orthotropic_job(Section("part", "to")))
    assert "*ORIENTATION" not in deck
    assert "ORIENTATION=" not in deck
    assert "*SOLID SECTION, ELSET=E_part, MATERIAL=M_to\n" in deck
    assert deck_names(orthotropic_job(Section("part", "to"))).orientations == {}


def test_orientation_names_are_safe_and_unique() -> None:
    job = tiny_job()
    one = np.array([0], dtype=np.int64)
    sets = {"part:a": one, "Part:A": one, "plain": one, "no frame": one}
    mesh = Mesh(job.mesh.nodes, job.mesh.elements, element_sets=sets)
    frame = MaterialFrame((1.0, 0.0, 0.0), (0.0, 1.0, 0.0))
    sections = (
        Section("part:a", "to", frame),
        Section("Part:A", "to", SIDE),
        Section("plain", "to", frame),
        Section("no frame", "to"),
    )
    job = SolveJob(mesh, (TO,), sections)
    names = deck_names(job)
    assert names.orientations == {"part:a": "O_part_a", "Part:A": "O_Part_A_2", "plain": "O_plain"}
    lines = write_deck(job).splitlines()
    for elset, ccx in names.orientations.items():
        assert f"*ORIENTATION, NAME={ccx}, SYSTEM=RECTANGULAR" in lines
        section = f"*SOLID SECTION, ELSET={names.element_sets[elset]}, MATERIAL=M_to"
        assert f"{section}, ORIENTATION={ccx}" in lines
    assert "*SOLID SECTION, ELSET=E_no_frame, MATERIAL=M_to" in lines
    assert sum(line.startswith("*ORIENTATION") for line in lines) == 3


def test_orientation_names_in_metadata() -> None:
    meta = deck_names(orthotropic_job()).as_metadata()
    assert meta["ccx_name.orientation.part"] == "O_part"
    assert meta["ccx_name.elset.part"] == "E_part"
    assert meta["ccx_name.material.to"] == "M_to"
    assert not any(
        k.startswith("ccx_name.orientation.") for k in deck_names(tiny_job()).as_metadata()
    )


def test_missing_executable_raises(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    monkeypatch.setenv(CCX_ENV_VAR, str(tmp_path / "no-such-ccx"))
    with pytest.raises(SolverError, match=CCX_ENV_VAR):
        CalculixAdapter().run(tiny_job(), tmp_path)


FRD_SNIPPET = """\
    1C
    1UVERSION           Version 2.21
    2C                             3                                     1
 -1         1 0.00000E+00 0.00000E+00 0.00000E+00
 -3
    1PSTEP                         1           1           1
  100CL  101 1.000000000           3                     0    1           1
 -4  DISP        4    1
 -5  D1          1    2    1    0
 -5  ALL         1    2    0    0    1ALL
 -1         1 9.00000E+00 9.00000E+00 9.00000E+00
 -3
    1PSTEP                         2           1           1
  100CL  101 1.000000000           3                     0    1           1
 -4  DISP        4    1
 -5  D1          1    2    1    0
 -1         1 0.00000E+00 0.00000E+00 0.00000E+00
 -1         3 1.25000E-02-3.00000E-03 4.50000E+01
 -1         2-1.00000E+00-2.00000E+00-3.00000E+00
 -3
  100CL  101 1.000000000           3                     0    1           1
 -4  STRESS      6    1
 -5  SXX         1    4    1    1
 -1         1 1.00000E+00 2.00000E+00 3.00000E+00 4.00000E+00 5.00000E+00 6.00000E+00
 -1         3-1.00000E+00-2.00000E+00-3.00000E+00-4.00000E+00-5.00000E+00-6.00000E+00
 -3
  100CL  101 1.000000000           3                     0    1           1
 -4  ERROR       1    1
 -1         1 1.00000E+00
 -3
 9999
"""


def test_parse_frd_last_step_and_packed_fields() -> None:
    res = parse_frd(FRD_SNIPPET)
    assert res.version == "Version 2.21"
    disp = res.blocks["DISP"]  # second DISP block (last step) wins
    assert disp[1] == [0.0, 0.0, 0.0]
    assert disp[2] == [-1.0, -2.0, -3.0]
    assert disp[3] == [1.25e-2, -3.0e-3, 45.0]
    assert res.blocks["STRESS"][3] == [-1.0, -2.0, -3.0, -4.0, -5.0, -6.0]
    assert 2 not in res.blocks["STRESS"]


@given(st.floats(allow_nan=False, allow_infinity=False, width=64))
def test_deck_numbers_fit_calculix_20_character_fields(x: float) -> None:
    # CalculiX reads numeric fields of at most 20 characters; longer reprs abort the run.
    text = calculix._num(x)
    assert len(text) <= 20
    assert float(text) == pytest.approx(x, rel=1e-12, abs=0.0)


def test_element_nodes_without_results_are_a_solver_error() -> None:
    mesh = tiny_job().mesh
    disp = np.zeros((mesh.n_nodes, 3))
    stress = np.zeros((mesh.n_nodes, 6))
    calculix.check_results_complete(mesh, disp, stress)  # complete: no error
    disp[4] = np.nan
    with pytest.raises(SolverError, match="1 mesh node"):
        calculix.check_results_complete(mesh, disp, stress)


def test_number_format_is_pinned() -> None:
    assert calculix._num(0.5) == "5.000000000000e-01"
    assert calculix._num(-8.526512829121204e-17) == "-8.526512829121e-17"


def test_user_set_named_eall_does_not_extend_the_print_set() -> None:
    # ccx appends to a set defined twice, which would print every element of it twice (ADR-008
    # finding 8). User sets are always prefixed, so EALL is defined once, on *ELEMENT.
    job = tiny_job()
    one = np.array([0], dtype=np.int64)
    mesh = Mesh(job.mesh.nodes, job.mesh.elements, {"part": one, "EALL": one, "eall": one})
    deck = write_deck(SolveJob(mesh, job.materials, job.sections, point_stresses=True))
    lines = deck.splitlines()
    defines = [line for line in lines if not line.startswith("**") and "ELSET=" in line.upper()]
    defined = [line.upper().split("ELSET=")[1].split(",")[0] for line in defines]
    assert defined.count("EALL") == 2  # *ELEMENT and the *EL PRINT card
    assert "*ELEMENT, TYPE=C3D10, ELSET=EALL" in lines
    assert not any(line.upper().startswith("*ELSET, ELSET=EALL") for line in lines)
    assert "*ELSET, ELSET=E_EALL" in lines


# A ccx 2.21 job.dat in its exact layout (element I10, 1X, point I3, values 1X,E13.6): one step,
# two increments; the second holds distinct values in every column. Elements 1 and 3 only.
DAT_SNIPPET: str = """\

                        S T E P       1


                                INCREMENT     1


 stresses (elem, integ.pnt.,sxx,syy,szz,sxy,sxz,syz) for set EALL and time  0.5000000E+00

         1   1  9.000000E+00  9.000000E+00  9.000000E+00  9.000000E+00  9.000000E+00  9.000000E+00
         1   2  9.000000E+00  9.000000E+00  9.000000E+00  9.000000E+00  9.000000E+00  9.000000E+00
         1   3  9.000000E+00  9.000000E+00  9.000000E+00  9.000000E+00  9.000000E+00  9.000000E+00
         1   4  9.000000E+00  9.000000E+00  9.000000E+00  9.000000E+00  9.000000E+00  9.000000E+00
         3   1  9.000000E+00  9.000000E+00  9.000000E+00  9.000000E+00  9.000000E+00  9.000000E+00
         3   2  9.000000E+00  9.000000E+00  9.000000E+00  9.000000E+00  9.000000E+00  9.000000E+00
         3   3  9.000000E+00  9.000000E+00  9.000000E+00  9.000000E+00  9.000000E+00  9.000000E+00
         3   4  9.000000E+00  9.000000E+00  9.000000E+00  9.000000E+00  9.000000E+00  9.000000E+00

 global coordinates (elem, integ.pnt.,x,y,z) for set EALL and time  0.5000000E+00

         1   1  9.000000E+00  9.000000E+00  9.000000E+00
         1   2  9.000000E+00  9.000000E+00  9.000000E+00
         1   3  9.000000E+00  9.000000E+00  9.000000E+00
         1   4  9.000000E+00  9.000000E+00  9.000000E+00
         3   1  9.000000E+00  9.000000E+00  9.000000E+00
         3   2  9.000000E+00  9.000000E+00  9.000000E+00
         3   3  9.000000E+00  9.000000E+00  9.000000E+00
         3   4  9.000000E+00  9.000000E+00  9.000000E+00


                                INCREMENT     2


 stresses (elem, integ.pnt.,sxx,syy,szz,sxy,sxz,syz) for set EALL and time  0.1000000E+01

         1   1  1.110000E+00  1.120000E+00  1.130000E+00  1.140000E+00  1.150000E+00  1.160000E+00
         1   2  1.210000E+00  1.220000E+00  1.230000E+00  1.240000E+00  1.250000E+00  1.260000E+00
         1   3  1.310000E+00  1.320000E+00  1.330000E+00  1.340000E+00  1.350000E+00  1.360000E+00
         1   4  1.410000E+00  1.420000E+00  1.430000E+00  1.440000E+00  1.450000E+00  1.460000E+00
         3   1 -2.110000E+01 -2.120000E+01 -2.130000E+01 -2.140000E+01 -2.150000E+01 -2.160000E+01
         3   2 -2.210000E+01 -2.220000E+01 -2.230000E+01 -2.240000E+01 -2.250000E+01 -2.260000E+01
         3   3 -2.310000E+01 -2.320000E+01 -2.330000E+01 -2.340000E+01 -2.350000E+01 -2.360000E+01
         3   4 -2.410000E+01 -2.420000E+01 -2.430000E+01 -2.440000E+01 -2.450000E+01 -2.460000E+01

 global coordinates (elem, integ.pnt.,x,y,z) for set EALL and time  0.1000000E+01

         1   1  1.110000E+02  1.120000E+02  1.130000E+02
         1   2  1.210000E+02  1.220000E+02  1.230000E+02
         1   3  1.310000E+02  1.320000E+02  1.330000E+02
         1   4  1.410000E+02  1.420000E+02  1.430000E+02
         3   1  3.110000E+02  3.120000E+02  3.130000E+02
         3   2  3.210000E+02  3.220000E+02  3.230000E+02
         3   3  3.310000E+02  3.320000E+02  3.330000E+02
         3   4  3.410000E+02  3.420000E+02  3.430000E+02
"""

# Element 1, point 2 of the last increment, in .dat order sxx syy szz sxy sxz syz.
DAT_E1_P2 = [1.21, 1.22, 1.23, 1.24, 1.25, 1.26]


def test_parse_dat_keeps_the_last_increment_in_neutral_order() -> None:
    dat = parse_dat(DAT_SNIPPET)
    assert set(dat.blocks) == {"stresses", "global coordinates"}
    stress = dat.blocks["stresses"]
    np.testing.assert_array_equal(stress.elements, [1, 3])
    assert stress.values.shape == (2, 4, 6)
    # .dat prints sxz before syz; the neutral order is xx yy zz xy yz zx (columns 0 1 2 3 5 4).
    s = DAT_E1_P2
    np.testing.assert_array_equal(stress.values[0, 1], [s[0], s[1], s[2], s[3], s[5], s[4]])
    np.testing.assert_array_equal(stress.values[1, 3], [-24.1, -24.2, -24.3, -24.4, -24.6, -24.5])
    coord = dat.blocks["global coordinates"]
    np.testing.assert_array_equal(coord.elements, [1, 3])
    assert coord.values.shape == (2, 4, 3)
    np.testing.assert_array_equal(coord.values[1, 2], [331.0, 332.0, 333.0])
    assert not (stress.values == 9.0).any()  # the first increment is gone


def test_parse_dat_of_an_empty_file_has_no_blocks() -> None:
    assert parse_dat("").blocks == {}
    assert parse_dat("\n                        S T E P       1\n").blocks == {}


def test_parse_dat_rows_may_come_in_any_order() -> None:
    head, rows = DAT_SNIPPET.rsplit("\n\n", 1)
    shuffled = "\n".join(reversed(rows.splitlines()))
    coord = parse_dat(f"{head}\n\n{shuffled}\n").blocks["global coordinates"]
    np.testing.assert_array_equal(coord.elements, [1, 3])
    np.testing.assert_array_equal(coord.values[0, 0], [111.0, 112.0, 113.0])


def _last_coord_rows() -> tuple[str, list[str]]:
    head, rows = DAT_SNIPPET.rsplit("\n\n", 1)
    return head, rows.splitlines()


def test_parse_dat_rejects_duplicate_rows() -> None:
    # A set defined twice makes ccx print its elements twice, without a warning (finding 8).
    head, rows = _last_coord_rows()
    text = "\n".join([head, "", *rows, *rows[:4]])
    with pytest.raises(SolverError, match="duplicate"):
        parse_dat(text)


def test_parse_dat_rejects_a_missing_point() -> None:
    head, rows = _last_coord_rows()
    text = "\n".join([head, "", *rows[:-1]])  # element 3 lacks point 4
    with pytest.raises(SolverError, match="points"):
        parse_dat(text)


def test_parse_dat_rejects_a_point_number_out_of_order() -> None:
    head, rows = _last_coord_rows()
    rows[-1] = rows[-1].replace("         3   4", "         3   5")  # points 1, 2, 3, 5
    with pytest.raises(SolverError, match="points"):
        parse_dat("\n".join([head, "", *rows]))


def test_parse_dat_pins_the_stress_column_order() -> None:
    # The xz/yz reorder is only right for the order ccx 2.21 prints; another order is refused.
    text = DAT_SNIPPET.replace("sxy,sxz,syz", "sxy,syz,sxz")
    with pytest.raises(SolverError, match="columns"):
        parse_dat(text)


def test_parse_dat_rejects_a_short_row() -> None:
    head, rows = _last_coord_rows()
    rows[0] = rows[0][:-14]  # one value missing
    with pytest.raises(SolverError, match="malformed"):
        parse_dat("\n".join([head, "", *rows]))


def test_parse_dat_rejects_point_numbers_below_one() -> None:
    head, rows = _last_coord_rows()
    for bad in ("   0", "  -1"):
        changed = [*rows]
        changed[0] = changed[0][:10] + bad + changed[0][14:]  # element 1: points 0, 2, 3, 4
        with pytest.raises(SolverError, match=r"point number -?[01]; point numbers start at 1"):
            parse_dat("\n".join([head, "", *changed]))


def test_parse_dat_reads_bytes_and_crlf_lines() -> None:
    windows = DAT_SNIPPET.replace("\n", "\r\n").encode()
    stress = parse_dat(windows).blocks["stresses"]
    s = DAT_E1_P2
    np.testing.assert_array_equal(stress.values[0, 1], [s[0], s[1], s[2], s[3], s[5], s[4]])


@pytest.mark.parametrize("word", ["inf", "Infinity", "-Infinity", "NaN"])
def test_parse_dat_keeps_non_finite_values_for_the_completeness_check(word: str) -> None:
    text = DAT_SNIPPET.replace("-2.410000E+01", word, 1)
    values = parse_dat(text).blocks["stresses"].values
    assert not np.isfinite(values[1, 3, 0])


def three_tets(offset: float = 0.0) -> Mesh:
    """Three copies of the unit tet10 (enough elements for ccx ids 1 to 3)."""
    return Mesh(
        nodes=TET_NODES + offset,
        elements=np.tile(np.arange(10, dtype=np.int64), (3, 1)),
    )


# C3D10 4-point rule: point k has barycentric coordinate B at corner k and A at the other three.
RULE_A = (5.0 - np.sqrt(5.0)) / 20.0
RULE_B = (5.0 + 3.0 * np.sqrt(5.0)) / 20.0


def test_rule_points_of_a_straight_tet() -> None:
    a, b = RULE_A, RULE_B
    expected = [[a, a, a], [b, a, a], [a, b, a], [a, a, b]]  # ccx's COORD for this tet
    points = rule_points(three_tets())
    assert points.shape == (3, 4, 3)
    np.testing.assert_allclose(points[0], expected, rtol=0, atol=1e-15)
    np.testing.assert_allclose(
        points[0],
        [
            [0.1381966] * 3,
            [0.5854102, 0.1381966, 0.1381966],
            [0.1381966, 0.5854102, 0.1381966],
            [0.1381966, 0.1381966, 0.5854102],
        ],
        atol=5e-8,
    )


def test_rule_points_keep_full_precision_far_from_the_origin() -> None:
    # COORD prints 7 significant digits of the absolute position (0.5 um at 1000 mm); the
    # points computed from the mesh do not lose that.
    near, far = rule_points(three_tets()), rule_points(three_tets(offset=1000.0))
    np.testing.assert_allclose(far - 1000.0, near, rtol=0, atol=1e-12)


def test_point_stresses_are_in_mesh_order_with_nan_for_missing_elements() -> None:
    mesh = three_tets()
    ps = point_stresses(parse_dat(DAT_SNIPPET), mesh)
    assert isinstance(ps, PointStresses)
    assert ps.stresses.shape == (3, 4, 6)
    s = DAT_E1_P2
    np.testing.assert_array_equal(ps.stresses[0, 1], [s[0], s[1], s[2], s[3], s[5], s[4]])
    np.testing.assert_array_equal(ps.stresses[2, 0], [-21.1, -21.2, -21.3, -21.4, -21.6, -21.5])
    assert np.isnan(ps.stresses[1]).all()  # ccx element 2 was not printed
    # The coordinates come from the mesh, not from a COORD block.
    np.testing.assert_array_equal(ps.points, rule_points(mesh))


def test_point_stresses_need_only_the_stress_block() -> None:
    dat = parse_dat(DAT_SNIPPET)
    without_coord = DatResults({"stresses": dat.blocks["stresses"]})
    assert point_stresses(without_coord, three_tets()).n_points == 4
    without_stress = DatResults({"global coordinates": dat.blocks["global coordinates"]})
    with pytest.raises(SolverError, match="'stresses' block"):
        point_stresses(without_stress, three_tets())


def test_point_stresses_check_the_point_rule() -> None:
    block = DatBlock(elements=np.array([1]), values=np.zeros((1, 10, 6)))
    with pytest.raises(
        SolverError, match=r"10 integration points per element; C3D10 \(tet4\) has 4"
    ):
        point_stresses(DatResults({"stresses": block}), three_tets())


def test_point_stresses_reject_unknown_elements() -> None:
    two = Mesh(TET_NODES, np.tile(np.arange(10, dtype=np.int64), (2, 1)))
    with pytest.raises(SolverError, match=r"element ids outside 1\.\.2"):
        point_stresses(parse_dat(DAT_SNIPPET), two)  # ccx element 3 does not exist


@pytest.mark.parametrize("content", [None, ""])
def test_missing_or_empty_dat_is_a_solver_error(tmp_path: Path, content: str | None) -> None:
    # ccx writes an empty job.dat and exits 0 when *EL PRINT is missing or names no set.
    path = tmp_path / "job.dat"
    if content is not None:
        path.write_text(content)
    with pytest.raises(SolverError, match=r"job\.dat is missing or empty") as info:
        read_point_stresses(path, three_tets(), log="ccx log")
    assert info.value.log == "ccx log"


def test_read_point_stresses_reads_the_file(tmp_path: Path) -> None:
    path = tmp_path / "job.dat"
    path.write_text(DAT_SNIPPET)
    ps = read_point_stresses(path, three_tets())
    np.testing.assert_array_equal(ps.stresses[2, 3], [-24.1, -24.2, -24.3, -24.4, -24.6, -24.5])


def test_dat_parse_errors_carry_the_log(tmp_path: Path) -> None:
    path = tmp_path / "job.dat"
    head, rows = _last_coord_rows()
    path.write_text("\n".join([head, "", *rows, *rows[:4]]))
    with pytest.raises(SolverError, match="duplicate") as info:
        read_point_stresses(path, three_tets(), log="ccx log")
    assert info.value.log == "ccx log"


def test_elements_without_point_results_are_a_solver_error() -> None:
    mesh = tiny_job().mesh
    disp = np.zeros((mesh.n_nodes, 3))
    stress = np.zeros((mesh.n_nodes, 6))
    ps = PointStresses(np.zeros((1, 4, 3)), np.zeros((1, 4, 6)))
    calculix.check_results_complete(mesh, disp, stress, ps)  # complete: no error
    for bad in (np.nan, np.inf, -np.inf):
        broken = PointStresses(ps.points, ps.stresses.copy())
        broken.stresses[0, 2, 5] = bad
        with pytest.raises(SolverError, match=r"no integration-point result for 1 element\(s\)"):
            calculix.check_results_complete(mesh, disp, stress, broken, "ccx log")
    with pytest.raises(SolverError, match=r"for 0 element\(s\); the mesh has 1 element\(s\)"):
        calculix.check_results_complete(
            mesh, disp, stress, PointStresses(np.zeros((0, 4, 3)), np.zeros((0, 4, 6)))
        )


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf])
def test_non_finite_nodal_results_are_a_solver_error(bad: float) -> None:
    mesh = tiny_job().mesh
    disp = np.zeros((mesh.n_nodes, 3))
    stress = np.zeros((mesh.n_nodes, 6))
    stress[4, 2] = bad
    with pytest.raises(SolverError, match=r"no result for 1 mesh node\(s\)") as info:
        calculix.check_results_complete(mesh, disp, stress, log="ccx log")
    assert info.value.log == "ccx log"


def test_infinity_in_job_dat_is_caught_by_the_completeness_check(tmp_path: Path) -> None:
    path = tmp_path / "job.dat"
    path.write_text(DAT_SNIPPET.replace("-2.410000E+01", "Infinity", 1))
    mesh = three_tets()
    ps = read_point_stresses(path, mesh)
    disp, stress = np.zeros((mesh.n_nodes, 3)), np.zeros((mesh.n_nodes, 6))
    ps.stresses[1] = 0.0  # fill the element ccx did not print, so only the infinity is left
    with pytest.raises(
        SolverError, match=r"no integration-point result for 1 element\(s\), e\.g\. element 2"
    ):
        calculix.check_results_complete(mesh, disp, stress, ps)
