"""CalculiX adapter unit tests that do not need the ccx executable."""

from pathlib import Path

import numpy as np
import pytest
from hypothesis import given
from hypothesis import strategies as st

from tessera.model.mesh import Mesh
from tessera.model.solver import (
    FixedSupport,
    IsotropicMaterial,
    NodalLoads,
    PointStresses,
    Section,
    SolveJob,
    SolverAdapter,
    SolverError,
)
from tessera.solvers import CalculixAdapter, calculix
from tessera.solvers.calculix import (
    CCX_ENV_VAR,
    DatResults,
    deck_names,
    parse_dat,
    parse_frd,
    point_stresses,
    read_point_stresses,
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
    for kw in ("*NODE", "*ELEMENT, TYPE=C3D10, ELSET=EALL", "*STEP", "*STATIC", "*END STEP"):
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
    # ADR-008: integration-point stresses and their coordinates go to job.dat, for every element.
    k = lines.index("*EL FILE") + 2
    print_cards = ["*EL PRINT, ELSET=EALL", "S", "*EL PRINT, ELSET=EALL", "COORD", "*END STEP"]
    assert lines[k : k + 5] == print_cards


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
    lines = write_deck(SolveJob(mesh, job.materials, job.sections)).splitlines()
    defines = [line for line in lines if not line.startswith("**") and "ELSET=" in line.upper()]
    defined = [line.upper().split("ELSET=")[1].split(",")[0] for line in defines]
    assert defined.count("EALL") == 3  # *ELEMENT and the two *EL PRINT cards
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


def test_point_stresses_are_in_mesh_order_with_nan_for_missing_elements() -> None:
    ps = point_stresses(parse_dat(DAT_SNIPPET), n_elements=3)
    assert isinstance(ps, PointStresses)
    assert ps.stresses.shape == (3, 4, 6)
    assert ps.points.shape == (3, 4, 3)
    s = DAT_E1_P2
    np.testing.assert_array_equal(ps.stresses[0, 1], [s[0], s[1], s[2], s[3], s[5], s[4]])
    np.testing.assert_array_equal(ps.points[2, 0], [311.0, 312.0, 313.0])  # ccx element 3
    assert np.isnan(ps.stresses[1]).all()  # ccx element 2 was not printed
    assert np.isnan(ps.points[1]).all()


def test_point_stresses_need_both_blocks() -> None:
    dat = parse_dat(DAT_SNIPPET)
    for name in ("stresses", "global coordinates"):
        partial = DatResults({k: v for k, v in dat.blocks.items() if k != name})
        with pytest.raises(SolverError, match=name):
            point_stresses(partial, n_elements=3)


def test_point_stresses_check_the_point_rule() -> None:
    with pytest.raises(SolverError, match="4 integration points"):
        point_stresses(parse_dat(DAT_SNIPPET), n_elements=3, n_points=10)


def test_point_stresses_reject_unknown_elements() -> None:
    with pytest.raises(SolverError, match="element"):
        point_stresses(parse_dat(DAT_SNIPPET), n_elements=2)  # ccx element 3 does not exist


@pytest.mark.parametrize("content", [None, ""])
def test_missing_or_empty_dat_is_a_solver_error(tmp_path: Path, content: str | None) -> None:
    # ccx writes an empty job.dat and exits 0 when *EL PRINT is missing or names no set.
    path = tmp_path / "job.dat"
    if content is not None:
        path.write_text(content)
    with pytest.raises(SolverError, match=r"job\.dat") as info:
        read_point_stresses(path, n_elements=3, log="ccx log")
    assert info.value.log == "ccx log"


def test_read_point_stresses_reads_the_file(tmp_path: Path) -> None:
    path = tmp_path / "job.dat"
    path.write_text(DAT_SNIPPET)
    ps = read_point_stresses(path, n_elements=3)
    np.testing.assert_array_equal(ps.points[0, 3], [141.0, 142.0, 143.0])


def test_dat_parse_errors_carry_the_log(tmp_path: Path) -> None:
    path = tmp_path / "job.dat"
    head, rows = _last_coord_rows()
    path.write_text("\n".join([head, "", *rows, *rows[:4]]))
    with pytest.raises(SolverError, match="duplicate") as info:
        read_point_stresses(path, n_elements=3, log="ccx log")
    assert info.value.log == "ccx log"


def test_elements_without_point_results_are_a_solver_error() -> None:
    mesh = tiny_job().mesh
    disp = np.zeros((mesh.n_nodes, 3))
    stress = np.zeros((mesh.n_nodes, 6))
    ps = PointStresses(np.zeros((1, 4, 3)), np.zeros((1, 4, 6)))
    calculix.check_results_complete(mesh, disp, stress, ps)  # complete: no error
    ps.stresses[0, 2, 5] = np.nan
    with pytest.raises(SolverError, match="1 element"):
        calculix.check_results_complete(mesh, disp, stress, ps)
    with pytest.raises(SolverError, match="1 element"):
        calculix.check_results_complete(
            mesh, disp, stress, PointStresses(np.zeros((0, 4, 3)), np.zeros((0, 4, 6)))
        )
