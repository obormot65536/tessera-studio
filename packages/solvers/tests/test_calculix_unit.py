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
    MaterialFrame,
    NodalLoads,
    OrthotropicMaterial,
    Section,
    SolveJob,
    SolverAdapter,
    SolverError,
)
from tessera.solvers import CalculixAdapter, calculix
from tessera.solvers.calculix import CCX_ENV_VAR, deck_names, parse_frd, write_deck

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
