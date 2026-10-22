# ADR-009: Orthotropic materials and print orientation (v0.5, slice 1)

**Status:** Accepted (maintainer, 2026-10-22)
**Deciders:** Maintainer

## Context

The v0.5 gate (`docs/development.md`) starts with "orthotropic patch tests exact" (TO-01), and
ADR-004 item 3 asks to extend the solve job additively with orientations and orthotropic
materials. This slice does only that: orthotropic linear elasticity with a frame taken from the
print setup, enough to pass TO-01. Off-axis checks (TO-02) are the next slice; failure criteria,
slicer import, infill homogenization and real filament data are out of scope.

What exists today:

- `tessera.model.solver` has `IsotropicMaterial` only; `Section` is `(element_set, material)`;
  `SolveJob.materials: tuple[IsotropicMaterial, ...]`; `FixedSupport` holds displacements at zero
  only (no prescribed non-zero values).
- `tessera.materials.library.get_isotropic(material_id)` holds two isotropic metals;
  `PartInstance.material_id` can name only library materials.
- `tessera.printsetup` is an empty package. `PrintSetup(process, slicer_profile,
  build_direction=(0, 0, 1))` exists in `tessera.model.project`, but nothing reads it.
- `physics.build_static_job(study, part, mesh, material: IsotropicMaterial)`; `write_deck` writes
  `*ELASTIC` (isotropic) and `*SOLID SECTION` without orientation.
- `SCHEMA_VERSION = 5`. Every additive field so far has bumped it; `extra="forbid"` means an
  older Tessera refuses newer fields.
- The catalog's axis convention (TO section): 1 = raster/extrusion, 2 = across the road in the
  layer, 3 = build (Z).

## Decision

### D1. Store engineering constants

`OrthotropicMaterial(name, e1, e2, e3, nu12, nu13, nu23, g12, g13, g23, density=None)` in
`tessera.model.solver`, with `Material = IsotropicMaterial | OrthotropicMaterial`.

- **Why not a stiffness matrix:** datasheets, FDM papers and coupon tests report engineering
  constants; CalculiX takes them directly (`*ELASTIC, TYPE=ENGINEERING CONSTANTS`), so the deck
  carries the user's numbers without a conversion; Voigt order and the shear factor of 2 are easy
  to get wrong in a reviewed matrix. A general anisotropic type (21 constants, from
  homogenization) can come later as a separate type.
- **Poisson convention:** nu_ij = -eps_j / eps_i under uniaxial sigma_i (the major ratios, as in
  Jones, Abaqus and CalculiX). Minor ratios are derived (nu_ji = nu_ij E_j / E_i), never stored, so
  the data cannot contradict itself. FDM papers mix nu_xz and nu_zx, so the curator records which
  one a paper measured.
- **Validation** (exact math, no threshold to approve): all E and G > 0, finite; the 6 × 6
  compliance S is positive definite (`eigvalsh(S).min() > 0`). The error names the usual bounds
  |nu_ij| < sqrt(E_i/E_j) and Delta = 1 - nu12 nu21 - nu23 nu32 - nu13 nu31 - 2 nu21 nu32 nu13 > 0.
  The catalog's TO set gives Delta = 0.7427.
- **Pure functions** in a new `tessera/model/elasticity.py`: `orthotropic_compliance(...)` (Voigt
  xx yy zz xy yz zx, the `SolveResult.stresses` order, engineering shear strains),
  `check_orthotropic(...)`, and `rotate_stress(stress, frame)` for the later failure criteria.
- **Lookup** in `tessera.materials`: `get_material(material_id, defined=()) -> Material`; project
  definitions (D3) win over the library. `get_isotropic` stays. No printed-material library
  entries in this slice.
- **XY/Z data:** datasets with only XY and Z values are the transversely isotropic case E1 = E2,
  nu13 = nu23, G13 = G23, G12 = E1 / (2 (1 + nu12)). A helper can come later; not needed here.

### D2. Orientation: one frame per section

- `MaterialFrame(axis_1, axis_2)` in `tessera.model.solver`: unit, orthogonal within 1e-9;
  `axis_3 = axis_1 × axis_2`; `rotation` (rows = material axes in part coordinates);
  `from_axes(a, b)` does Gram-Schmidt.
- `Section.orientation: MaterialFrame | None = None` (additive; `None` = global axes). CalculiX
  orientations are per element set anyway, so region-wise frames (TO-05) or curved rasters become
  more element sets, not a new table. No CalculiX keyword enters `model` (ADR-004).
- `printsetup.material_frame(build_direction, raster_angle=0.0) -> MaterialFrame`:
  1. e3 = normalized build direction (the part-coordinate vector pointing up on the printer).
  2. R_min = the smallest rotation taking printer +Z onto e3 (Rodrigues); e3 = -Z uses 180° about
     printer X. A zero vector raises `ValueError`.
  3. e1 = R_min · (1, 0, 0), then rotated about e3 by `raster_angle` (degrees, counter-clockwise
     looking down e3); e2 = e3 × e1.

  Any printer placement splits into a minimal tilt plus a yaw about the build axis, so build
  direction plus raster angle loses nothing. Example: build (1, 0, 0) gives e1 = (0, 0, -1),
  e2 = (0, 1, 0), e3 = (1, 0, 0).
- `core.study.solve`: an orthotropic material requires `part.print_setup`, else
  `InputError("part 'x' uses orthotropic material 'm' and needs a print setup (build
  direction)")`; silently assuming global axes is a wrong-orientation trap. Isotropic parts get
  `orientation=None` even with a print setup.
- Deck (`write_deck`):

  ```
  *MATERIAL, NAME=M_...
  *ELASTIC, TYPE=ENGINEERING CONSTANTS
  E1, E2, E3, nu12, nu13, nu23, G12, G13
  G23
  *ORIENTATION, NAME=O_..., SYSTEM=RECTANGULAR
  a1, a2, a3, b1, b2, b3
  *SOLID SECTION, ELSET=E_..., MATERIAL=M_..., ORIENTATION=O_...
  ```

  a = axis_1, b = axis_2. Names come from `_safe_names(..., "O_")`, recorded in `DeckNames` and
  the metadata (`ccx_name.orientation.<elset>`). Engineering constants rather than `ORTHO`
  (D1111…), which would put an untested compliance inversion inside the adapter. Decks without
  orientation stay byte-identical. `calculix-specialist` confirms that the frd stresses come back
  in global axes with an `*ORIENTATION` active, and the card order; the tilted TO-01 cases catch
  either way.

### D3. Schema 6: materials defined in the project

TO-01 must run through the public API like TS-03, and a project has no way to name orthotropic
constants today.

- `Project.materials: list[MaterialDefinition] = []`, a union discriminated on `kind`:
  `IsotropicDefinition(kind="isotropic", id, youngs_modulus, poisson_ratio, density=None)` and
  `OrthotropicDefinition(kind="orthotropic", id, e1, …, g23, density=None)` (validated with
  `check_orthotropic`). Ids are unique; a project id may not reuse a built-in id
  (`edit.add_material` refuses it), and at load time a project definition wins, so a file keeps its
  meaning if a later library adds the same id.
- `PrintSetup.raster_angle: float = 0.0` (degrees).
- `SCHEMA_VERSION = 6`, both fields in one PR so the schema bumps once. Additive with defaults:
  the migration is the identity, pinned by `test_schema_5_project_loads_unchanged` (ADR-003).
- `edit`: `add_material`, `set_print_setup`; `_check_material`, `set_material` and
  `_check_gravity_densities` resolve project materials too.
- Hard rule 6 covers numbers Tessera ships (library, criteria, tolerances), not values a user
  types in. This is where "X/Y/Z coupons tested" enters: users type their own coupon data.

### D4. Results: von Mises stays, labelled

Von Mises is not a strength check for a printed part (interlayer strength is often under half the
in-plane strength), but it still shows where stress concentrates, and ADR-007's flags and two-mesh
check work on any scalar.

- Peaks keep their meaning. `StudyResult.orthotropic: bool = field(default=False, kw_only=True)`,
  set by `post.summarize` when any job material is orthotropic.
- `tessera solve` prints, before the stress lines (wording to approve with this ADR):
  "note: this part uses an orthotropic (layer-wise) material; von Mises stress shows where the
  stress is, but it is not a strength check: printed parts are weaker across layers. Layer-wise
  failure checks are planned for v0.5."
- ADR-007's "clear/settled = compare with a strength" does not apply to orthotropic parts; the
  MCP tools must expose `orthotropic` and say so.
- The failure-criteria slice depends on ADR-008 (integration-point stresses, proposed: averaging
  across material or orientation boundaries is wrong), on `rotate_stress`, and on per-direction
  strength data from the curator.

## Test-first plan

Unit tests (fail before the code):

- `model/tests/test_elasticity.py`: compliance equals the hand-built matrix and is symmetric
  (reciprocity nu12/E1 = nu21/E2); the isotropic limit inverts to Lamé λ, μ; rejections with
  message tests (nu12 = 1.2 with E1 = E2; nu12 = nu13 = nu23 = 0.6 with equal moduli; G ≤ 0;
  E ≤ 0; NaN); `MaterialFrame` validation, `from_axes`, right-handedness; `rotate_stress`
  round-trips and leaves an isotropic state unchanged.
- `model/tests/test_project.py`: orthotropic definition round-trip; invalid definition and
  duplicate ids rejected; `raster_angle` round-trip; `test_schema_5_project_loads_unchanged`.
- `printsetup/tests/test_orientation.py`: Z gives identity; (1, 0, 0) gives the example above;
  -Z gives 180° about X; Z with raster 90 gives e1 = y; an arbitrary tilt is orthonormal,
  right-handed, e3 = b; zero vector raises.
- `solvers/tests/test_calculix_unit.py`: exact `*ELASTIC, TYPE=ENGINEERING CONSTANTS` text (8 + 1
  values); `*ORIENTATION` card and `ORIENTATION=` on its section; none when `orientation is None`
  (`test_deck_contents` unchanged); safe unique names; metadata.
- `solvers/tests/test_calculix_orthotropic.py` (ccx): a direct `SolveJob` on a small cube, load
  along material 1, identity frame and the (1, 0, 0) frame. The early answer on the stress output
  frame, before any core wiring.
- `materials`, `physics`, `core`, `post`, `cli`: project definitions resolve before the library;
  `Section(orientation=frame)` and orthotropic density in gravity; orthotropic part without a
  print setup raises; `add_material` refuses built-in ids; the CLI note appears only for
  orthotropic parts.

TO-01 benchmark (`benchmarks/cases/test_orthotropic_patch_test.py`, by `fea-verifier`; reference
`orthotropic_patch_test.toml` proposed in the PR, created by the maintainer). Cube L = 10 mm from
`box_step`, unstructured mesh (`_assert_unstructured`), gate 3.0 mm, sweep 3.0 / 1.3 / 0.9 mm as
TS-03. σ = τ = 10 MPa as a `ForceLoad` of σ L². Catalog TO constants: E1 = 3000, E2 = 2400,
E3 = 1800, G12 = 1000, G13 = 700, G23 = 650 MPa, nu12 = 0.35, nu13 = nu23 = 0.30 (so nu21 = 0.28,
nu31 = 0.18, nu32 = 0.225).

(a) Axial, rollers as TS-03 (u_x = 0 on x = 0, u_y = 0 on y = 0, u_z = 0 on z = 0), load on the
opposite face; exact field u = eps · x at every node, σ_load = 10 and the other stresses 0.

| Case | Build | Load | eps (x, y, z) × 1e-3 | Far-face u (mm) |
|---|---|---|---|---|
| A1 | Z, raster 0 | +x, material 1 | 3.333333, -1.166667, -1.000000 | 0.0333333333 |
| A2 | Z | +y, material 2 | -1.166667, 4.166667, -1.250000 | 0.0416666667 |
| A3 | Z | +z, material 3 | -1.000000, -1.250000, 5.555556 | 0.0555555556 |
| A1-side | (1, 0, 0) | +x, material 3 | 5.555556, -1.250000, -1.000000 | 0.0555555556 |
| A2-side | (1, 0, 0) | +y, material 2 | -1.250000, 4.166667, -1.166667 | 0.0416666667 |
| A3-side | (1, 0, 0) | +z, material 1 | -1.000000, -1.166667, 3.333333 | 0.0333333333 |

A2's eps_x = -nu21 σ/E2 = -nu12 σ/E1 checks that the deck's Poisson convention matches the model
(the catalog's reciprocity check). The side cases (a part printed on its side) are 90°
permutations, still exact; they catch swapped axes, card order and stresses returned in local
axes. Arbitrary angles are TO-02.

(b) Shear, as simple shear. The catalog's displacement patch test needs prescribed non-zero
displacements, which Tessera does not have; this form uses a fixed face and face forces only, and
is still exact (the linear field meets every traction and support condition, forces and moments
balance). Other faces free, A = L².

| Case | Fixed face (all dofs) | Face forces (N) | Exact field | Tip u (mm) |
|---|---|---|---|---|
| S12 | y = 0 | y = L: (τA, 0, 0); x = L: (0, τA, 0); x = 0: (0, -τA, 0) | u = (γy, 0, 0), γ = τ/G12 = 0.01 | 0.1 |
| S13 | z = 0 | z = L: (τA, 0, 0); x = L: (0, 0, τA); x = 0: (0, 0, -τA) | u = (γz, 0, 0), γ = 0.0142857143 | 0.142857142857 |
| S23 | z = 0 | z = L: (0, τA, 0); y = L: (0, 0, τA); y = 0: (0, 0, -τA) | u = (0, γz, 0), γ = 0.0153846154 | 0.153846153846 |

τ_ij = 10 MPa at every node, the other components 0.

Precision: `test_reference_matches_formula` recomputes every table value from S to 1e-9. The
finite-element comparison is limited by the frd output, as measured in TS-03 (E12.5 rounding about
5e-6, nodal extrapolation bias 1e-5), so the proposed tolerances are TS-03's approved ones:
`rel_tolerance_stress = 2e-5`, `abs_tolerance_other_stress_mpa = 1e-6`,
`rel_tolerance_displacement = 1e-5`. Tightening to the catalog's 1e-6 needs ADR-008 first. About
27 solves of 1 to 8 s, solver tier only.

## Pull requests

1. **PR-0 (docs):** this ADR and the TO-01 catalog edit.
2. **PR-1 (model):** `elasticity.py`, `OrthotropicMaterial`, `MaterialFrame`, `Material`,
   `Section.orientation`, `SolveJob.materials: tuple[Material, ...]`, exports, tests.
3. After PR-1, in parallel worktrees (different packages, no file overlap):
   - **PR-2 (solvers):** deck writer, `DeckNames.orientations`, unit and ccx tests
     (`calculix-specialist`).
   - **PR-3 (printsetup):** `material_frame` and tests.
   - **PR-4 (model, schema 6):** `Project.materials`, `PrintSetup.raster_angle`, tests;
     `docs/user/project-files.md` via `docs-writer`.
   - **PR-5 (post + CLI):** `StudyResult.orthotropic` and the CLI note.
4. **PR-6 (wiring):** `materials.get_material`; `build_static_job(..., material: Material,
   frame: MaterialFrame | None = None)`; `core.study.solve` resolves material and frame and checks
   the print setup; `edit.add_material`, `edit.set_print_setup`. Needs PR-2, PR-3, PR-4.
5. **PR-7 (TO-01):** drafted in parallel once PR-4's API is merged (fails until PR-6); opened
   after PR-6 with the proposed toml in the description.
6. Optional: `tessera material add-orthotropic` and `tessera part print-setup` CLI commands; not
   needed for TO-01.

## Approvals needed

1. The Poisson and frame conventions (D1, D2) and the CLI wording (D4), with this ADR.
2. The TO-01 reference toml, including gating at the frd precision floor rather than 1e-6.
3. The made-up TO constants going into that reference file (the catalog already marks them as not
   user-trusted data).
4. No real filament data in this slice. PLA/PETG XY/Z values go through `materials-curator`, each
   number with a cited source and print parameters, in their own issue.

## Consequences

- Easier: TO-02 to TO-05 reuse a frame per section and region; failure criteria plug in through
  `rotate_stress` and the `orthotropic` flag.
- Harder: schema 6 must be read forever; project files carry user-entered stiffnesses, validated
  but not reviewed like library data; von Mises stays on screen for orthotropic parts until failure
  criteria land (the CLI note limits the risk).
- Revisit when slicer import lands: build direction and raster angle should come from the 3MF
  transform.

## Action items

1. [ ] Maintainer approves the conventions and CLI wording (approvals 1).
2. [ ] PR-1 model types.
3. [ ] PR-2 to PR-5 in parallel.
4. [ ] PR-6 wiring.
5. [ ] PR-7 TO-01 with the maintainer-created reference file.
6. [ ] ADR-004 item 3: mark orientations and orthotropic materials done after PR-6.
