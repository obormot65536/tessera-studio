# Changelog

## Unreleased

- Project files are now schema 7 (ADR-010 D3): `IsotropicDefinition` and
  `OrthotropicDefinition` take an optional `strength` (`IsotropicStrengthDefinition`: `tensile`,
  `compressive`, `brittle`; `OrthotropicStrengthDefinition`: `t1`-`t3`, `s12`, `s13`, `s23`,
  optional `c1`-`c3`; MPa, positive), and `IsotropicMaterial` / `OrthotropicMaterial` take
  `strength=` (`IsotropicStrength`, `OrthotropicStrength`, keyword only; a strength of the other
  kind is refused). A missing compressive strength is taken equal to the tensile one and named
  in `assumed`; a ductile material refuses a compressive strength. The solver ignores strengths
  (decks unchanged). Result types for the coming strength check: `FailureCriterion`,
  `FailureMode`, `StrengthSpot`, `StrengthResult`, `StudyResult.strength` and
  `MeshCheck.lowest_strength_ratio` (both `None` until it lands). Schema 6 files open unchanged;
  older Tessera versions refuse schema 7 files.
- Solve parts made of your own orthotropic (layer-wise) material (ADR-009): define a material
  in the project with `edit.add_material` (an `IsotropicDefinition` or
  `OrthotropicDefinition`, or its fields; an id that is built in or already defined is
  refused), give the part a print setup with `edit.set_print_setup(project, part_id, process,
  build_direction, raster_angle)`, and `solve` orients the material axes from the build
  direction and raster angle. An orthotropic part without a print setup is refused ("needs a
  print setup (build direction)") rather than solved in the part axes; an isotropic part
  ignores its print setup. `set_material`, `add_part` and the gravity density check accept
  project materials. API: `tessera.materials.get_material(material_id, defined=())` (project
  definitions win over the library) and `build_static_job(..., material, frame=None)`.
  Invalid orthotropic constants on these paths raise `InputError`. No command-line support
  yet.
- Material frame from a print setup (ADR-009 D2):
  `tessera.printsetup.material_frame(build_direction, raster_angle=0.0)` gives the `MaterialFrame` of a part printed with that direction up (axis 3)
  and rasters turned by `raster_angle` degrees about it (axis 1 starts at printer X carried by
  the smallest tilt; upside down uses a half turn about X).
- Project files are now schema 6 (ADR-009 D3): `Project.materials` holds materials defined in
  the project (`IsotropicDefinition`, `OrthotropicDefinition`, discriminated on `kind`, ids
  unique, orthotropic constants checked with `check_orthotropic`), and `PrintSetup` gains
  `raster_angle` (degrees, default 0). Schema 5 files open unchanged; older Tessera versions
  refuse schema 6 files.
- Orthotropic results are labelled (ADR-009 D4): `StudyResult.orthotropic` is true when any
  material in the study is an `OrthotropicMaterial` (set by `tessera.post.summarize`), and
  `tessera solve` then prints a note before the stress lines that von Mises shows where the
  stress is but is not a strength check for a layer-wise part. Peaks, flags and mesh-check
  verdicts are unchanged.
- Thin-wall guard for the two-mesh check (ADR-007 Amendment 2): a marked peak that changed by
  5 % or less is no longer called `settled` when the coarse mesh had fewer than 5 elements across
  the wall there, too few to show the stress rise, or when the wall could not be measured; it
  reads `not judged` with the elements across and a mesh size that gives 5. The wall is the
  width of the largest ball inside the part touching the surface at the peak, never doubled
  across a support; at a peak marked by a face support, that support's faces are left out, as
  the part continues into the fixture. Without `--check-mesh`, a marked peak on such a wall gets a `thin wall` line
  with the same advice, and marked face peaks say `thin wall`. Unmarked numbers are judged as
  before. API: `StressSpot.wall` (`WallSpan`), `Verdict.THIN_WALL`, `SpotCheck.elements_across`,
  `MeshCheck.wall_elements_min`; `WALL_ELEMENTS_MIN`, `outer_surface`, `wall_span`,
  `samples_on` and `node_edge_length` in `tessera.post`, `WALL_ELEMENTS_MIN` in `tessera.core`.
- Orthotropic material types (ADR-009, model only): `OrthotropicMaterial` (engineering
  constants, major Poisson ratios, validated to a positive definite compliance), `Material`,
  `MaterialFrame` and `Section.orientation` in `tessera.model`, and
  `tessera.model.elasticity` (`orthotropic_compliance`, `check_orthotropic`, `rotate_stress`).
  `SolveJob.materials` now takes either kind.
- CalculiX deck for orthotropic materials (ADR-009): `write_deck` writes an
  `OrthotropicMaterial` as `*ELASTIC, TYPE=ENGINEERING CONSTANTS` and a section's
  `MaterialFrame` as `*ORIENTATION, SYSTEM=RECTANGULAR` with `ORIENTATION=` on its
  `*SOLID SECTION` (names in `DeckNames.orientations` and `ccx_name.orientation.<elset>`
  metadata). Results stay in global axes. Decks without orientations are unchanged.
- API (ADR-008): `SolveJob(point_stresses=True)` asks for each element's own stresses at its
  integration points, returned in `SolveResult.point_stresses` (`PointStresses`). Opt-in;
  no reported number changes.
- CalculiX fills it from `job.dat` when asked; a missing, empty or incomplete `job.dat` is
  then a `SolverError`.
- `tessera.post` (ADR-008): `element_nodal_stresses`, `point_von_mises`, `flag_points`;
  `von_mises` keeps the leading axes, so (M, Q, 6) gives (M, Q).
- Solver timeout (ADR-004): a solve started through `tessera.core` now stops CalculiX after
  1800 s per solver run instead of waiting forever, together with the processes it started
  (except one that detaches itself with `setsid`). Ctrl+C, and in `tessera solve` on Linux
  and macOS closing the terminal or SIGTERM, also stop it; if Tessera itself is killed, only
  Linux stops ccx. `tessera solve --timeout SECONDS` and `solve(project, timeout=...)` change the limit; with
  `--check-mesh` each of the two solves gets the full limit, and the message names the one
  that was stopped. A timeout is a `SolverTimeoutError` (a `SolverError`, in
  `tessera.model`) and exit code 2 with `solver failed: CalculiX did not finish within ... s
  and was stopped`. `tessera.core.SOLVE_TIMEOUT` holds the default. `CalculixAdapter` itself
  still defaults to no limit, and refuses a timeout that is not a positive number.
- Two-mesh check (ADR-007): `tessera solve --check-mesh` solves the study again with every
  element about 1/√2 as long and marks each printed number `settled` (changed by 5 % or less),
  `mesh-dependent`, or `not judged` where the second mesh did not get finer (its elements there
  more than 0.85 times their first size). All numbers come from the finer mesh; it takes 3 to 5
  times as long. API: `solve(project, check_mesh=True)` fills `StudyResult.mesh_check`
  (`MeshCheck`, `SpotCheck`, `Verdict` in `tessera.model`; `compare_meshes` in
  `tessera.post`); `mesh_step` and `mesh_surface` take a `size_factor`.
- Where the peak stress is, and whether it depends on the mesh (ADR-007): `tessera solve` now
  leads with the largest von Mises stress away from supports and sharp inside edges, then the
  overall peak with its location, face and, when it sits near a fixed support or a sharp
  inside edge, a note that it may keep rising as the mesh is refined, then the highest stress on
  the top five faces. The API adds `StudyResult.peak_von_mises`, `clear_peak_von_mises`,
  `face_peaks` and `node_flags` (`StressSpot`, `StressFlag`, `NodeFlags` in `tessera.model`).
  `max_von_mises` is unchanged. Where the highest stress is at a support (a clamped base under
  bending), the value away from supports reads low; the user guide explains how to read both.
- Line supports (ADR-006): hold a part along an edge, named by the two faces that meet there,
  so it can still turn about the line (a pinned or roller line, e.g. a simply supported beam).
  CLI: `tessera support add PATH PART --edge 3-5 [--edge 4-5] --dofs z`, and `tessera edges
  PATH PART [--face N]` lists a part's edges; `tessera info` shows them. API:
  `edit.add_edge_support`, `describe_edges`, `EdgeSet`, `EdgeInfo`, `Support.edges` and
  `Support.part_id`. Peak stress right at a line support depends on the mesh.
- Project files are now schema 5. Schema 4 files open unchanged; older Tessera versions refuse
  schema 5 files.
- A force, pressure or moment of zero is refused when it is added, since it does nothing and is
  almost always a typo. Project files that already contain one still open.
- Gravity and materials without a density are refused when you set them up (turning gravity
  on, or choosing such a material while gravity is on), not only when you solve.
- Faster study setup: a part's geometry is read once, not again for every support and load
  added to it (it is read again if the file changes).
- CLI: set up a study from the command line with `tessera part add`, `part material`, `support add`, `load force|pressure|moment`, `study`, `remove` and `materials`; `tessera info` summarises the project.
- Clearer errors: a support or load on a part that does not exist is refused instead of
  being applied to the only part, and a missing geometry file names the part.
- API: `save_project` raises `InputError` for an invalid project (it raised pydantic's
  `ValidationError`; both are `ValueError`s).
- Moment loads: apply a torque (N mm) to a face. It is spread over the face as forces that
  add no net push or pull (`MomentLoad`).
- Pressure loads (MPa, normal to the face, positive pushes into the part; works on curved
  faces) and gravity (mm/s²; uses the material density).
- `.tessera` project files: one file that holds the project and a copy of its STEP or STL
  geometry, so it still opens after the original files move. `.tessera.json` stays the
  hand-editable form.
- `tessera pack` bundles a `.tessera.json` project and its geometry into a `.tessera` file.
- Project files carry a schema version (now 4). Older files still open; files from a newer
  Tessera are refused with a message that names both versions.
- CLI: `tessera faces` lists a part's faces for choosing supports and loads; `tessera solve`
  meshes and solves a one-part study and prints the largest displacement and von Mises stress.
- Linear static analysis of one part with CalculiX (needs `ccx` on the PATH, or set
  `TESSERA_CCX` to its location).
- Built-in materials: `steel-generic` and `aluminium-6061-t6`, both isotropic.
- Validation: eight textbook cases now pass (patch test, cantilever, plate with a hole, bar under
  its own weight, thick-walled cylinder, thick-walled sphere, clamped circular plate, torsion); see
  `docs/user/validation.md`.
- STL import (binary and ASCII): faces are split at edges sharper than 40°, the file must be a
  single closed solid, and lengths are taken as mm (ADR-005).
- Roller supports: a support can hold only chosen directions (x, y or z), for symmetry planes
  and sliding faces.
- Python API in `tessera.core` (`new_project`, `load_project`, `save_project`,
  `describe_faces`, `solve`); see `docs/api/core.md`.
- `tessera` CLI commands `new`, `info` and `version`. The model can describe assemblies, but
  assemblies cannot be solved yet.
