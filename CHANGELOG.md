# Changelog

## Unreleased

- Orthotropic results are labelled (ADR-009 D4): `StudyResult.orthotropic` is true when any
  material in the study is an `OrthotropicMaterial` (set by `tessera.post.summarize`), and
  `tessera solve` then prints a note before the stress lines that von Mises shows where the
  stress is but is not a strength check for a layer-wise part. Peaks, flags and mesh-check
  verdicts are unchanged. Projects cannot define an orthotropic material yet, so no solve shows
  the note until the wiring lands.
- Orthotropic material types (ADR-009, model only): `OrthotropicMaterial` (engineering
  constants, major Poisson ratios, validated to a positive definite compliance), `Material`,
  `MaterialFrame` and `Section.orientation` in `tessera.model`, and
  `tessera.model.elasticity` (`orthotropic_compliance`, `check_orthotropic`, `rotate_stress`).
  `SolveJob.materials` now takes either kind. Nothing uses them yet: projects cannot define an
  orthotropic material, and the CalculiX adapter refuses one with `NotImplementedError`.
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
