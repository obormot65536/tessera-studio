# Tessera Studio documentation

Tessera Studio tells you whether a 3D-printed or machined part will hold under its loads. It
solves linear stress problems with finite elements (CalculiX) and reports displacement and
stress. The current release, v0.1, handles one part at a time.

## For makers: the user guide

- [Project files](user/project-files.md): the `.tessera` container and the `.tessera.json`
  form, how to create, pack, share and open them, and what schema versions mean.
- [Set up a study](user/study-setup.md): from a STEP or STL file to a solved result, step by
  step. Covers materials, face numbers, supports, loads, gravity, mesh size and how to read the
  output. Start here if you are new.
- [How we know the results are right](user/validation.md): the textbook checks Tessera passes,
  their tolerances, and what they do not show.

## For developers: the API

- [`tessera.core`](api/core.md): the public Python functions (`new_project`, `load_project`,
  `save_project`, `describe_faces`, `solve`), what they raise, and the model types they use.

## For contributors

- [Development cycle](development.md): issue, plan, test first, implement, review, merge.
- Architecture decision records in [`adr/`](adr/):
  - [ADR-001](adr/0001-package-layering.md): package layering
  - [ADR-002](adr/0002-units.md): internal unit system (mm, N, MPa, t, s)
  - [ADR-003](adr/0003-project-file-format.md): project file format
  - [ADR-004](adr/0004-solver-adapter-interface.md): solver adapter interface
  - [ADR-005](adr/0005-stl-import-and-faces.md): STL import and face definition
- Validation plan and status: [`benchmarks/CATALOG.md`](../benchmarks/CATALOG.md).
- Release notes: [`CHANGELOG.md`](../CHANGELOG.md).
