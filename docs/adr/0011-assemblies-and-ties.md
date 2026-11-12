# ADR-011: Assemblies of placed parts joined by bonded (tied) connections (v0.5)

**Status:** Proposed (draft for the maintainer, issue #69)
**Deciders:** Maintainer

## Context

The v0.5 gate asks for Hertz and bolted-joint tests (TC-04, TC-07, TC-08). Every one of them has
two or more parts. This ADR covers the first step: several parts, each placed in space, bonded
where they touch. Contact (TC-03, TC-04) and bolts and pins (TC-06 to TC-09) get their own ADRs
and build on this one.

What exists today:

- **Model.** `Assembly.parts` holds `PartInstance`s, each with a `Transform` (row-major rotation,
  then translation in mm). `Assembly.connections` holds `Connection(id, kind, a, b, friction,
  preload_n)` with `kind` bonded / contact / bolt / pin and two `FaceSet`s. Schema 7 stores
  both, and `Assembly` checks that connection part ids exist. Supports and loads already name
  their part (`FaceSet.part_id`, `EdgeSet.part_id`).
- **Solve.** `tessera.core.solve` refuses more than one part, any connection, and any
  non-identity placement ("v0.1 solves a single part; assemblies arrive in v0.5").
- **Mesh.** Each part is meshed on its own (`mesh_step`, `mesh_surface`). Its `Mesh` names its
  sets `part:<id>` (elements) and `face:<i>` (TRI6 surfaces and their node sets).
- **Job.** `SolveJob` already takes several materials and sections. `physics.build_static_job`
  builds the job of one part. Supports, forces, moments and gravity are in global axes.
- **Results.** `StressSpot.face` and `StrengthSpot.face` give a face number only. Face peaks
  and strength face minima are keyed by face number. ADR-007 flags nodes near singular supports
  and sharp inside edges, found on the part's outer surface.
- **Orientation.** The material frame comes from the part's print setup, in part axes
  (ADR-009 D2).
- **Edit and CLI.** No function or command places a part or adds a connection.

Constraints:
- ADR-001: domain packages do not import each other, and core orchestrates.
- ADR-002: mm, N, MPa.
- ADR-004: no CalculiX names in `tessera.model`, and changes are additive.
- ADR-007 decision 5: flagged spots never become the headline unless settled.
- Hard rule 6: the tie tolerance, the gap a bond accepts, the flag rule at a tie and the TC-01
  and TC-02 tolerances are numbers users trust.

## Decision

### D1. Placement

- Each part is meshed in its own axes, as today. Core then moves its nodes to global axes:
  x_global = R x_part + t. Faces, edges and sets keep their part-local numbering.
- Supports, forces, moments and gravity stay in global axes, as now. A support's `dofs` and a
  force vector mean the same thing for a placed part as for one at the origin.
- The material frame turns with the part. The print setup stays in part axes, because it
  describes how that part was printed. The frame used for the solve is R applied to the
  ADR-009 frame (rows: axes_global = axes_part @ R^T).
- `Transform` gets a validator: R must be a rotation (R R^T = I and det R = +1, to 1e-9). A
  mirror is refused, because it would flip the part's face normals and handedness. The schema
  number does not change: the edit API never wrote another value, so only a hand-edited file
  can fail, and it fails with a clear message.

### D2. One job from several meshes

- `physics` merges the parts' meshes into one `Mesh` with node and element offsets.
- Set names gain the part id: `part:<id>` (unchanged), `face:<id>:<i>` and `surface` names
  likewise. A single-part study uses the same names, so there is one code path.
- Nodes of different parts are never merged, even where they coincide. Parts are joined only by
  connections. Coincident nodes without a connection mean the parts are not joined.
- Each part gets its own `Section` (material and frame). Gravity uses each part's density.
- A part that no support reaches through bonded connections is refused before the solve:
  "part 'cap' is not held: add a support on it or bond it to a held part". Without it CalculiX
  fails with a singular matrix, and the message says nothing useful.

### D3. Bonded connection → tie

*(Waiting for the `calculix-specialist` measurements. The choices below are the draft
proposal.)*

- **Solver-neutral type.** `tessera.model.solver.Tie(name, slave, master, tolerance)` names two
  mesh surfaces. `SolveJob.ties: tuple[Tie, ...] = ()` is additive. The CalculiX adapter
  writes it as `*TIE`.
- **Which side.** The slave is the side with the finer mesh on the bonded faces (smaller median
  edge length). On a tie the slave nodes follow the master surface, so the finer side as slave
  avoids gaps opening between coarse master nodes.
- **Tolerance and gap check.** Before the solve, `physics` measures how far each slave node is
  from the master surface. If any node is farther than the bond tolerance, the solve is refused:
  "connection 'glue' does not touch: faces 6 of 'beam' and 1 of 'cap' are up to 0.42 mm apart
  (allowed 0.05 mm)". The same tolerance goes to CalculiX, so no node is silently left
  untied. Proposed tolerance: 5 % of the smaller median edge length on the two faces, at least
  1e-3 mm.
- **Refused.**
  - A node that is a slave in two ties. The connections overlap.
  - A slave node that is also held by a support. Tessera swaps the sides for that tie if the
    other side is free, and refuses otherwise.
  - Contact, bolt and pin connections, with "contact connections arrive in a later v0.5
    release".
- **Model checks.** A bonded `Connection` with `friction` or `preload_n` set is refused. So is
  one whose two sides are the same part.

### D4. Naming spots in results

- `StressSpot` and `StrengthSpot` gain `part: str`. Face peaks and strength face minima are
  keyed by (part, face).
- CLI wording: with one part, every line stays as it is today. With several parts, the face
  becomes "face 3 of 'cap'" and the per-face lists are grouped by part. Docs and tests for
  single parts do not change.
- `MeshCheck` matches spots and faces by (part, face).

### D5. Trust at a tie

- **Sharp inside edges** are found on the assembly's outer surface. Bonded faces are inside the
  assembly and are left out. Without this, the inside corner where a rib meets a plate is
  invisible: each part on its own has no inside corner there.
- **Tie edges.** A new flag kind, `FlagKind.TIE_EDGE`, flags nodes within the flag distance of
  the boundary of each bonded face pair. Where a bond ends, the stress depends on the mesh, as
  at a support, so these spots never become the headline unless the two-mesh check settles
  them. The line reads "near the edge of bond 'glue'".
- **Inside a bond**, away from its edge, nothing is flagged. TC-01 shows the tie carries a
  uniform stress across non-matching meshes to its tolerance.

### D6. Edit API and CLI

- `edit.place_part(project, part_id, *, translation=(0, 0, 0), rotations=())`. Rotations are
  (axis, degrees) pairs about the global axes, applied in order, then the translation. It
  replaces the whole placement.
- `edit.add_bond(project, part_a, faces_a, part_b, faces_b, *, connection_id=None)` checks that
  the faces exist and the parts differ. Whether the faces touch is checked at solve time, when
  the meshes exist.
- CLI:
  - `tessera part place FILE PART --move X Y Z --turn z 90 [--turn x 90]`
  - `tessera connect bond FILE PART_A FACES... --to PART_B --to-faces FACES... [--id ID]`
  - `tessera info` lists placements and connections. `tessera remove` learns to remove a
    connection; it already refuses to remove a part that a connection names.
- A command that lists which faces of different parts touch would help makers find the face
  numbers. It is left out of this ADR and noted as a follow-up.

### D7. Mesh check

Both meshes refine every part by `MESH_CHECK_REFINEMENT`. The gap check and slave choice are
redone on each mesh.

## Options considered

### Joining the parts

| Option | For | Against |
|---|---|---|
| A. `*TIE` (chosen, pending D3 measurements) | Linear and cheap. Meshes need not match. Standard CalculiX | Small error at non-matching quadratic faces (TC-01 measures it) |
| B. Merge coincident nodes | Exact when meshes match | Needs matching meshes, which gmsh does not give for separately meshed parts |
| C. `*CONTACT PAIR` with tied behaviour | Same input as later contact | Nonlinear solve and penalty stiffness for a linear bond. More settings users must trust |
| D. One gmsh model with shared faces (fragment) | Conforming mesh, exact | Couples the meshing of all parts. STL parts cannot be fragmented. Changes the face numbering users already see |

### Where placement is applied

Meshing in part axes and then moving the nodes (chosen) keeps face numbers, meshes and STL
handling per part, as now. Moving the geometry before meshing would need gmsh transforms for
STEP and a separate path for STL.

### Results naming

Always printing the part (consistent, but changes every single-part line and its docs) versus
printing it only with several parts (chosen: no churn for single parts). The maintainer may
prefer the first.

## Trade-off analysis

Ties keep the analysis linear and each part's mesh independent, which keeps STEP and STL on one
path and face numbers stable. The price is a small interface error and a new kind of flagged
spot at the edges of a bond. A conforming mesh (D) would remove both but would tie every part's
meshing to every other's, and could not handle STL. That cost grows with each part. Contact
(C) is needed later anyway, but bringing it in for a bond would make the first assembly release
nonlinear for no gain.

## Test-first plan

Unit tests (fail before the code):
- `model`: `Transform` refuses a non-rotation and a mirror. A bonded connection with friction or
  preload is refused, and so is one with the same part on both sides. `Tie` and
  `SolveJob.ties` default to empty. The spot `part` field.
- `physics`: merging two meshes gives offset nodes, prefixed sets and one section per part.
  Placement moves nodes and turns the material frame. The gap check refuses faces 0.1 mm apart
  and accepts touching ones. The slave is the finer side, and sides are swapped away from a
  supported face. A part with no path to a support is refused.
- `solvers`: the deck writes `*TIE` with element-face surfaces, and a job without ties gives the
  same deck as today.
- `post`: bonded faces are not outer surface. `TIE_EDGE` flags the bond boundary and not its
  middle. Face peaks are keyed by part.
- `core`: placement invariance. TS-01 turned 90° about z and moved gives the same tip deflection,
  turned. An orthotropic cube turned so that axis 1 lies along global x gives the TO-01
  `A1-side` result. Contact, bolt and pin connections are refused with the message.
- `cli`: `part place`, `connect bond`, `info` lines, and two-part result wording.

Benchmarks (`fea-verifier` writes them first; references proposed in the PR and created by the
maintainer):
- **TC-01**, tied patch test with non-matching meshes: z tension and x tension. The tolerance
  comes from the `calculix-specialist` measurement (D3), replacing the catalog's placeholder
  1e-4.
- **TC-02**, the TS-01 cantilever cut at mid-span and tied: within 2 % of 1.526724 mm and within
  0.5 % of the single-part result at the same mesh size.

## Pull requests

1. **PR-0 (docs):** this ADR, and the TC-01 and TC-02 catalog entries updated with the measured
   tolerance.
2. **PR-1 (model):** `Transform` and `Connection` validators, `Tie`, `SolveJob.ties`, spot
   `part`, `FlagKind.TIE_EDGE`.
3. **PR-2 (physics):** merge, placement, frame rotation, bonded → `Tie`, gap check, slave
   choice, held-part check.
4. **PR-3 (solvers):** `*TIE` in the deck (`calculix-specialist`).
5. **PR-4 (post):** keys by part, assembly outer surface, tie-edge flags.
6. **PR-5 (core + edit):** multi-part solve, `place_part`, `add_bond`, the refusal for other
   connection kinds.
7. **PR-6 (CLI and docs):** commands and wording, and the user guide section on assemblies.
8. **PR-7 (benchmarks):** TC-01 and TC-02 with proposed reference files.

PR-2 to PR-4 can run in parallel after PR-1. PR-3 and PR-7's decks can be tried early against
hand-built jobs.

## Approvals needed

1. Ties (`*TIE`) for bonded connections, with the slave rule and the refusals (D3).
2. The bond tolerance (5 % of the smaller median edge, at least 1e-3 mm) and refusing the solve
   when a gap exceeds it (D3).
3. The `TIE_EDGE` flag and leaving bonded faces out of the outer surface (D5).
4. Placement rules: mirrors refused, the print setup in part axes (D1).
5. Results and CLI wording for several parts (D4, D6).
6. The TC-01 and TC-02 tolerances (in PR-0 and PR-7).

## Consequences

- Easier: makers can check glued or welded parts, and parts printed separately and bonded. The
  merge and naming work is what contact and bolts need too.
- Harder: each result needs a part as well as a face. A bond's edge adds flagged spots that the
  two-mesh check must settle. The gap check makes a solve fail where a slightly misplaced part
  would once have been silently loose.
- Later: contact ADR (TC-03, TC-04), bolts and pins ADR (TC-06 to TC-09), a "touching faces"
  helper, and the GUI.

## Action items

1. [ ] Maintainer approves the items above.
2. [ ] PR-0 to PR-7.
3. [ ] ADR-004 item 3: record ties as done after PR-5.
4. [ ] Tick ADR-009's action items, which shipped with TO-01 (#43).
