# ADR-011: Assemblies of placed parts joined by bonded connections (v0.5)

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
- Hard rule 6: the penalty stiffness, the gap a bond accepts, the flag rule at a bond and the
  TC-01 and TC-02 tolerances are numbers users trust.

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

### D3. Bonded connection → tied penalty interface

`calculix-specialist` measured the candidates in CalculiX 2.21 on two 10 mm steel cubes meshed
separately (sizes 3.0 and 2.1 mm, C3D10), with 10 MPa uniform tension. The main session
re-ran the `*TIE` and K = 1e10 cases and got the same numbers.

| Join | Tension across the interface: largest stress error | Largest displacement error | In-plane tension |
|---|---|---|---|
| `*TIE`, finer mesh as slave | 16 % | 0.8 % | exact |
| `*TIE`, coarser mesh as slave | 45 % | 2.6 % | exact |
| `*TIE`, both meshes halved | 16.5 % | 0.4 % | not run |
| Tied penalty, K = 1e8 MPa/mm | 6.6e-6 | 1.1e-4 | 1e-6 |
| Tied penalty, K = 1e10 MPa/mm | 1.0e-7 | 1.6e-6 | 1e-6 |
| Matching meshes, `*TIE` | 0 | 5e-7 | exact |

`*TIE` fails the patch test whenever load crosses the interface, and refining the mesh does not
fix it. Each slave node is placed exactly on the master face, but forces are not passed on
consistently. The 16 % error sits in the elements at the interface, so it would show up as a
false stress peak in every bonded part. For the TS-01 cantilever cut in two (TC-02) every
method is within 0.03 % of the single part, because the error stays local.

Decision:

- **Method.** A bonded connection becomes a surface-to-surface contact pair with tied
  behaviour, a penalty stiffness K, and a stick stiffness equal to K
  (`*SURFACE BEHAVIOR, PRESSURE-OVERCLOSURE=TIED`, `*FRICTION`, `*CONTACT PAIR, TYPE=SURFACE TO
  SURFACE`). CalculiX needs the friction line, or it stops with an error. The interface acts
  as glued in every direction, and the step converges in two iterations. A solve takes about
  twice as long as one with `*TIE`.
- **Penalty stiffness.** K = 1e5 × E_max / h, where E_max is the largest Young's modulus of the
  two parts' materials (any direction for an orthotropic one) and h is the smaller median edge
  length on the two bonded faces. For the steel cubes this gives K ≈ 1e10 MPa/mm, the row
  above. The error scales with E / (K h), so the rule keeps it near 1e-7 for any material and
  mesh. PR-3 checks that the matrix stays well conditioned with SPOOLES for PLA and for steel.
- **Solver-neutral type.** `tessera.model.solver.Bond(name, slave, master, stiffness)` names two
  mesh surfaces. `SolveJob.bonds: tuple[Bond, ...] = ()` is additive. The CalculiX adapter
  writes the contact pair. A job with bonds runs a static step that CalculiX solves
  iteratively, which the adapter allows for.
- **Surfaces.** Each bonded face becomes an element-face surface (`*SURFACE, TYPE=ELEMENT`, one
  `element, Sk` line per TRI6). The face id comes from which corners of its tetrahedron the
  triangle uses: {1,2,3} S1, {1,2,4} S2, {2,3,4} S3, {1,3,4} S4. The specialist checked this
  mapping with a uniform pressure. The adapter has no surfaces today (loads are nodal), so PR-3
  adds this.
- **Which side.** The finer mesh on the bonded faces (smaller median edge) is the slave. The
  penalty method passes either way (9e-6 with the coarse side as slave at K = 1e8), but the
  finer side is the safer choice.
- **Gap check.** Before the solve, `physics` measures how far each slave node is from the master
  surface. If any is farther than the bond tolerance, the solve is refused: "connection 'glue'
  does not touch: faces 6 of 'beam' and 1 of 'cap' are up to 0.42 mm apart (allowed 0.05 mm)".
  Proposed tolerance: 5 % of the smaller median edge length on the two faces, at least 1e-3 mm.
  This matters: with `*TIE` and a 0.02 mm gap, CalculiX left every node unjoined, still
  exited normally, and gave displacements 1e14 times too large. PR-3 finds out how the contact
  pair behaves at gaps below the tolerance, and adds a check of the solver output if it can
  also leave nodes unjoined.
- **Refused.**
  - Two bonds that share faces. They would overlap.
  - Contact, bolt and pin connections, with "contact connections arrive in a later v0.5
    release".
- **Model checks.** A bonded `Connection` with `friction` or `preload_n` set is refused. So is
  one whose two sides are the same part.
- **Solver.** SPOOLES, as now. PARDISO is not linked in this CalculiX build, and the iterative
  Cholesky solver converged slowly in the tests (84 iterations).
- **Later.** Two refinements are left for later, if solve time ever matters: merging nodes when
  two meshes match exactly, and `*TIE` for interfaces known to carry only in-plane load. Mortar
  contact (`TYPE=MORTAR`) crashed with no output in CalculiX 2.21 and is not used.

### D4. Naming spots in results

- `StressSpot` and `StrengthSpot` gain `part: str`. Face peaks and strength face minima are
  keyed by (part, face).
- CLI wording: with one part, every line stays as it is today. With several parts, the face
  becomes "face 3 of 'cap'" and the per-face lists are grouped by part. Docs and tests for
  single parts do not change.
- `MeshCheck` matches spots and faces by (part, face).

### D5. Trust at a bond

- **Sharp inside edges** are found on the assembly's outer surface. Bonded faces are inside the
  assembly and are left out. Without this, the inside corner where a rib meets a plate is
  invisible: each part on its own has no inside corner there.
- **Bond edges.** A new flag kind, `FlagKind.BOND_EDGE`, flags nodes within the flag distance of
  the boundary of each bonded face pair. Where a bond ends, the stress depends on the mesh, as
  at a support, so these spots never become the headline unless the two-mesh check settles
  them. The line reads "near the edge of bond 'glue'".
- **Inside a bond**, away from its edge, nothing is flagged. TC-01 shows the bond carries a
  uniform stress across non-matching meshes to about 1e-7 (D3).

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

Both meshes refine every part by `MESH_CHECK_REFINEMENT`. The gap check, slave choice and K
are redone on each mesh.

## Options considered

### Joining the parts

| Option | For | Against |
|---|---|---|
| A. `*TIE` | Linear and cheap. Meshes need not match. Standard CalculiX | Fails the patch test when load crosses the interface: 16 % stress error, a false peak in every bonded part (D3) |
| B. Merge coincident nodes | Exact when meshes match | Needs matching meshes, which gmsh does not give for separately meshed parts |
| C. Tied penalty contact pair (chosen) | Passes the patch test to 1e-7. Either side can be the slave. Same input as later contact | Iterative step, about twice as slow as `*TIE`. A penalty stiffness to choose and approve |
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

The tied penalty pair keeps each part's mesh independent, which keeps STEP and STL on one path
and face numbers stable, and it passes the patch test. `*TIE` would be faster, but its false
16 % peak at every bond is exactly the kind of number Tessera exists to get right. The price
is a penalty stiffness (a number users trust, set by a rule rather than by the user), an
iterative step about twice as slow as a tie, and a new kind of flagged spot at a bond's edge. A
conforming mesh (D) would remove the penalty and the flags, but would tie every part's meshing
to every other's and could not handle STL.

## Test-first plan

Unit tests (fail before the code):
- `model`: `Transform` refuses a non-rotation and a mirror. A bonded connection with friction or
  preload is refused, and so is one with the same part on both sides. `Bond` and
  `SolveJob.bonds` default to empty. The spot `part` field.
- `physics`: merging two meshes gives offset nodes, prefixed sets and one section per part.
  Placement moves nodes and turns the material frame. The gap check refuses faces 0.1 mm apart
  and accepts touching ones. The slave is the finer side, and sides are swapped away from a
  supported face. A part with no path to a support is refused.
- `solvers`: the deck writes the tied contact pair with element-face surfaces (face ids checked
  against the table in D3), and a job without bonds gives the same deck as today.
- `post`: bonded faces are not outer surface. `BOND_EDGE` flags the bond boundary and not its
  middle. Face peaks are keyed by part.
- `core`: placement invariance. TS-01 turned 90° about z and moved gives the same tip deflection,
  turned. An orthotropic cube turned so that axis 1 lies along global x gives the TO-01
  `A1-side` result. Contact, bolt and pin connections are refused with the message.
- `cli`: `part place`, `connect bond`, `info` lines, and two-part result wording.

Benchmarks (`fea-verifier` writes them first; references proposed in the PR and created by the
maintainer):
- **TC-01**, bonded patch test with non-matching meshes: z tension and x tension. The catalog's
  1e-4 relative is now realistic: D3 measured 1e-7 (stress) and 1.6e-6 (displacement) with the
  proposed K. The catalog note that ties are not patch-test-consistent is replaced with the D3
  numbers.
- **TC-02**, the TS-01 cantilever cut at mid-span and bonded (measured: -0.40 %, the same as
  the single part): within 2 % of 1.526724 mm and within
  0.5 % of the single-part result at the same mesh size.

## Pull requests

1. **PR-0 (docs):** this ADR, and the TC-01 and TC-02 catalog entries updated with the D3
   measurements.
2. **PR-1 (model):** `Transform` and `Connection` validators, `Bond`, `SolveJob.bonds`, spot
   `part`, `FlagKind.BOND_EDGE`.
3. **PR-2 (physics):** merge, placement, frame rotation, bonded → `Bond` with K, gap check,
   slave choice, held-part check.
4. **PR-3 (solvers):** element-face surfaces and the tied contact pair in the deck, behaviour at
   small gaps, conditioning (`calculix-specialist`).
5. **PR-4 (post):** keys by part, assembly outer surface, bond-edge flags.
6. **PR-5 (core + edit):** multi-part solve, `place_part`, `add_bond`, the refusal for other
   connection kinds.
7. **PR-6 (CLI and docs):** commands and wording, and the user guide section on assemblies.
8. **PR-7 (benchmarks):** TC-01 and TC-02 with proposed reference files.

PR-2 to PR-4 can run in parallel after PR-1. PR-3 and PR-7's decks can be tried early against
hand-built jobs.

## Approvals needed

1. Tied penalty contact pairs for bonded connections, not `*TIE`, with the slave rule and the
   refusals (D3).
2. The penalty stiffness rule K = 1e5 × E_max / h (D3).
3. The bond tolerance (5 % of the smaller median edge, at least 1e-3 mm) and refusing the solve
   when a gap exceeds it (D3).
4. The `BOND_EDGE` flag and leaving bonded faces out of the outer surface (D5).
5. Placement rules: mirrors refused, the print setup in part axes (D1).
6. Results and CLI wording for several parts (D4, D6).
7. The TC-01 and TC-02 tolerances (in PR-0 and PR-7).

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
3. [ ] ADR-004 item 3: record bonded connections (ties) as done after PR-5.
4. [ ] Tick ADR-009's action items, which shipped with TO-01 (#43).
