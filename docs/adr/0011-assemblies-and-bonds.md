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
- **Solve.** `tessera.core.solve` refuses more than one part or any connection ("v0.1 solves a
  single part; assemblies arrive in v0.5"), and any non-identity placement ("v0.1 does not
  apply part placements yet"). It asks for point stresses only when the part's material has a
  strength.
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

- Each part is meshed in its own axes, as today. `physics`, which also merges the meshes (D2),
  then moves its nodes to global axes: x_global = R x_part + t. Faces, edges and sets keep
  their part-local numbering.
- Supports, forces, moments and gravity stay in global axes, as now. A support's `dofs` and a
  force vector mean the same thing for a placed part as for one at the origin.
- The material frame turns with the part. The print setup stays in part axes, because it
  describes how that part was printed. The frame used for the solve is R applied to the
  ADR-009 frame (rows: axes_global = axes_part @ R^T). `MaterialFrame` then holds global axes,
  so PR-1 updates its docstring ("rows = material axes in part coordinates") and PR-0 adds a
  note to ADR-009 D2.
- `Transform` gets a validator: R must be a rotation (R R^T = I and det R = +1, to 1e-9). A
  mirror is refused, because it would flip the part's face normals and handedness. The schema
  number does not change: the edit API never wrote another value, so only a hand-edited file
  can fail, and it fails with a clear message.

### D2. One job from several meshes

- `physics` merges the parts' meshes into one `Mesh` with node and element offsets.
- Set names gain the part id: `part:<id>` stays, and `face:<i>` (both the TRI6 surface and its
  node set) becomes `face:<id>:<i>`. A single-part study uses the same names, so there is one
  code path.
- Nodes of different parts are never merged, even where they coincide. Parts are joined only by
  connections. Coincident nodes without a connection mean the parts are not joined.
- Each part gets its own `Section` (material and frame). Gravity uses each part's density, and
  refuses as today if any part's material has none.
- Point stresses are asked for when any part's material has a strength. The strength check
  already works per section (ADR-010), so a part without strengths is listed as unchecked
  while the others are checked.
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
  two parts' materials (for an orthotropic one, max(E1, E2, E3)) and h is the smaller median
  edge length on the two bonded faces. For the steel cubes this gives K ≈ 1e10 MPa/mm, the
  row above.
- **What error to expect.** The penalty lets the faces separate by σ / K. Over a part of length
  L, the relative displacement error is about E / (K L) = 1e-5 × h / L, so at most about 1e-5
  and smaller for long parts. The 1e-7 stress error was measured on a uniform patch of one
  material. Every measurement so far bonds steel to steel, so PR-3 and TC-01 add PLA bonded to
  steel and an orthotropic part before the rule is relied on for other materials. PR-3 also
  checks that the matrix stays well conditioned with SPOOLES.
- **Solver-neutral type.** `tessera.model.solver.Bond(name, slave, master, stiffness)` names two
  mesh surfaces. `SolveJob.bonds: tuple[Bond, ...] = ()` is additive. The CalculiX adapter
  writes the contact pair. A job with bonds runs a static step that CalculiX solves
  iteratively, which the adapter allows for.
- **Surfaces.** Each bonded face becomes an element-face surface (`*SURFACE, TYPE=ELEMENT`, one
  `element, Sk` line per TRI6). The face id comes from which corners of its tetrahedron the
  triangle uses: {1,2,3} S1, {1,2,4} S2, {2,3,4} S3, {1,3,4} S4. The specialist checked this
  mapping with a uniform pressure. The adapter has no surfaces today (loads are nodal), so PR-3
  adds this.
- **Faces that only partly cover each other.** Bonded faces often differ in size: a rib's end
  face sits on part of a plate's face. Before the solve, `physics` measures how far each node of
  one side is from the other side's surface, both ways. The side whose nodes all lie within the
  bond tolerance of the other is the **covered** side, and the bond is its whole face set.
  - If only one side is covered (the rib), it is the slave, whatever the mesh sizes. The
    plate's nodes outside the rib's footprint are simply not bonded.
  - If both are (same footprint), the finer mesh (smaller median edge) is the slave. The
    penalty method passes either way (9e-6 with the coarse side as slave at K = 1e8), but the
    finer side is the safer choice.
  - If neither is, the faces are apart or only partly overlap, and the solve is refused. Where
    any node of the covered candidate is far off, the message gives the gap: "connection 'glue'
    does not touch: faces 6 of 'beam' and 1 of 'cap' are up to 0.42 mm apart (allowed
    0.05 mm)". For two faces that overlap only in part (a staggered lap), it says so and asks
    for the faces to be split in CAD so that one covers the other. Staggered overlaps are a
    v0.5 limit.
  - Proposed bond tolerance: 5 % of the smaller median edge length on the two faces, at least
    1e-3 mm.
  This matters: with `*TIE` and a 0.02 mm gap, CalculiX left every node unjoined, still
  exited normally, and gave displacements 1e14 times too large. PR-3 finds out how the contact
  pair behaves at gaps below the tolerance, and adds a check of the solver output if it can
  also leave nodes unjoined.
- **Refused.**
  - Two bonds that share faces. They would overlap.
  - Nothing else about supports: a supported node on a bonded face is not a known problem for
    the penalty pair (the conflict the specialist found was in `*TIE`'s equations). PR-3 checks
    it, and this ADR is amended if it turns out to be one.
  - Contact, bolt and pin connections, with "contact connections arrive in a later v0.5
    release".
- **Model checks.** A bonded `Connection` with `friction` or `preload_n` set is refused. So is
  one whose two sides are the same part. As for `Transform` (D1), the schema number does not
  change: no edit function has ever written a connection, so only a hand-edited file can
  fail, with a clear message.
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

The **bonded region** of a connection is its covered face set and, on the other side, the
triangles whose centroids lie on it (within the bond tolerance). Its **boundary** is the
boundary of the covered face set.

- **Not a free surface.** The bonded region is left out of the outer surface used for sharp
  inside edges and for the wall span (ADR-007 Amendment 2): the wall continues into the other
  part there. The rest of each face, such as the plate around a rib, stays outer surface.
- **Bond edges.** A new flag kind, `FlagKind.BOND_EDGE`, flags nodes within the flag distance of
  each bond's boundary. This is also where the inside corner of a rib on a plate lies, which
  neither part shows on its own, so this flag covers it. Where a bond ends the stress depends
  on the mesh, as at a support, so these spots never become the headline unless the two-mesh
  check settles them. The line reads "near the edge of bond 'glue'".
- **Inside a bond**, away from its boundary, nothing is flagged. This is a trust decision (it
  is in approval 4): it rests on TC-01, a uniform patch of one material, which the bond carries
  to about 1e-7 (D3). The PLA-to-steel and orthotropic cases (D3) test it further.

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
| C. Tied penalty contact pair (chosen) | Passes the patch test (1e-7 stress, steel on steel). Either side can be the slave. Same input as later contact | Iterative step, about twice as slow as `*TIE`. A penalty stiffness to choose and approve |
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
iterative step about twice as slow as `*TIE`, and a new kind of flagged spot at a bond's edge. A
conforming mesh (D) would remove the penalty and the flags, but would tie every part's meshing
to every other's and could not handle STL.

## Test-first plan

Unit tests (fail before the code):
- `model`: `Transform` refuses a non-rotation and a mirror. A bonded connection with friction or
  preload is refused, and so is one with the same part on both sides. `Bond` and
  `SolveJob.bonds` default to empty. The spot `part` field.
- `physics`: merging two meshes gives offset nodes, prefixed sets and one section per part.
  Placement moves nodes and turns the material frame. The gap check refuses faces 0.1 mm apart
  and accepts touching ones. A small face on a larger one (rib on plate) is bonded with the
  small face as slave, and a staggered overlap is refused. With the same footprint, the finer
  side is the slave. K uses max(E1, E2, E3) for an orthotropic part. A part with no path to a
  support is refused.
- `solvers`: the deck writes the tied contact pair with element-face surfaces (face ids checked
  against the table in D3), and a job without bonds gives the same deck as today.
- `post`: the bonded region is not outer surface (sharp edges, wall span), but the plate around
  a rib is. `BOND_EDGE` flags the bond boundary and not its middle. Face peaks are keyed by
  part.
- `core`: point stresses are asked for when any part's material has a strength. Placement
  invariance. TS-01 turned 90° about z and moved gives the same tip deflection,
  turned. An orthotropic cube turned so that axis 1 lies along global x gives the TO-01
  `A1-side` result. Contact, bolt and pin connections are refused with the message.
- `cli`: `part place`, `connect bond`, `info` lines, and two-part result wording.

Benchmarks (`fea-verifier` writes them first; references proposed in the PR and created by the
maintainer):
- **TC-01**, bonded patch test with non-matching meshes. Proposed tolerances, each relative (stress
  to sigma, displacement to the largest exact displacement):
  - z tension, steel on steel: 1e-4 stress and displacement (measured 1.0e-7 and 1.6e-6).
  - x tension, steel on steel: 1e-5 stress and displacement (measured about 1e-6).
  - z tension, PLA on steel, and an orthotropic block on steel: 1e-4 (new cases, measured in
    PR-3 before the reference is proposed).
  - The catalog's old 1e-6 case for matching meshes is dropped: Tessera never merges nodes and
    meshes parts separately, so a matching interface goes through the same bond as any other.
- **TC-02**, the TS-01 cantilever cut at mid-span and bonded (measured: -0.40 %, the same as
  the single part): within 2 % of 1.526724 mm and within 0.5 % of the single-part result at
  the same mesh size. Unchanged from the catalog.
- **TC-01b**, partial cover: a 10 mm block bonded to the middle of a 30 mm plate face. Checked
  for: the bond is accepted with the block as slave, the force balance holds, and the corner
  is flagged `BOND_EDGE`. No closed-form stress is claimed for it.

## Pull requests

1. **PR-0 (docs):** this ADR, the TC-01 and TC-02 catalog entries updated with the D3
   measurements, and the ADR-009 D2 note on frames after placement.
2. **PR-1 (model):** `Transform` and `Connection` validators, `Bond`, `SolveJob.bonds`, spot
   `part`, `FlagKind.BOND_EDGE`, the `MaterialFrame` docstring.
3. **PR-2 (physics):** merge, placement, frame rotation, bonded → `Bond` with K, the covered
   side and gap check, slave choice, held-part check.
4. **PR-3 (solvers):** element-face surfaces and the tied contact pair in the deck
   (`calculix-specialist`). It also measures: behaviour at small gaps, conditioning, a
   supported node on a bonded face, and the PLA-on-steel and orthotropic TC-01 cases.
5. **PR-4 (post):** keys by part, bonded regions left out of the outer surface and wall span,
   bond-edge flags.
6. **PR-5 (core + edit):** multi-part solve, point stresses when any material has a strength,
   `place_part`, `add_bond`, the refusal for other connection kinds.
7. **PR-6 (CLI and docs):** commands and wording, and the user guide section on assemblies.
8. **PR-7 (benchmarks):** TC-01, TC-01b and TC-02 with proposed reference files.

PR-2 to PR-4 can run in parallel after PR-1. PR-3 and PR-7's decks can be tried early against
hand-built jobs.

## Approvals needed

1. Tied penalty contact pairs for bonded connections, not `*TIE`; the covered-side and slave
   rules; the refusals, including staggered overlaps (D3).
2. The penalty stiffness rule K = 1e5 × E_max / h, with E_max = max(E1, E2, E3) for an
   orthotropic material (D3).
3. The bond tolerance (5 % of the smaller median edge, at least 1e-3 mm) and refusing the solve
   when the faces do not touch (D3).
4. Trust at a bond (D5): the `BOND_EDGE` flag, leaving the bonded region out of the outer
   surface and wall span, and flagging nothing inside a bond.
5. Placement rules: mirrors refused, the print setup in part axes (D1).
6. Results and CLI wording for several parts (D4, D6).
7. TC-01 tolerances: z 1e-4, x 1e-5, the new PLA-on-steel and orthotropic cases at 1e-4, and
   dropping the matching-mesh case. TC-02 unchanged. TC-01b checks without a closed form.
   References come in PR-7.

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
3. [ ] ADR-004 item 3: record bonded connections as done after PR-5.
4. [ ] Tick ADR-009's action items, which shipped with TO-01 (#43).
