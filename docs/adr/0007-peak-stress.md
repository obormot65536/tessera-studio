# ADR-007: Reporting peak stress near supports and sharp edges

**Status:** Accepted (maintainer; defaults below approved)
**Deciders:** Maintainer

## Context

`tessera solve` reports `max von Mises stress`: the largest nodal von Mises stress anywhere in the
part (`tessera.post.summarize`, `StudyResult.max_von_mises`). In linear elasticity the stress has no
finite value at the edges and corners of a fully fixed face, at sharp inside (re-entrant) edges and
corners, and at line and point supports or loads (ADR-006). There the nodal peak keeps rising as the
mesh is refined. The user guide shows it on the example beam: from 2.5 mm to 1.25 mm the largest
displacement changes by 0.04 % while the peak von Mises rises by 19 %, at the clamped root.

Who reads this number:

- **Makers** compare it with a datasheet strength to decide whether a print will hold. A spurious
  peak makes them over-build; a mesh-dependent number they don't know about may be trusted.
- **MCP agents** pass it on as a fact unless the result says, machine-readably, that it is not one.
- **Future safety factors and failure criteria** (v0.5 orthotropic FDM, layer-wise criteria) use
  stress tensors, not von Mises, so whatever marks a spot untrustworthy must work for any field.

Today the only guard is prose: the user guide's two-mesh check (Step 7), which is manual and does
not say where the peak is or why it moved.

Constraints: ADR-001 (post imports only `tessera.model`; only core runs several solves), ADR-002
(mm, MPa), ADR-005 (places are named by face index), hard rule 6 (the thresholds decide which number
a user trusts, so they need the maintainer's approval; the defaults below are approved).

## Decision

The stress result reports *where* the peak is and whether that spot is known to be mesh-dependent.
The raw peak is kept but is no longer the only number shown. A two-mesh check is available on
request.

1. **Flagging nodes (post, every solve).** `tessera.post` builds a per-node mask from the mesh and
   the solve job:
   - **Support:** a node within distance d of a node of any support that holds two or three
     directions. Supports holding one direction (symmetry, rollers) are not flagged: they are exact
     symmetry planes. Every line support (ADR-006) is flagged, whatever its directions.
   - **Sharp inside edge:** a node within d of an edge of the outer surface where the surface turns
     *inward* (concave) by at least θ_edge, measured on the mesh's outer triangles. Smooth fillets,
     holes and outward (convex) edges are not flagged.
   - **Distance:** d = k_d × h, where h is the mesh's median corner-to-corner edge length (works
     with an automatic mesh size too).
   - The mask is independent of the scalar reported; safety-factor and failure-criteria code use
     the same mask.
2. **Reporting (post, every solve).** `StudyResult` gains, additively:
   - `peak_von_mises: StressSpot`: the raw peak with value, node, location (mm), face (or `None`
     inside the part) and its flags (support id, or "sharp inside edge").
   - `clear_peak_von_mises: StressSpot | None`: the peak over unflagged nodes; `None` if every node
     is flagged.
   - `face_peaks: tuple[StressSpot, ...]`: the peak on each geometry face, sorted by value.
   - `node_flags`: the per-node mask.
   - `max_von_mises` keeps its meaning (the raw peak), so no existing number changes.
3. **Two-mesh check, on request (core runs, post compares).** `tessera.core.solve(project,
   check_mesh=True)` / `tessera solve --check-mesh` solves at h and again at h/√2. Each reported
   spot gets a verdict: `settled` if it changes by at most X between the meshes, otherwise
   `mesh-dependent`. The verdict overrides the flags. All reported numbers come from the fine mesh.
   *Amendment (maintainer-approved):* the fine mesh scales every element size the mesher uses
   (`mesh_step(..., size_factor=1/√2)`, which also raises gmsh's minimum points per circle and
   curve), not only the size cap: a cap alone leaves the elements at small holes unchanged (0.96
   times), so a peak there would "settle" untested. A spot is judged only where its local element
   edge length shrank to at most R times the coarse one; elsewhere it is `not-refined` and its
   flags still apply.
4. **What clients show.** The CLI leads with the clear peak, labelled "away from supports and sharp
   inside edges", then the raw peak with location, face and a one-line explanation if flagged, then
   the top face peaks, then the mesh-check verdicts if run. MCP carries the same data as structured
   fields (`flags`, `trust: "clear" | "flagged" | "settled" | "mesh-dependent"`; a
   `not-refined` spot keeps "clear" or "flagged"), and its tool
   descriptions tell agents to quote only `clear` or `settled` values against a strength.
5. **Downstream rule.** A future safety factor or failure index is reported from unflagged or
   `settled` spots; flagged spots that are not settled are listed separately and never become the
   headline.
6. **Defaults (approved by the maintainer):**

| Parameter | Value |
|---|---|
| Supports flagged | Holding 2 or 3 directions; every line support |
| Inside-edge threshold θ_edge | 30° inward turn |
| Flag distance | d = 2 × h (median corner edge length) |
| Mesh-check refinement | h → h/√2 |
| Mesh-check threshold X | 5 % change |
| Mesh-check local refinement needed for a verdict R | 0.85 (fine / coarse edge length) |
| Face peaks shown by the CLI | Top 5 (all in the API) |
| Mesh check default | Off; opt in with `--check-mesh` |
| Thin-wall minimum N (Amendment 2) | 5 elements across the wall at a flagged spot, on the coarse mesh (approved by the maintainer) |

## Options considered

- **A: keep the peak, flag spots near supports and sharp inside edges, report the peak away from
  them.** Free per solve and says *why* a peak is suspect; heuristic, so it misses singularities it
  has no rule for.
- **B: report an averaged stress instead of the nodal peak.** Depends on an arbitrary radius, cuts
  real stress raisers (about 19 % at 0.5 mm from a 10 mm hole edge), and is only principled with a
  material length the library does not have. Rejected as the headline.
- **C: peak per face with location.** Cheap; tells *where* in face numbers users know; says nothing
  about trust on its own.
- **D: automatic two-mesh check.** Measures convergence instead of guessing, but costs 3-5× per
  solve if always on.
- **E: A + C every solve, D on request (chosen).** Location and reason always shown for free;
  measured trust one flag away, overriding the guess.

## Trade-off analysis

Flagging costs nothing and catches the common cases makers hit (clamped mounting faces, sharp inside
corners), but its clear peak still drifts slowly with h because d scales with h. The two-mesh check
is the honest test but multiplies solve time; making it opt-in with its verdict overriding the flags
gets most of its honesty without the cost on every iteration. The raw peak stays, with its location
and a warning, so a real stress raiser inside a flagged band (a hole next to a clamp) is not hidden.

## Consequences

- Easier: makers see where the peak is and why it may be wrong; agents get a trust label; future
  safety factors reuse one mask and one verdict; the manual two-mesh advice becomes one option.
- Harder: concave-edge detection needs outward triangle normals and a distance query; every new
  kind of singular input (line supports, point loads, contact) needs a flag rule.
- The clear peak is not mesh-independent; docs say "away from supports and sharp inside edges", not
  "true stress".
- Not covered yet: edges of loaded faces (weak singularities; needs `load:<id>` node sets).
- Revisit the thresholds once TS-08 and a re-entrant-corner benchmark give measured growth rates.

## Action items

1. [x] Maintainer approves the defaults.
2. [x] `model.results`: `StressSpot` (and flag/verdict types); extend `StudyResult` additively.
3. [x] `post`: node flags, stress spots, face peaks; `summarize` takes the solve job.
4. [x] CLI output: clear peak first, raw peak with location and warning, top 5 faces.
5. [x] `core.solve(..., check_mesh=False)` and `tessera solve --check-mesh` (second PR).
6. [x] Benchmarks: TS-01 root peak flagged; TS-02 hole peak not flagged.
7. [x] Docs: user guide Steps 7 and 9, API reference.
8. [ ] MCP: `flags` and `trust` fields, and tool descriptions that quote only `clear` or `settled`
       values (the MCP server does not exist yet).
9. [x] Amendment 2: maintainer approves N = 5 and the flagged-only scope.
10. [ ] Amendment 2: `model.results` adds `WallSpan`, `StressSpot.wall`, `Verdict.THIN_WALL`,
        `SpotCheck.elements_across` and `MeshCheck.wall_elements_min` (all additive).
11. [ ] Amendment 2: `post` measures the wall at reported spots (`wall.py`, largest inscribed
        ball) and `compare_meshes` applies the guard; CLI wording; benchmark assertions on TS-01,
        TS-08A and TS-02.
12. [ ] Follow-up: local refinement near flagged spots, so N = 5 does not force a fine global mesh.

## Amendment 2: thin walls

**Status:** Accepted (maintainer, 2026-10-22), including N = 5. Revised the same day after review of
the first implementation: the maintainer chose the largest-ball wall measure over the ray cast, no
doubling across one-direction supports, and no verdict where no wall can be measured; then that the
faces of the flagging support are left out of the wall measure.

**Problem.** A mesh too coarse to show the stress rise at a singularity can read `settled`: the
peak there barely moves between the two meshes because neither resolves it. The user guide's
example beam does this at 5 mm, two elements through its 10 mm depth (`settled, changed 1.12 %`
at the clamp); at 2.5 mm the same peak is `mesh-dependent`. Today only the user guide warns.

**Decision.**

1. **Measure the wall at each reported spot** as the diameter of the largest ball inside the part
   that touches the spot's surface node (the "shrinking ball": centre on the inward normal, grown
   until another surface sample would fall inside it; samples on the spot's own faces or within
   about one edge length are ignored, so an inside edge does not collapse the ball). This reads the
   beam depth at a node in the middle of a clamped end face, where a straight inward ray would read
   the bar's length. At a flagged spot, the faces of the support that flags it are left out (the
   part continues into the fixture there), so a peak beside a clamp reads the section depth, not
   twice its distance from the clamp. A flagged spot inside the part is measured at its nearest
   surface node. Nothing is doubled across a one-direction support: such a support may be a symmetry
   plane or a real roller, and doubling across a roller would hide a thin wall, so the guard errs
   toward judging less. Elements across = thickness / local median corner edge (`local_edge_length`,
   as for the refinement ratio). Post imports only `tessera.model`, so ADR-001 holds and STL parts
   work. One KD-tree over the surface per mesh; a few milliseconds per spot (3 ms at 20k nodes, 30
   ms at 400k), measured only at the raw peak, the clear peak and flagged spots.
2. **Guard flagged spots only.** If a spot is flagged and its change is within the threshold, but
   the coarse mesh had fewer than N elements across the wall there, the verdict is the new
   no-verdict value `thin-wall`. Like `not-refined`, the spot's flags still apply and MCP's `trust`
   falls back to `flagged`. A `mesh-dependent` reading stays `mesh-dependent`; a measured failure
   needs no guard. A flagged spot whose wall cannot be measured (for example an open STL) is treated
   as thin: it gets no verdict rather than `settled`. Unflagged spots are never guarded: the sweep
   below shows truly settled unflagged spots with 0.5 to 3.3 elements across (TS-02 hole, TS-05
   midspan, TS-01 and TS-08B clear peaks), because a smooth bending or membrane stress is
   represented well by quadratic tets even with one or two layers.
3. **Warn without the check too.** A flagged spot with fewer than N elements across gets one line
   in `tessera solve` saying how many there are and the mesh size that would give N, so the user
   refines before running `--check-mesh`.
4. **N = 5** (approved by the maintainer). Every number stored and shown is the measured count,
   not a pass/fail, so a later N changes no stored data.

**Evidence.** A sweep of `check_mesh=True` solves through the benchmarks' own helpers, recording
each flagged spot's coarse elements across, change and verdict, and the same spot's verdict at
each finer size:

| Case | Coarse h (mm) | Across (flagged peak) | Raw peak | Flagged spots |
|---|---|---|---|---|
| TS-01 cantilever, 10 mm deep | 5 / 3.54 / 2.5 / 1.77 | 2.0 / 2.9 / 3.9 / 6.0 | settled +1.4 % / settled −2.9 % / mesh-dep. +8.5 % / mesh-dep. +9.2 % | all settled at 2.0; mixed at 2.9; all mesh-dep. from 3.9 |
| TS-05 beam, line supports | 3.33 / 2.5 / 1.77 | 3.2–4.0 / 4.1–4.5 / 5.9 | mesh-dep. +27 % / +16 % / +25 % | face 5 settled +2.0 % at 3.21; face 3 settled +3.2 % at 4.50; all mesh-dep. at 5.9 |
| TS-08A disc, clamped rim, t = 2 | 2 / 1 / 0.667 / 0.471 / 0.333 | 1.0 / 1.9 / 3.0 / 4.1 / 5.8 | settled +3.7 % / settled +3.9 % / mesh-dep. +5.6 % / +10.9 % / +11.3 % | face 4 settled +4.5, +4.1, +4.4 % up to 4.56 across; mesh-dep. +10.1 % at 5.8 |
| TS-08B disc, ring line support | 2 / 1 / 0.667 | 0.9 / 1.8 / 3.0 | mesh-dep. at every size (+8.7 to +25 %) | no settled flagged reading |
| TS-02 strip with hole (control) | 2 | 0.25 modelled (0.5 full) | settled −1.3 %, unflagged | no flags; correctly settled |

Every flagged spot at 5 or more elements across reads `mesh-dependent` (18 of 18). Of the 28
flagged readings below 5 that read `settled`, none is truly settled: each reads `mesh-dependent`
at a finer size, 15 of them at the next one. N = 3 would pass false settles at 2.91 (TS-01 at
3.54 mm) and 3.21 (TS-05 face 5); N = 4 at 4.50 (TS-05 face 3) and 4.56 (TS-08A face 4, which
read +4.40 % and then +10.07 %). The benchmark helper loads TS-01 slightly differently from the
user guide's example beam, so it reads +1.44 % at 5 mm where the guide's beam reads 1.12 %; both
are re-run and reproduce, with the same verdict.

**Limits of the evidence.** Only three geometries (TS-01, TS-05, TS-08A) were measured at 5 or more
elements across, and there is no measurement between 4.56 and 5.8 across: N = 5 clears the
weakest case seen (TS-08A face 4, hovering at 4.1–4.5 % before rising) by one √2 step. A
weaker singularity could still hover under 5 % beyond 5 across. Revisit N when a re-entrant
corner benchmark exists.

**Cost to makers.** Printed walls are often 1.2 to 2 mm, so N = 5 means about 0.25 to 0.4 mm
elements near a flagged spot; with today's global mesh size that is minutes per solve for a
small part (TS-08A at 0.333 mm: about 19 min for both solves). The guard will therefore fire
often on thin printed parts. Local refinement near flagged spots (action item 12) is the
follow-up that makes resolving them affordable; until then the guard tells the truth, which is
that the check cannot judge that spot at that mesh.
