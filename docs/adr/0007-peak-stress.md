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
4. **What clients show.** The CLI leads with the clear peak, labelled "away from supports and sharp
   inside edges", then the raw peak with location, face and a one-line explanation if flagged, then
   the top face peaks, then the mesh-check verdicts if run. MCP carries the same data as structured
   fields (`flags`, `trust: "clear" | "flagged" | "settled" | "mesh-dependent"`), and its tool
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
| Face peaks shown by the CLI | Top 5 (all in the API) |
| Mesh check default | Off; opt in with `--check-mesh` |

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
2. [ ] `model.results`: `StressSpot` (and flag/verdict types); extend `StudyResult` additively.
3. [ ] `post`: node flags, stress spots, face peaks; `summarize` takes the solve job.
4. [ ] CLI output: clear peak first, raw peak with location and warning, top 5 faces.
5. [ ] `core.solve(..., check_mesh=False)` and `tessera solve --check-mesh` (second PR).
6. [ ] Benchmarks: TS-01 root peak flagged; TS-02 hole peak not flagged.
7. [ ] Docs: user guide Steps 7 and 9, API reference.
