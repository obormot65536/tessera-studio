# ADR-010: Failure criteria and strength ratios (v0.5)

**Status:** Accepted (maintainer, 2026-10-25)
**Deciders:** Maintainer

## Context

Tessera's question is whether a part will hold, and today it answers with von Mises stress and
leaves the comparison with a strength to the user. For printed parts that comparison is wrong
(ADR-009 D4: a print is often less than half as strong across its layers as along them), and
the CLI says so. ADR-008 action item 7 asks for a failure-criteria ADR built on
`element_nodal_stresses`; this is it.

What exists today:

- **Stresses.** `SolveJob(point_stresses=True)` returns each element's own stresses at its four
  integration points (ADR-008), and `tessera.post.element_nodal_stresses` gives the unaveraged
  affine field at the element's ten nodes, exact for linear fields. Asking for them costs about
  nothing in the solve and 8 s of reading at 1e6 elements (ADR-008 item 6).
- **Frames.** `Section.orientation: MaterialFrame | None` and
  `tessera.model.elasticity.rotate_stress(stress, frame)`, which returns 11 22 33 12 23 31 in
  the material axes (1 = raster, 2 = across the road in the layer, 3 = build) (ADR-009).
- **Trust.** ADR-007 flags nodes near singular supports and sharp inside edges, gives a two-mesh
  verdict on request, and fixes the rule this ADR must follow: "a future safety factor or failure
  index is reported from unflagged or `settled` spots; flagged spots that are not settled are
  listed separately and never become the headline" (decision 5).
- **No strength data anywhere.** The library holds two metals with stiffness and density only;
  project definitions (`IsotropicDefinition`, `OrthotropicDefinition`, schema 6) hold no
  strengths. `materials-curator` is set up to record per-direction tensile, compressive and
  shear strengths, but no entry exists.
- **One part per study** (v0.1), so in practice one material and one frame per solve; the design
  below works per section so assemblies need no change.

Constraints: ADR-001 (post imports only `tessera.model`; core orchestrates), ADR-002 (MPa),
ADR-004 (no CalculiX names in `model`), ADR-007 decision 5, ADR-008 decision 4 (reported peaks
stay at nodes; no switch to integration-point maxima), hard rule 6 (the criterion, its
interaction coefficient and any default that fills in a missing strength decide numbers users
trust).

## Decision

### D1. Report a strength ratio, not a pass/fail or a safety factor

For each point Tessera reports the **strength ratio** R: the factor by which every load of the
study (gravity included) could be multiplied before the criterion predicts first failure there.
The analysis is linear, so R is exact for the criterion: R < 1 means failure is predicted at the
applied loads, R = 2 means twice the loads would reach it. The headline is the lowest R in the
part.

- **Why not a "safety factor".** A safety factor is a design target that also covers load,
  material and print-quality uncertainty; that choice is the user's. Tessera reports the ratio
  and does not judge it. No pass/fail threshold is introduced, so none needs approval.
- **Why not a failure index (1/R).** For Tsai-Wu the index is not proportional to load (the
  linear terms), so "index 0.5" does not mean "half way"; R does.

### D2. Criteria by material kind

| Material | Criterion | R |
|---|---|---|
| Orthotropic (FDM) | **Tsai-Wu**, interaction F*_ij = -1/2 | Positive root of a R^2 + b R - 1 = 0 |
| Isotropic, ductile (metals) | **von Mises** against the yield strength | yield / sigma_vm |
| Isotropic, brittle (resin) | **Maximum principal stress** (Rankine) | min(T / sigma_1, C / -sigma_3) over the positive terms |

Tsai-Wu, in the material axes (Voigt 11 22 33 12 23 31, compression strengths positive):

- F_i = 1/T_i - 1/C_i (i = 1, 2, 3), F_ii = 1/(T_i C_i), F_44 = 1/S12^2, F_55 = 1/S23^2,
  F_66 = 1/S13^2, F_ij = F*_ij sqrt(F_ii F_jj) for i != j in 1..3, no linear shear terms.
- a = sigma . F_quad . sigma, b = F_lin . sigma, R = (-b + sqrt(b^2 + 4a)) / (2a); R = inf
  where the stress is zero.

Why Tsai-Wu for FDM: it is the criterion most FDM strength studies fit (e.g. Ahn et al., Rapid
Prototyping Journal 8 (2002) 248-257, ABS, `unverified`), it uses exactly the nine strengths the
curator records, it handles combined in-layer tension and interlayer shear (the usual way a
printed bracket breaks) in one smooth envelope, and it is a single expression with no mode
switching. Why F*_ij = -1/2: with T = C = X and S = X/sqrt(3) in every direction it reduces
**exactly** to von Mises (checked numerically for random stress states while drafting this
ADR), so isotropic and orthotropic parts are judged consistently; it is the common default (Tsai
and Hahn, *Introduction to Composite Materials*, 1980, `unverified`). Real F12 values are rarely
measured for FDM.

**Mode label.** Tsai-Wu alone does not say *how* a part fails, and makers need "across the
layers" to fix a print. At the lowest-R spot Tessera also reports the term with the largest
maximum-stress ratio |sigma_k| / strength_k, named:

| Term | Label |
|---|---|
| sigma_11 | along the raster (tension / compression) |
| sigma_22 | across the roads in a layer (tension / compression) |
| sigma_33 | across the layers (tension / compression) |
| tau_12 | in-layer shear |
| tau_13, tau_23 | interlayer shear |

The label explains R; it never changes it.

### D3. Strength data: optional, per material, in model

New frozen types in `tessera/model/strength.py`, exported from `tessera.model`:

- `OrthotropicStrength(t1, t2, t3, s12, s13, s23, c1=None, c2=None, c3=None)`, MPa, all > 0.
- `IsotropicStrength(tensile, compressive=None, brittle=False)`: `tensile` is the yield strength
  for a ductile material and the tensile strength for a brittle one.

`IsotropicMaterial` and `OrthotropicMaterial` gain `strength: ... | None = field(default=None,
kw_only=True)` of the matching kind. The solver ignores it; decks stay byte-identical (tested).
Carrying it on the material, rather than in a parallel table, means post already finds it through
`Section -> material` together with the frame.

**Missing compressive strengths.** FDM datasheets and most coupon sets give tensile values only.
A missing C_i is taken equal to T_i, and the result lists every assumed value
(`StrengthResult.assumed`) and the CLI prints them. For FDM this is conservative in the common
case, since measured compressive strengths are at or above tensile ones, especially across layers
(to be confirmed by the curator's data). Shear strengths are required: there is no safe default
for interlayer shear.

**Project files (schema 7).** `IsotropicDefinition.strength` and `OrthotropicDefinition.strength`
(optional, defaults `None`), validated as above; identity migration pinned by
`test_schema_6_project_loads_unchanged`. Strengths a user types in are their own data (ADR-009
D3); library strengths are hard rule 6 data and arrive only through `materials-curator` in their
own PR. This ADR adds no library numbers.

### D4. Where R is evaluated

- **Input.** `element_nodal_stresses` (ADR-008): each element's own stress at its ten nodes,
  rotated into that element's material frame with `rotate_stress` (the section's orientation;
  global axes when `None`). No averaging across elements, so a node on an interface between two
  materials or frames is judged in each material separately.
- **Per node.** R at a mesh node is the lowest over the elements meeting there. Nodes are the same
  places ADR-007 flags and the same places the von Mises peaks are reported, so flags, faces,
  walls and the two-mesh check apply unchanged.
- **Not at integration points** for the headline: their minimum lies inside the elements and is
  biased high at surfaces, exactly where printed parts fail (ADR-008 decision 4).
- **Refuse rather than fall back.** If the solution has no point stresses, the strength check is
  not computed (`None`) and says why; nodal `.frd` stresses averaged across materials are never
  used for it.

### D5. Results and trust (ADR-007 decision 5)

New frozen types in `tessera.model.results`:

- `StrengthSpot(ratio, node, element, location, face, mode, flags=(), wall=None)`, with
  `mode: FailureMode` (D2 labels as a `StrEnum`, plus `VON_MISES`, `MAX_PRINCIPAL_TENSION`,
  `MAX_PRINCIPAL_COMPRESSION` for isotropic materials).
- `StrengthResult(criteria: dict[str, str] (material -> criterion), lowest: StrengthSpot,
  clear_lowest: StrengthSpot | None, face_lowest: tuple[StrengthSpot, ...], assumed:
  tuple[str, ...], unchecked: tuple[str, ...])`. `unchecked` names materials without strength
  data; their elements are left out and the CLI says so.
- `StudyResult.strength: StrengthResult | None = field(default=None, kw_only=True)`.
- `MeshCheck.lowest_strength_ratio: SpotCheck | None = field(default=None, kw_only=True)`, judged
  with the approved 5 % threshold, the refinement limit and the thin-wall guard, exactly as the
  von Mises peaks.

Headline: `clear_lowest` (unflagged), or `lowest` when the mesh check rates it `settled`.
`lowest` is always shown with its flags; a flagged, unsettled spot is never the headline.

### D6. Core and clients

- `core.study.solve`: when any material of the study has a strength, the job asks for point
  stresses (`point_stresses=True`) and `post.summarize` adds the `StrengthResult`. Studies without
  strengths are unchanged (no `.dat`, no extra time).
- `tessera solve` prints, before the von Mises lines (wording to approve):

  ```
  lowest strength ratio 2.41 (Tsai-Wu, across the layers, tension), away from supports and
    sharp inside edges, at (12.0, 4.5, 30.0) mm, face 7
    The loads could be about 2.4 times higher before the first failure is predicted there.
  lowest strength ratio 1.62 at (0.0, 2.0, 1.0) mm, face 3: near support 'base'; mesh-dependent
  assumed: compressive strengths c1, c2, c3 equal to the tensile ones (material 'petg-x')
  ```

  ADR-009's orthotropic note is kept only when no strength is given.
- MCP (when it exists) carries `strength` with the same `trust` field as the peaks.
- Results cache (ADR-003 item 3) stores the per-node R with the stresses.

## Options considered

### Criterion for FDM

| Option | For | Against |
|---|---|---|
| **Tsai-Wu, F* = -1/2 (chosen)** | Interaction of in-layer and interlayer stresses; used in FDM papers; reduces to von Mises; one expression | F12 rarely measured; R needs a quadratic; mode needs a separate label |
| Maximum stress | Mode is explicit; no interaction parameter; simple | Ignores interaction: tension across the roads plus interlayer shear is judged by each alone, unconservative for the bracket case |
| Tsai-Hill | No linear terms, one fewer choice | Tension and compression must be picked per sign, which makes it piecewise; no advantage over Tsai-Wu once compressive strengths exist |
| Hashin / Puck | Physically based modes for fibre composites | Built for fibre/matrix failure; FDM has no fibres; needs data FDM tests do not give |

### Where to evaluate

| Option | For | Against |
|---|---|---|
| **Element-nodal, lowest per node (chosen)** | Per-material, unaveraged; same places as flags and peaks; exact for linear fields | Extrapolated values at nodes can overshoot in steep gradients (as do the `.frd` peaks) |
| Integration points | Exact solver output | Biased high R at surfaces; ADR-008 decision 4 forbids it for headlines |
| `.frd` nodal stresses | No `.dat` | Averaged across elements and materials; wrong frame mixing at interfaces |

### Storage of strengths

On the material dataclass (chosen) versus a separate `strengths` map in `SolveJob` or passed to
`summarize`. The separate map keeps solver types pure but duplicates the section-to-material
lookup and can disagree with it.

## Trade-off analysis

The choice that matters most to makers is reporting a ratio with a mode, not a verdict: it tells
them how far from failure the part is and what to change (orientation, wall count) without
Tessera picking a safety factor for them. Tsai-Wu costs one interaction coefficient, which is
the main number needing approval; -1/2 is defensible because it makes the isotropic limit equal
von Mises, so a resin or metal part and an FDM part with equal strengths in every direction get
the same answer. Defaulting missing compressive strengths to tensile ones makes the feature usable
with typical coupon data, and the assumption is listed in every result rather than hidden.

## Test-first plan

Unit tests (fail before the code):

- `model/tests/test_strength.py`: validation (non-positive, NaN, wrong kind for the material);
  `strength` defaults to `None`; schema 7 round-trip and `test_schema_6_project_loads_unchanged`.
- `post/tests/test_failure.py`: uniaxial stress along each material axis gives R = T_i/sigma or
  C_i/|sigma| exactly; pure shear gives S/tau; Tsai-Wu with isotropic strengths equals von Mises
  for random states (to 1e-12); R(k sigma) = R(sigma) / k for every k > 0, linear terms
  included (corrected in PR-2: the drafted "only when b = 0" was wrong), and R sigma lies on the
  envelope; reversing the stress changes R when tension and compression strengths differ; zero
  stress gives inf; brittle isotropic picks tension or compression;
  mode labels for each term; missing C_i fills in and is listed.
- `post/tests/test_strength_result.py`: per-node minimum over elements (two elements, two frames
  meeting at a node); flagged lowest never becomes `clear_lowest`; `None` without point stresses;
  `unchecked` materials are skipped.
- `solvers/tests/test_calculix_unit.py`: decks with and without strength are identical.
- `core` and `cli`: point stresses are requested only when a strength exists; CLI lines and the
  assumed-values line.

Benchmark **TO-06, off-axis strength** (new catalog entry in this PR; `fea-verifier`; reference
proposed in the implementing PR, created by the maintainer): a cube under uniform uniaxial stress
sigma = 10 MPa, made-up TO strengths (not user data), 3-2-1 nodal supports and equal and opposite
face tractions so that the off-axis shear coupling is free and the stress is exactly uniform:

| Case | Frame | Closed-form R |
|---|---|---|
| raster 0, 30, 45, 60, 90 deg | build Z | Tsai-Wu of (sigma cos^2, sigma sin^2, 0, -sigma sin cos, 0, 0) |
| printed on its side | build (1, 0, 0) | T3 / sigma (across the layers) |

With T1 40, C1 50, T2 30, C2 45, T3 15, C3 40, S12 20, S13 = S23 12 MPa this gives R = 4.000000,
3.380006, 3.096069, 2.982983, 3.000000 and 1.5 (`computed`). The field is exact, so the
tolerance proposed is the element-nodal precision floor already approved for ADR-008 (1.5e-6
relative). One end-to-end case (printed on its side) runs through `tessera.core.solve` with
rollers, which are exact when the load is along a material axis.

## Pull requests

1. **PR-0 (docs):** this ADR and the TO-06 catalog entry.
2. **PR-1 (model):** `strength.py`, `strength` on the materials, schema 7 definitions,
   `StrengthSpot`, `StrengthResult`, `FailureMode`, `StudyResult.strength`,
   `MeshCheck.lowest_strength_ratio`.
3. **PR-2 (post):** `failure.py` (criteria, R, modes) and the `StrengthResult` in `summarize`;
   mesh-check comparison.
4. **PR-3 (core + materials):** request point stresses when a strength exists; definitions carry
   strengths into materials; `edit.add_material` accepts them.
5. **PR-4 (CLI):** the lines in D6; the orthotropic note only without strengths.
6. **PR-5 (TO-06):** benchmark and proposed reference toml.
7. **Separate issue:** `materials-curator` drafts PLA and PETG X/Y/Z strengths with sources (the
   v0.5 gate's "X/Y/Z coupons tested"), approved number by number.

## Approvals needed

1. Tsai-Wu for orthotropic materials with F*_ij = -1/2, and von Mises / maximum principal stress
   for ductile / brittle isotropic ones (D2).
2. Missing compressive strengths taken equal to tensile ones, listed in every result (D3).
3. Strength ratio as the reported quantity, its definition (all loads scaled, gravity included)
   and the CLI wording (D1, D6).
4. The two-mesh check of the lowest ratio reusing the approved 5 % threshold (D5).
5. The TO-06 made-up strengths and tolerance (in the implementing PR).

## Consequences

- Easier: makers get "how far from failure, where, and which way it breaks" for printed parts;
  agents get a number with a trust label; region-wise materials (TO-05) and assemblies need no
  new design, because R is per section.
- Harder: schema 7 must be read forever; every material entry gains up to nine more numbers to
  source and approve; a strength check makes the solve read the `.dat` (8 s at 1e6 elements).
- First-ply onset only: no progressive damage, no fatigue, no creep, no buckling, and no
  account of print defects beyond what the coupon strengths contain. The user guide must say so
  next to the first strength ratio it shows.
- Revisit F*_ij when FDM biaxial data exists, and the criterion when contact (v0.5) adds bearing
  failure at bolt holes.

## Action items

1. [x] Maintainer approves the items above (approvals 1-4; 5 comes with the TO-06 PR).
2. [ ] PR-1 to PR-5.
3. [ ] Curator issue for PLA/PETG strengths.
4. [ ] ADR-009 D4: replace the orthotropic note with the strength lines once PR-4 lands.
5. [ ] ADR-007: record that strength ratios follow decision 5 (cross-reference).
