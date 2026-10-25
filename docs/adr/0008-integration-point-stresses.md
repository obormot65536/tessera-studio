# ADR-008: Stresses at integration points

**Status:** Accepted (maintainer, 2026-10-22). Amended the same day: integration-point output is
opt-in per solve job (decision 2).
**Deciders:** Maintainer

## Context

`SolveResult.stresses` holds one stress tensor per node, read from CalculiX's `.frd` file. To
produce it, ccx extrapolates each element's integration-point stresses to its nodes and averages
over the elements that share a node. ADR-004 left "stress recovery point (nodal vs integration
point)" to be revisited when failure criteria land, and its action item 7 records the first
symptom: nodal `.frd` stresses carry about 1e-5 relative bias even for an exact constant field
(TS-03). Three things make this more than a rounding issue:

- **Interfaces.** Nodal averaging mixes the stresses of elements that meet at a node. In v0.5,
  materials are assigned per region and FDM layer orientations per element, so at a region
  boundary or an orientation change the averaged tensor belongs to neither side.
- **Pointwise failure criteria.** Tsai-Wu and other layer-wise criteria (v0.5) are evaluated on one
  material's tensor in that material's axes. They need each element's own stress, not an average
  across materials.
- **Precision.** `.frd` stores `E12.5` (6 significant digits); the `.dat` file stores `E13.6`
  (7 digits).

Integration points lie inside the element. With the 4-point rule each point's smallest barycentric
coordinate is a = 0.138, so the point sits 0.138 of an element height from the nearest face (about
0.11 h for a regular tetrahedron with edge h). A maximum over integration points therefore
**underestimates** a surface peak: for a bending root it is lower by roughly the stress gradient
times that distance (TS-01 root roughly 5 % low; an estimate, not a measurement).

Constraints: ADR-001 (solvers and post import only `tessera.model`; domain packages never import
each other), ADR-002 (mm, MPa), ADR-004 (neutral, additive extensions of `SolveResult`), ADR-007
(peaks, flags and the two-mesh check are nodal), and hard rule 6 (no number a user trusts changes
without the maintainer's approval).

### What CalculiX 2.21 writes (verified)

Decks: one C3D10 element and a 2 x 1 x 1 block (12 C3D10 elements), with a homogeneous strain
prescribed on the boundary nodes and the three interior nodes solved for (9 equations). The block's
strain gives a rotated uniaxial stress σ = 14 n⊗n with n = (1, 2, 3)/√14, so that (sxx, syy, szz,
sxy, sxz, syz) = (1, 4, 9, 2, 3, 6) MPa: all six components are distinct and the xz/yz order is
visible. E = 1000 MPa, ν = 0.3. Output requests:

```
*EL PRINT, ELSET=EALL
S
*EL PRINT, ELSET=EALL
COORD
```

`.dat`, 2 x 1 x 1 block, linear `*STATIC` (excerpt):

```

                        S T E P       1


                                INCREMENT     1


 stresses (elem, integ.pnt.,sxx,syy,szz,sxy,sxz,syz) for set EALL and time  0.1000000E+01

         1   1  1.000000E+00  4.000000E+00  9.000000E+00  2.000000E+00  3.000000E+00  6.000000E+00
         1   2  1.000000E+00  4.000000E+00  9.000000E+00  2.000000E+00  3.000000E+00  6.000000E+00
         1   3  1.000000E+00  4.000000E+00  9.000000E+00  2.000000E+00  3.000000E+00  6.000000E+00
...
 global coordinates (elem, integ.pnt.,x,y,z) for set EALL and time  0.1000000E+01
```

The `.frd` of the same run, for comparison (nodal, extrapolated and averaged):

```
 -4  STRESS      6    1
 -5  SXX         1    4    1    1
 -5  SYY         1    4    2    2
 -5  SZZ         1    4    3    3
 -5  SXY         1    4    1    2
 -5  SYZ         1    4    2    3
 -5  SZX         1    4    3    1
 -1         1 9.99990E-01 3.99996E+00 8.99991E+00 1.99998E+00 5.99994E+00 2.99997E+00
```

`.dat`, one straight-sided C3D10 with corners at the origin and the unit axes (COORD block):

```
 global coordinates (elem, integ.pnt.,x,y,z) for set EALL and time  0.1000000E+01

         1   1  1.381966E-01  1.381966E-01  1.381966E-01
         1   2  5.854102E-01  1.381966E-01  1.381966E-01
         1   3  1.381966E-01  5.854102E-01  1.381966E-01
         1   4  1.381966E-01  1.381966E-01  5.854102E-01
```

Findings:

1. **Format.** Each block has a one-line header (`stresses (elem, integ.pnt.,sxx,...) for set
   <SET> and time <t>` or `global coordinates (elem, integ.pnt.,x,y,z) for set <SET> and time
   <t>`), a blank line, then one row per element and integration point: element id `I10`,
   point number `I3`, then values `1X,E13.6` (stress rows are 98 characters, COORD rows 56).
   Every value is preceded by a space, also negative ones (`-1.000000E+00`), so splitting on
   whitespace is safe. Ids up to 120000 were checked. A stress of -2e-101 was printed as
   `0.000000E+00`; no exponent without an `E` was seen.
2. **Component order.** `.dat` is **sxx syy szz sxy sxz syz**. The rotated uniaxial case prints
   `1 4 9 2 3 6`, while `.frd` (SXX SYY SZZ SXY SYZ SZX) prints the same state as `1 4 9 2 6 3`.
   The neutral order (xx yy zz xy yz zx, ADR-004) is `.dat` columns `[0, 1, 2, 3, 5, 4]`.
3. **Accuracy.** All 48 integration-point rows of the block equal (1, 4, 9, 2, 3, 6) exactly;
   the `.frd` nodal values of the same run are about 1e-5 low (0.99999, 3.99996, ...). This is
   the bias in ADR-004 item 7, and it is absent at integration points.
4. **COORD.** Accepted under `*EL PRINT` for C3D10 (no warning). The points are the 4-point rule
   with natural coordinates a = (5 - √5)/20 = 0.1381966 and **b = 1 - 3a = (5 + 3√5)/20 =
   0.5854102**. Point k is the one nearest corner node k (barycentric coordinate b at corner k,
   a at the other three):
   (ξ, η, ζ) = (a, a, a), (b, a, a), (a, b, a), (a, a, b), with ξ, η, ζ along the edges from
   node 1 to nodes 2, 3, 4. On a distorted tetrahedron with curved (perturbed) midside nodes,
   COORD equals the isoparametric tet10 map of these points to 3.5e-7 mm (the print rounding),
   while the straight-sided map is up to 0.033 mm off. COORD is the physical position.
5. **Repetition.** A linear `*STATIC` step writes one increment, even with an increment of 0.25
   requested. Each further step adds another `S T E P` section with both blocks. A nonlinear
   step (`NLGEOM`, `DIRECT`, 0.25) writes both blocks for every increment (4 times). A reader
   that keeps the last block per name gets the final state.
6. **Size.** S + COORD is 4 x (99 + 57) = **624 bytes per C3D10 element**: 74,880,268 bytes for
   120,000 elements, so about 0.6 GB at 1e6 elements. At 120,000 elements the run took about
   2 s longer (38.2 s against 36.4 s) on a mesh whose solve is trivial. This is not yet a measurement
   on a real part.
7. **Silent failures.** Without `*EL PRINT`, ccx still creates `job.dat`, with 0 bytes. An
   `*EL PRINT` without `ELSET=` gives only `*WARNING reading *EL PRINT: no set was defined`,
   exit code 0, and an empty `.dat`. The adapter checks only `*ERROR` and the exit code, so it must
   check the `.dat` contents.
8. **Set names.** `*ELEMENT, TYPE=C3D10, ELSET=EALL` works alongside user sets, which
   `_safe_names` always prefixes with `E_` (a user set named `EALL` becomes `E_EALL`; tested). If a
   second `*ELSET, ELSET=EALL` (any case) were ever written, ccx would silently **append** to it,
   and the printed blocks would contain duplicate rows (element 1 printed twice, no warning).
   The prefix rules this out, and the parser must reject duplicates as a second guard.
9. `*EL PRINT, ELSET=EALL` with `S, COORD` on one data line gives a byte-identical `.dat` to two
   separate `*EL PRINT` cards.

## Decision

Read integration-point stresses and their coordinates from the `.dat` file and carry them as a
new, optional, solver-neutral result. Peaks, flags and the two-mesh check stay nodal.

1. **Model (`tessera.model.solver`).** New frozen type
   `PointStresses(points: (M, Q, 3) mm, stresses: (M, Q, 6) MPa, neutral order xx yy zz xy yz zx)`,
   validated for shape. Here M is the number of elements and Q the number of points per element
   (4 for C3D10). `SolveResult` gains `point_stresses: PointStresses | None = None` as an
   additive last field. Metadata gets `point_rule` (e.g. `tet4`). `None` means "not available".
   Code that needs point stresses (failure criteria) must refuse rather than fall back to nodal
   values.
2. **Solvers (`tessera.solvers`, the only place that knows `.dat`).** The output is **opt-in**:
   `SolveJob.point_stresses: bool = False`. Only a job that asks for it gets `*EL PRINT` and has
   its `.dat` read; every other solve is unchanged and cannot fail on `.dat`. (Maintainer
   decision after review: reading it on every solve cost about 1.7 s and 300 MB at 100k elements,
   about 18 s and 3 GB at 1e6, for a result nothing used yet. Failure criteria will ask for it.)
   - Deck: `*ELEMENT, TYPE=C3D10, ELSET=EALL`, then `*EL PRINT, ELSET=EALL` with `S` and a
     second `*EL PRINT, ELSET=EALL` with `COORD` (finding 9 allows either form).
   - `parse_dat(text) -> DatResults` reads the `stresses` and `global coordinates` blocks,
     keeps the last block per name (finding 5) and parses with vectorized code. It reorders
     stresses with `[0, 1, 2, 3, 5, 4]` (finding 2). It requires exactly Q rows per element,
     numbered 1..Q, with no duplicates (finding 8). Missing entries become NaN.
   - When the job asks for point stresses, `run()` reads `job.dat`. A missing, empty or
     incomplete file is then a `SolverError` (finding 7), and `check_results_complete` is
     extended to point results.
3. **Post (`tessera.post`, points only, no ccx knowledge).** New `post/points.py`:
   - `element_nodal_stresses(mesh, ps) -> (M, 10, 6)`: the affine field through each element's Q
     points, in physical coordinates, evaluated at its 10 nodes. It is not averaged and is exact
     for linear stress fields. This is the input for v0.5 failure criteria.
   - `point_von_mises` for diagnostics.
   - `von_mises` in `summary.py` indexes `stresses[..., k]` instead of unpacking `.T`. Today it
     returns a transposed (Q, M) array for (M, Q, 6) input.
   - `flags.py` is split into `_flag_sources` and `flag_points(job, points)`, with
     `flag_nodes(job)` built on `flag_points(job, mesh.nodes)`. Its output stays identical.
4. **What stays nodal.** All reported peaks (`max_von_mises`, `peak_von_mises`,
   `clear_peak_von_mises`, face peaks), the node flags and the two-mesh check (ADR-007) keep using
   the `.frd` nodal stresses. **No reported number changes.** Replacing a peak with a maximum over
   integration points is forbidden (it is biased low at surfaces, see Context). Any such change
   needs a new ADR and the maintainer's approval (hard rule 6).
5. **Layering.** model (types), solvers (`.dat` parsing), post (point maths); core and the
   clients are unchanged. Storing point stresses in the results cache is part of ADR-003 action
   item 3.

## Options considered

### Option A: `.dat` via `*EL PRINT` S + COORD (chosen)
| Dimension | Assessment |
|---|---|
| Complexity | Medium: a second parser, one new model type |
| Accuracy | Exact at the points for linear fields; 7 significant digits |
| Interfaces | Each element keeps its own tensor |
| Cost | 624 bytes per element; about 0.6 GB at 1e6 elements |

### Option B: integration-point output in `.frd`
No integration-point output to `.frd` was found. `*EL FILE` gave only the nodal (extrapolated,
averaged) `STRESS` block in these runs; other `*EL FILE` options were not tried. Rejected.

### Option C: recompute stresses from displacements in Tessera
Duplicates the constitutive law (orthotropic, oriented in v0.5) outside the solver, and starts
from 6-digit `.frd` displacements. Rejected for production; kept as a test oracle.

### Option D: external readers (ccx2paraview, pygccx, meshio)
They read `.frd`, so they inherit the nodal averaging; they add a dependency for no gain.
Rejected.

## Trade-off analysis

Option A is the only one that returns CalculiX's own tensors, unaveraged and at full printed
precision, without re-implementing the material law. Its costs are a second parser that depends
on a text format and a large `.dat` file on big meshes. The format risk is limited by pinning
the observed layout and the xz/yz order in tests. The size is measured before release and can be
made opt-in. Keeping every reported peak nodal costs nothing now and avoids the low bias at
surfaces.

## Consequences

- Easier: v0.5 failure criteria get per-element tensors in each element's material; patch tests
  check stresses to 1e-6 instead of 1e-5; interface stresses are no longer averaged across materials.
- Harder: a second results file to parse and validate; `.dat` format changes across ccx versions
  are caught only by tests; disk use grows by about 624 bytes per element (396 as built, S only;
  measured in action item 6).
- Risks:
  - xz/yz swap. Pinned by the rotated uniaxial test.
  - `.dat` drift across ccx versions. Pinned by the parser snippet test and the ccx tests.
  - COORD missing in some ccx version. Fallback: map the rule points through the mesh (finding 4
    gives the rule and order).
  - `.dat` size at 1e6 elements. Vectorize, measure; if needed, opt in with a `SolveJob` flag.
    Measured: 402 MB and 8 s to read at 1e6 elements (action item 6).
  - A later "simplification" of peaks to integration-point maxima. Forbidden by Decision 4.
- Revisit: under `NLGEOM` (v0.5 contact), check whether COORD gives the deformed or the reference
  position; this has not been tested.

## Action items

1. [x] Maintainer accepts this ADR (including the verified `.dat` findings above).
2. [x] Tests first (solver tests marked `@pytest.mark.solver`):
   - model shape validation;
   - deck contents;
   - a `parse_dat` snippet test with distinct columns, two increments and a duplicate-row
     rejection;
   - ccx patch test: every integration point equals σ to 1e-6, element-nodal values are exact,
     and the `.frd` bias is absent;
   - rotated uniaxial tension with all six components distinct;
   - pure bending: integration points and element-nodal values exact for the linear field; the
     maximum over integration points lies below the surface peak;
   - COORD against the rule on a distorted, curved tet10;
   - post: affine-reproduction oracle, no averaging, `von_mises` shape, `flag_points` regression.
3. [x] `model`: `PointStresses`, `SolveResult.point_stresses`, `point_rule` metadata.
4. [x] `solvers`: deck output, `parse_dat`, `.dat` checks in `run()`.
5. [x] `post`: `points.py`, `von_mises` fix, `flag_points`. When items 3-5 have landed, tick
       ADR-004 action item 7.
6. [x] Benchmark gate: no existing reference changes. Propose
   `benchmarks/references/patch_test_points.toml` (TS-03 at integration points, 1e-6) in the PR;
   the maintainer creates it. Measure `.dat` size and run time on a 1e6-element part.
   Done: the reference was approved at 1e-6 (points) and 1.5e-6 (element-nodal, the E13.6
   floor times the tet4 fit weight); existing references unchanged. Measured (2026-10-25,
   CalculiX 2.21, 4 cores, 15 GB) on a clamped steel bar with an axial end force, a
   structured C3D10 mesh (6 tetrahedra per cube) through `CalculixAdapter`'s own deck and
   readers, each size run with `point_stresses` off and on:

   | elements | solver | ccx off / on | `.dat` | `.dat` read | `.frd` | `.frd` read |
   |---|---|---|---|---|---|---|
   | 10,800 | SPOOLES | 5.5 / 5.0 s | 4.3 MB | 0.06 s | 5.1 MB | 0.12 s |
   | 103,680 | SPOOLES | 249 / 251 s | 41 MB | 0.57 s | 46 MB | 1.2 s |
   | 103,680 | iterative Cholesky | 57.5 / 57.9 s | 41 MB | 0.65 s | 46 MB | 1.3 s |
   | 1,014,000 | iterative Cholesky | 970 / 991 s | 402 MB | 8.2 s | 433 MB | 19 s |

   The adapter prints S only (the points come from the mesh, not COORD), so the `.dat` costs
   **396 bytes per element** (4 rows of 99), not 624; it is slightly smaller than the `.frd`.
   ccx time and memory do not change measurably (+2 % at 1e6, within run-to-run noise; ccx
   peak 4.7 GB / 5.4 GB unchanged). Reading the `.dat` takes about half the time of reading the
   `.frd` and does not raise Python's peak memory (2.9 GB at 1e6, set by the `.frd` read); the
   arrays kept are 288 bytes per element (292 MB at 1e6). The 1e6 run used
   `SOLVER=ITERATIVE CHOLESKY` because SPOOLES needs 4.7 GB already at 1e5 elements and would
   not fit in 15 GB; the solver does not change what the `.dat` holds. Opt-in (decision 2)
   stays as it is; nothing here argues for more.
7. [x] Failure-criteria ADR (v0.5) built on `element_nodal_stresses`. Accepted as ADR-010
   (`0010-failure-criteria.md`).
