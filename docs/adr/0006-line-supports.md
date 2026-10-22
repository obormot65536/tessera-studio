# ADR-006: Line supports (supports on geometry edges)

**Status:** Accepted (maintainer)
**Deciders:** Maintainer

## Context

Supports hold faces today: `Support.faces` is a `FaceSet` of geometry face indices, physics turns
the face's tri6 surface into a node set, and CalculiX fixes its DOFs. Several validation cases need
a support on a *line* instead:

- TS-05 (simply supported beam): one pinned line and one roller line at the bottom of the ends.
- TS-08 variant B (simply supported disc): u_z = 0 on the bottom outer edge (a ring).
- Later: TS-12 / NAFEMS LE10 (u_z = 0 only on the mid-plane line of the outer face) and TM-03.

Makers want the same thing: a shelf resting on a ledge, a part sitting on a rod. A face support is
a clamp; a line support lets the part rotate about the line.

Constraints from earlier decisions:

- ADR-001: `geometry` and `mesh` are siblings and cannot share code except through `model`.
- ADR-005: STL faces are Tessera's own (sharp-edge components of triangles). gmsh rebuilds STL
  solids in patches that are not Tessera faces, and its curves are not Tessera edges either.
- ADR-003: `project.json` carries a `schema_version`; older files must keep loading.
- ADR-004: turning user intent into node sets is physics' job; the solver sees node sets only.
- The mesh already exposes, for every face `i`, a tri6 surface and a node set `face:<i>`, for STEP
  and STL alike. Neighbouring faces share edge nodes, corner and midside, because the mesh is
  conforming.

## Decision

1. **An edge is named by the two faces it separates.** New model type
   `EdgeSet(part_id, edges: tuple[tuple[int, int], ...])`. Each pair is two different face numbers
   of the part, stored ascending; pairs are unique. The JSON is `"edges": [[3, 5], [5, 6]]`. A pair
   stands for everything the two faces share; if they share several separate curves, all of them
   are held.
2. **`Support` holds faces or edges.** `Support.faces: FaceSet | None` and
   `Support.edges: EdgeSet | None`; exactly one is set. A `Support.part_id` property replaces
   `support.faces.part_id` in callers. `dofs` keeps its meaning (xyz: pinned line; z: roller line).
3. **Schema 5.** `SCHEMA_VERSION = 5`. Files up to 4 load unchanged because they only have
   `faces`, so no migration function is needed. Older Tessera refuses a schema-5 file with its
   existing "reads up to 4" message.
4. **Listing edges.** `tessera.geometry.describe_edges(path) -> list[EdgeInfo]`, with `EdgeInfo` in
   `tessera.model.geometry_info`: `faces (a, b)`, `kind`, `length` (mm), `start`/`end` (`None` for
   a closed edge such as a circle), `point` (a point on the edge), `curves` (number of separate
   curves the pair shares).
   - STEP: gmsh curves whose upward adjacency is two different faces, grouped by face pair. Seams
     (one face on both sides) and zero-length curves are skipped.
   - STL: triangle edges whose two triangles have different face ids, grouped by face pair. These
     are exactly the ADR-005 face boundaries, so they use the same feature angle.
   - Cached with the faces under the same file-content key.
5. **Meshing does not change.** Physics derives an edge's node set from the existing surfaces: the
   nodes of the tri6 element edges (2 corners + midside) that occur in both `face:<a>` and
   `face:<b>`. Using shared element edges rather than intersecting node sets leaves out isolated
   points where two faces only touch at a corner. An empty result is an `InputError`. The support
   node set keeps the name `support:<id>`.
6. **API surface.**
   - Core: `edit.add_edge_support(project, part_id, edges, *, dofs, support_id)` and
     `describe_edges(project, part_id)`. Every pair is checked against `describe_edges`, so faces
     that do not meet are refused with a pointer to `tessera edges`.
   - CLI: `tessera edges PATH PART [--face N]`, and
     `tessera support add PATH PART --edge 3-5 [--edge 5-6] --dofs z`. Give face numbers or
     `--edge`, not both. `tessera info` prints `edges 3-5, 5-6`.
7. **Lines inside a face** are not modelled directly. The user adds a split line in CAD so the
   line becomes the edge between two faces.

## Options considered

### Option A1: geometry edge indices (gmsh curve tags; Tessera chains for STL)
| Dimension | Assessment |
|---|---|
| Complexity | Medium-high: new numbering for STL plus a gmsh-curve-to-edge mapping in `mesh` |
| STL | Feasible but costly; gmsh curves on rebuilt STL solids are not Tessera edges |
| Stability | A second numbering that can drift independently of face numbers |
| Selection UX | Long lists on filleted parts; users still describe an edge by its faces |

### Option A2: edges named by their adjacent face pair (chosen)
| Dimension | Assessment |
|---|---|
| Complexity | Low: no mesh change; one listing function per format; small physics helper |
| STEP / STL | Both, the same way |
| Stability | Exactly as stable as face numbers |
| Selection UX | Users and agents already know the faces (`tessera faces`) |

Limits: a pair sharing two curves cannot pick just one; seams cannot be chosen.

### Option B: no new type; split a narrow strip in CAD and use a face support
Ships today, but a strip also holds rotation (TS-08B becomes partly clamped), changes the span, and
STL users often cannot make the strip.

### Option C: node set by geometric predicate (points plus tolerance)
Independent of face numbering, but mesh-dependent, needs tolerance tuning, and silently selects the
wrong nodes on curved or coarse geometry.

## Trade-off analysis

Every option that works for STL ends up computing "the boundary between Tessera faces a and b". A2
uses that as the identity; A1 adds a numbering on top that geometry and mesh must reproduce in the
same way. A2 gives up selecting one of several curves shared by the same pair, and seams; no
catalogued benchmark or typical maker use needs either.

**Singular stresses.** A line support carrying a finite force per length is a stress singularity,
so the peak stress at the line grows without bound as the mesh is refined. Displacements away from
the line converge, so line-supported benchmarks gate on deflection. ADR-007 flags every line
support when reporting peak stress.

## Consequences

- Easier: simply supported beams and plates (TS-05, TS-08B, later TM-03); no mesh change.
- Easier: the representation can be reused for line loads (N/mm) later.
- Harder: every user of `support.faces` must handle `faces is None` (use `support.part_id`).
- Limits: one of several curves shared by two faces cannot be selected alone; lines inside a face
  need a CAD split line; seams are not selectable.

## Action items

1. [x] `model`: `EdgeSet`, `EdgeInfo`, `Support.faces | edges` (exactly one), `Support.part_id`,
       `SCHEMA_VERSION = 5`; schema-4 load test.
2. [x] `geometry`: `describe_edges` for STEP and STL, cached with faces.
3. [x] `physics`: edge node sets from shared tri6 element edges.
4. [x] `core`: `edit.add_edge_support`, `describe_edges`.
5. [x] `cli`: `tessera edges`, `support add --edge A-B`, `info` output.
6. [x] Benchmarks: TS-05 (supports on the bottom edges of the end faces) and TS-08B.
7. [x] Docs: API reference, user guide, catalog "Features" lines.
