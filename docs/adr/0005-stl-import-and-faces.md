# ADR-005: STL import and face definition

**Status:** Accepted (maintainer, iteration 2)
**Deciders:** Maintainer

## Context

Makers mostly have STL files: a triangle soup with no units, no solid and no faces. Supports and
loads are placed on faces, and the face numbers a user or agent sees (`describe_faces`) must be
exactly the faces the mesher loads (`face:<i>` surfaces). gmsh can rebuild and remesh a solid
from an STL (`classifySurfaces` + `createGeometry`), but to reparametrize it cuts the surface
into patches: a box becomes 10 patches, a plate with a hole 20. Those patches are not faces a
person would name.

## Decision

1. Tessera defines STL faces itself: connected groups of triangles separated by sharp edges,
   where the dihedral angle between neighbours exceeds a feature angle (default 40°).
2. `tessera.geometry.read_stl` returns a new model type, `TriangleSurface` (vertices, outward
   triangles, face id per triangle). It checks the part is a single closed solid and fixes the
   winding so normals point out.
3. `tessera.mesh.mesh_surface` takes that `TriangleSurface`, lets gmsh rebuild and remesh the
   solid, and maps every gmsh patch back to a Tessera face through the original triangles.
4. `tessera.core` dispatches on the file extension (`.step`/`.stp` vs `.stl`); clients see no
   difference.
5. STL lengths are taken as millimetres, the 3D-printing convention (ADR-002: a unit option
   can later be added at the import edge).

## Options considered

### Option A: Tessera-defined faces passed as a model type (chosen)
**Pros:** one definition of faces, shared by describe and mesh through `tessera.model`; stable
face numbers independent of gmsh internals; respects ADR-001.
**Cons:** own feature-edge segmentation to maintain; coarse STLs (e.g. an 8-facet cylinder,
45° between facets) split a curved wall into several faces.

### Option B: expose gmsh's reparametrization patches as faces
**Pros:** no own segmentation. **Cons:** faces nobody would name (a box face split in two), and
geometry and mesh would each have to reproduce gmsh's tagging; one option drift away from
loads landing on the wrong face.

### Option C: convert STL to a STEP/B-rep solid first
**Pros:** reuses the STEP path. **Cons:** one B-rep face per triangle (thousands of faces), slow,
and the same face-naming problem.

## Consequences

- Easier: STL, STEP and later 3MF meshes share the same face contract.
- Harder: feature angle is a heuristic; expose it as a parameter and report face counts so users
  notice odd splits.
- Face numbers follow the triangle order in the file: stable for one file, but re-exporting
  the STL can renumber faces and move saved supports/loads. A later release should re-locate
  saved faces by centroid/normal when the geometry file changes.
- `describe_faces` and meshing must use the same feature angle. Both use the default today;
  when the angle becomes a project setting, it is passed to both.
- Feature angles below 1° are rejected: near-flat neighbours would become face boundaries
  that floating-point noise cannot separate reliably.
- Triangles with a repeated corner are dropped; zero-area triangles with three distinct
  collinear corners are rejected with a repair hint (dropping them would open a crack).
- Revisit: a unit option for STL and a bounding-box sanity warning; healing of slightly open
  or sliver-containing STLs (currently rejected); internal cavities (currently rejected as
  disconnected shells).
