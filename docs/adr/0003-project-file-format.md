# ADR-003: Project file format

**Status:** Accepted (maintainer, iteration 1)
**Deciders:** Maintainer

## Context

A project holds an assembly (parts, placements, materials, print setups, connections), a study
(supports, loads, analysis settings) and, once solved, a mesh and results with 10^5–10^6 nodes.
Users share projects with friends and attach them to bug reports, so a project must be one
self-contained file that still opens when the original STEP files have moved. GPL-3.0 and a
Python-first stack favour formats readable without extra native dependencies.

## Decision

A `.tessera` file is a **zip container**:

```
project.json        # Project model (pydantic), schema_version, relative paths only
geometry/<part>.step  # copies of imported geometry (original path kept as metadata)
mesh/<study>.npz      # nodes, elements, sets (numpy, optional)
results/<study>.npz   # displacements, stresses, safety factors (optional)
```

`project.json` is the only required entry. Mesh and results are caches: a file without them
still opens and re-solves. A `schema_version` integer plus a migration function per version
handles format changes from v0.5 on.

During v0.1, a plain `*.tessera.json` (the `project.json` part alone, geometry by path) is the
development format; the zip container lands with result storage in v0.5.

## Options considered

### Option A: Single JSON file, geometry by path
| Dimension | Assessment |
|---|---|
| Complexity | Low |
| Self-contained | No: breaks when STEP files move |
| Large arrays | Poor (text, huge) |
| Dependencies | None |

**Pros:** human-readable, diff-friendly. **Cons:** not shareable; no place for mesh/results.

### Option B: Zip container (JSON + STEP + NPZ) — chosen
| Dimension | Assessment |
|---|---|
| Complexity | Low–medium |
| Self-contained | Yes |
| Large arrays | Good (binary NPZ, compressed) |
| Dependencies | stdlib `zipfile` + numpy |

**Pros:** one shareable file; JSON part stays readable after unzip; caches optional.
**Cons:** whole-file rewrite on save (fine at our sizes).

### Option C: HDF5
**Pros:** partial I/O for huge results. **Cons:** h5py native dependency, opaque to users,
overkill below ~10^7 nodes.

### Option D: SQLite
**Pros:** transactional, partial updates. **Cons:** arrays as blobs anyway; poor fit for a
document-style model.

## Trade-off analysis

Shareability and robustness to moved files outweigh partial I/O. Results are recomputable, so
they are treated as caches, which keeps migrations limited to `project.json`.

## Consequences

- Easier: sharing, bug reports, MCP agents receiving one file.
- Harder: very large result sets (>1 GB) load slowly; revisit (HDF5 for `results/`) if that
  becomes real.
- Must write a migration for every `schema_version` bump from v0.5 on.

## Amendment (iteration 4): container brought forward to v0.1

STL faces are numbered by Tessera from the file's triangles (ADR-005), so saved supports and
loads are only safe while the exact geometry bytes travel with the project. The container
therefore lands in v0.1, without the mesh/results caches:

- `project.json` is the `Project` model; each part's `geometry` is `geometry/<name>`, relative
  to the container root. Parts that share one source file share one entry.
- `*.tessera.json` stays as the hand-editable form; relative geometry paths resolve against the
  JSON file's folder. An unzipped `.tessera` is therefore a valid editable project
  (`project.json` + `geometry/`).
- Loading reads the file by content (a zip is a container, anything else is JSON), copies the
  geometry into a private temporary folder that lives until the process exits, and never
  extracts by member name: only `geometry/<name>` entries with a plain name are read, and each
  is size-capped before decompression.
- The original import path is **not** stored: shared projects and bug reports would leak the
  sender's folder layout and user name.
- Saving writes a temporary file next to the target and renames it over the target, so a failed
  save never leaves a half-written project.

## Action items

1. [x] v0.1: keep `*.tessera.json`; store geometry paths relative to the project file.
2. [x] v0.1: zip container with `project.json` and geometry (amendment above).
3. [ ] v0.5: mesh/results caches with a staleness check, and the first migration test.
