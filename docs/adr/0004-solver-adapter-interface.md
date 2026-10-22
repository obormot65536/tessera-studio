# ADR-004: Solver adapter interface

**Status:** Accepted (maintainer, iteration 1)
**Deciders:** Maintainer

## Context

v1 solves with CalculiX; FEniCSx (thermal, CFD) and optionally Code_Aster come later. Physics
modules must not depend on a specific solver, and ADR-001 forbids domain packages from importing
each other. CalculiX and Gmsh are GPL; Tessera is GPL-3.0 too, but running external solvers as
separate processes also gives crash isolation (a segfault or a hung contact solve must not take
the GUI down) and lets a later licence change stay open.

## Decision

1. The interface lives in `tessera.model.solver` as solver-neutral, immutable types:
   - `SolveJob`: analysis kind, `Mesh` (nodes in mm; quadratic tetrahedra in the
     CalculiX/Abaqus C3D10 node order; named node sets, element sets and tri6 surfaces),
     materials, sections (element set → material), supports (node set + fixed DOFs) and
     nodal loads (N).
   - `SolveResult`: nodal displacements (mm), nodal stresses (MPa, xx yy zz xy yz zx), solver
     log and metadata.
   - `SolverAdapter` protocol: `name`, `capabilities() -> frozenset[AnalysisKind]`,
     `run(job, workdir) -> SolveResult`; raises `SolverError` carrying the solver log.
2. Converting user intent (a force on a face, a fixed face) into node sets and nodal loads is
   the physics package's job, so every adapter receives the same discrete problem.
3. Each adapter in `tessera.solvers` runs its solver as a **subprocess** in a per-job working
   directory, exchanging files (CalculiX: `.inp` in, `.frd`/`.dat` out). In-process Python
   solvers (FEniCSx) run in a worker process behind the same protocol.
4. `tessera.core` picks the adapter from an internal registry by analysis kind. A public plugin
   entry point is deferred (spec: no third-party plugin API in v1).

## Options considered

### Option A: Neutral job types + subprocess adapters — chosen
| Dimension | Assessment |
|---|---|
| Complexity | Medium |
| Isolation | Strong (process boundary) |
| Solver independence | High |
| Performance | File I/O overhead, negligible next to solve time |

### Option B: Solver-specific APIs called directly from physics
**Pros:** least code now. **Cons:** CalculiX keywords leak into physics; adding FEniCSx means
rewriting physics; violates ADR-001.

### Option C: Exchange through meshio/generic FE formats only
**Pros:** reuses a library. **Cons:** meshio carries mesh + fields, not BCs, materials,
orientations or contact; we would still need our own job types.

## Trade-off analysis

Option A costs one translation layer per solver but keeps physics, mesh and solvers testable
in isolation and protects the app from solver crashes. The discrete `SolveJob` is the natural
seam: it is exactly what a benchmark or a bug report needs to reproduce a run.

## Consequences

- Easier: adding solvers; reproducing failures (`SolveJob` + workdir is a complete repro).
- Harder: orthotropic orientations, contact pairs and couplings must be expressed neutrally;
  the job types grow with v0.5 (orientations, connections) — extend them additively.
- Revisit: stress recovery point (nodal vs integration point) when failure criteria land.

## Action items

1. [ ] Add `Mesh`, `SolveJob`, `SolveResult`, `SolverAdapter`, `SolverError` to `tessera.model`.
2. [ ] Implement `CalculixAdapter` (static, isotropic) with a patch test.
3. [ ] Extend additively in v0.5: orientations, orthotropic materials, ties, contact, couplings.
4. [ ] Registry by analysis kind once a second adapter exists (v0.1 core uses CalculiX
       directly but checks `capabilities()`).
5. [x] Default solver timeout before the GUI/MCP server can start solves (a hung contact solve
       must not block them forever). Done: `tessera.core.solve` gives the built-in solver
       `SOLVE_TIMEOUT` = 1800 s per solver run, overridable with `solve(timeout=...)` and
       `tessera solve --timeout SECONDS`. The limit is per run, so `check_mesh` (two solves) can
       take up to twice as long; meshing is not counted. On a timeout the adapter kills ccx
       and the processes it started (its own process group on POSIX, `taskkill /F /T` on
       Windows) and raises `SolverTimeoutError`, a `SolverError`; with `check_mesh` the
       message names the coarse or fine solve. An exception while waiting (Ctrl+C, a cancel)
       kills them the same way and re-raises that exception. On POSIX, SIGTERM or SIGHUP to
       a main thread that is waiting kills them before the signal takes effect, and on Linux
       `PR_SET_PDEATHSIG` kills the direct ccx process if Tessera dies in any way. Not
       covered: a process that leaves the group itself (`setsid`), and SIGTERM, SIGHUP or
       SIGKILL on macOS while a non-main thread waits (or SIGKILL at all on macOS/Windows).
       1800 s leaves room for a large fine mesh of the two-mesh check on a slow laptop, where
       killing a solve that would have finished costs more than waiting for a rare hang; the
       GUI should still offer a cancel button.
6. [ ] High-order mesh quality: enable gmsh high-order optimization or report inverted curved
       elements clearly before CalculiX rejects them with a Jacobian error.
7. [ ] Recover stresses at integration points for failure criteria: CalculiX's nodal frd
       stresses carry about 1e-5 relative extrapolation bias.
