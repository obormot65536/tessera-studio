# ADR-001: Package layering

- Status: Accepted (maintainer, iteration 1)
- Date: 2026-10-19

## Context

The spec puts "project model, units, selection, job runner" in `tessera.core`, and has
`tessera.core` orchestrate the domain packages (geometry → printsetup → materials → mesh →
physics → solvers → post). The domain packages also need the data model, so keeping it in
`core` would create import cycles between `core` and every domain package.

## Decision

Split the shared data model out of `core` into a bottom-layer package `tessera.model`:

1. `tessera.gui | tessera.cli | tessera.mcp`: clients, independent of each other
2. `tessera.core`: Project API and orchestration
3. `tessera.geometry | printsetup | materials | mesh | physics | solvers | post`: independent
   domain siblings
4. `tessera.model`: data model, units, solver-adapter and physics interfaces

Domain packages never import each other; `core` passes `model` objects between them.
`import-linter` enforces this in CI.

## Consequences

- No import cycles; each domain package can be tested alone.
- Interfaces (e.g. the solver adapter protocol) live in `model`, so a new solver or physics
  module touches one domain package plus, at most, `model`.
- Small change to the spec's architecture diagram: "project model" moves from `core` to `model`.
