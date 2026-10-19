# Tessera Studio

GPL-3.0 desktop app that tells makers whether a 3D-printed part or assembly will hold, using
layer-wise orthotropic material models for FDM prints (isotropic for resin). Python-first:
PySide6 GUI, Gmsh meshing, CalculiX solver. Solo maintainer + Claude.

## Package map (`packages/`, all under the `tessera.` namespace)

| Layer | Packages | Role |
|---|---|---|
| Clients | `gui`, `cli`, `mcp` | Thin front ends over `tessera.core` |
| Core | `core` | Public Project API, orchestration, save/load, job runner |
| Domain | `geometry`, `printsetup`, `materials`, `mesh`, `physics`, `solvers`, `post` | Independent siblings |
| Model | `model` | Data model, units, interfaces shared by everyone |

Each package has its own `CLAUDE.md` with scope and allowed imports.

## Commands

- Install: `uv sync`
- Tests: `uv run pytest` (one package: `uv run pytest packages/<name>`; solver tests need `ccx`)
- Lint: `uv run ruff check . && uv run ruff format --check .`
- Types: `uv run pyright`
- Layering: `uv run lint-imports`
- CLI: `uv run tessera --help`

## Hard rules

1. Layering (ADR-001): a package imports only from layers below it; domain packages never import
   each other. `lint-imports` enforces this in CI.
2. Units (ADR-002): mm, N, MPa, tonne, s everywhere inside Tessera; convert only at the edges.
3. Never edit `benchmarks/references/` (reference values and tolerances). A hook blocks it.
   Propose changes in the PR description; the maintainer applies them.
4. Every behaviour change lands with a test or benchmark that failed before it.
5. GUI, CLI and MCP contain no physics; they call `tessera.core`.
6. Any change to a number a user will trust (material data, failure criteria, tolerances) needs
   the maintainer's explicit approval in the PR.

## Workflow

Follow the development cycle in `docs/development.md`: pick issue → plan (ADR if a real choice)
→ test/benchmark first → implement → verify → review (`/code-review`, `/simplify`,
`/security-review`, `reviewer` agent) → maintainer merges → `docs-writer`.

Work on a branch, one issue per PR, Conventional Commits. Never push to `main` directly.

Subagents live in `.claude/agents/`; use `fea-verifier` before physics work and
`calculix-specialist` for anything solver-related.
