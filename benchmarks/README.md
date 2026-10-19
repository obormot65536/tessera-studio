# Benchmarks

Analytical validation cases that prove the solver pipeline is correct.

- `cases/`: one module per case (geometry, loads, closed-form solution with its source,
  mesh-convergence run). Written by the `fea-verifier` subagent.
- `references/`: reference values and tolerances. **Maintainer-owned**: a Claude Code hook
  blocks edits here; changes are proposed in the PR description.

Run: `uv run pytest benchmarks -m benchmark`.
