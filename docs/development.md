# Development cycle

The full plan (roles, subagents, skills, milestones, CI tiers, maintenance cadence) is the
"Tessera Studio — Development Plan & Cycle" document. This page is the short version that lives
with the code.

## Per issue

1. **Pick issue** (maintainer): one issue = one PR.
2. **Plan**: the `Plan` agent drafts; a real technical choice becomes an ADR in `docs/adr/`.
   The maintainer approves the plan.
3. **Test first**: `fea-verifier` writes the benchmark, `test-writer` the unit tests. They must
   fail before the code exists.
4. **Implement** on a branch; `calculix-specialist` for solver work; parallel git worktrees for
   independent packages.
5. **Verify**: hooks run ruff on every edit; then `uv run pytest`, `pyright`, `lint-imports`.
6. **Review**: `/code-review`, `/simplify`, `/security-review` (file import, MCP), `reviewer`.
7. **Approve and merge** (maintainer): reads the PR, signs off physics, squash-merges.
8. **Document**: `docs-writer` updates docs and CHANGELOG; close the issue.

## Milestones

| Milestone | Gate |
|---|---|
| v0.1 Foundations | Beam and plate benchmarks within 2 % in CI |
| v0.5 Print-aware + assemblies | Orthotropic patch tests exact; Hertz and bolted-joint tests pass; X/Y/Z coupons tested |
| v1.0 Studio MVP | 3 reference parts: failure spot matches, load within ±30 % |
