---
name: feature
description: Run Tessera Studio's development cycle for one issue - plan, test or benchmark first, implement, verify, review, PR. Use when starting work on any issue or feature request in this repo.
---

# /feature — one issue through the cycle

Input: an issue number or a one-paragraph feature description.

1. **Branch**: `git checkout -b feat/<issue>-<slug>` from an up-to-date `main`.
2. **Plan**: read CLAUDE.md, the relevant package CLAUDE.md files and `docs/adr/`. Draft a short
   plan: packages touched, public API changes, model (`tessera.model`) changes, tests.
   If the plan involves a real technical choice, write an ADR (`engineering:architecture`).
   Show the plan to the maintainer and wait for approval.
3. **Test first** (in parallel where independent):
   - physics behaviour → `fea-verifier` agent writes a benchmark under `benchmarks/cases/`
     and proposes reference + tolerance text (the maintainer writes `benchmarks/references/`);
   - other behaviour → `test-writer` agent writes failing tests in the package's `tests/`.
   Confirm the new tests fail for the right reason.
4. **Contract first**: if packages must agree on new types, add them to `tessera.model` and
   commit before splitting work.
5. **Implement**: one package per agent when packages are independent (e.g.
   `calculix-specialist` for `packages/solvers`); tell each agent exactly which folder it may
   write and which public API to code against. Main session integrates through `tessera.core`.
6. **Verify**: `uv run ruff check . && uv run ruff format --check . && uv run pyright &&
   uv run lint-imports && uv run pytest` (add `-m solver` / `benchmarks -m benchmark` when
   solver code changed).
7. **Review**: `/code-review`, `/simplify`, `/security-review` when file import or MCP changed,
   and the `reviewer` agent. Fix findings.
8. **PR**: Conventional Commit title; body lists what changed, numbers needing maintainer
   approval (materials, tolerances), and test results. The maintainer merges.
9. **After merge**: `docs-writer` updates docs and CHANGELOG.
