# tessera.core

Project API and job orchestration used by every client.

- **Scope:** The public Python API (`tessera.core.Project`), selection, save/load, job runner.
- **May import:** tessera.model and all domain packages. Enforced by import-linter (see ADR-001).
- **Tests:** `packages/core/tests/`; run `uv run pytest packages/core`.
