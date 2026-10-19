# tessera.solvers

Solver adapters (CalculiX first).

- **Scope:** Write solver input, run the external process, read results back into model types.
- **May import:** tessera.model. Enforced by import-linter (see ADR-001).
- **Tests:** `packages/solvers/tests/`; run `uv run pytest packages/solvers`.
