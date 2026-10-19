# tessera.model

Shared data model, units and interfaces (bottom layer).

- **Scope:** Project, Assembly, PartInstance, Connection, unit system, solver-adapter and physics protocols.
- **May import:** nothing else in tessera. Enforced by import-linter (see ADR-001).
- **Tests:** `packages/model/tests/`; run `uv run pytest packages/model`.
