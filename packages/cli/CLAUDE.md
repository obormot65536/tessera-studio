# tessera.cli

Command-line interface.

- **Scope:** Thin Typer CLI over tessera.core; no physics here.
- **May import:** tessera.core and tessera.model. Enforced by import-linter (see ADR-001).
- **Tests:** `packages/cli/tests/`; run `uv run pytest packages/cli`.
