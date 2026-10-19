# tessera.mcp

MCP server for AI agents.

- **Scope:** Thin MCP tool layer over tessera.core; no physics here.
- **May import:** tessera.core and tessera.model. Enforced by import-linter (see ADR-001).
- **Tests:** `packages/mcp/tests/`; run `uv run pytest packages/mcp`.
