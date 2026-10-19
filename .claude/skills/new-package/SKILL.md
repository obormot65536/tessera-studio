---
name: new-package
description: Create a new workspace package in Tessera Studio with pyproject, CLAUDE.md, tests and import-layer registration. Use when a new tessera.* package is needed.
---

# /new-package <name> "<one-line description>" <layer>

`<layer>` is one of: client, core, domain, model (see ADR-001).

1. Create `packages/<name>/` with:
   - `pyproject.toml`: name `tessera-<name>`, `requires-python = ">=3.12"`,
     `license = "GPL-3.0-or-later"`, hatchling build with `packages = ["src/tessera"]`,
     workspace sources for every `tessera-*` dependency;
   - `src/tessera/<name>/__init__.py` (docstring = the description) and `py.typed`;
     **no** `src/tessera/__init__.py` (namespace package);
   - `CLAUDE.md` with Scope, May import, Tests (copy the format of an existing package);
   - `tests/test_<name>_import.py` smoke test.
2. Root `pyproject.toml`: add the dependency and its `[tool.uv.sources]` entry, add
   `tessera.<name>` to `[tool.importlinter] root_packages` and to the right layer in the
   layers contract (siblings separated by ` | `).
3. Update the package map table in the root CLAUDE.md.
4. `uv sync && uv run lint-imports && uv run pytest packages/<name>`.
5. If the layer placement is not obvious, write or amend an ADR.
