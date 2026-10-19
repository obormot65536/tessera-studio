# tessera.printsetup

Slicer project import and print-region tagging.

- **Scope:** Parse PrusaSlicer/OrcaSlicer/Bambu 3MF/INI; derive wall/skin/infill regions.
- **May import:** tessera.model. Enforced by import-linter (see ADR-001).
- **Tests:** `packages/printsetup/tests/`; run `uv run pytest packages/printsetup`.
