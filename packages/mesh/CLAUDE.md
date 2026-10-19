# tessera.mesh

Gmsh-based tetrahedral meshing.

- **Scope:** Tet meshes with physical groups per part/region and per-element material orientation.
- **May import:** tessera.model. Enforced by import-linter (see ADR-001).
- **Tests:** `packages/mesh/tests/`; run `uv run pytest packages/mesh`.
