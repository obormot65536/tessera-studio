# Tessera Studio

Will your 3D-printed part hold? Tessera Studio is an open-source (GPL-3.0) desktop app that
analyzes 3D-printed parts and assemblies for stress and strain, taking into account that FDM
prints are much weaker between layers than along them.

**Status:** early development (v0.1 Foundations). Nothing usable yet.

## Planned v1

- Import parts and assemblies (STEP, STL, 3MF) and print settings from PrusaSlicer, OrcaSlicer
  or Bambu Studio
- Region-wise orthotropic material model for FDM (walls, skins, infill); isotropic resin
- Bonded joints, frictional contact, bolts and pins between parts
- Safety-factor map that says whether a part fails within a layer or between layers
- Desktop GUI with Beginner and Expert modes, a CLI, and an MCP server for AI assistants

## Development

```bash
uv sync
uv run pytest
uv run tessera --help
```

Solver tests need CalculiX (`ccx`) on `PATH`. See `CLAUDE.md` for the package map and rules,
and `docs/development.md` for the development cycle.

## License

GPL-3.0-or-later. See `LICENSE`.
