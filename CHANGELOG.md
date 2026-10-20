# Changelog

## Unreleased

- STL import (binary and ASCII): faces are split at sharp edges (40° default), parts must be
  a single closed solid, lengths are taken as mm (ADR-005).
- Directional (roller) supports; plate-with-hole benchmark; validation catalog.
- Repository scaffold: uv workspace with 12 packages, data model for assemblies, `tessera` CLI
  (`new`, `info`, `version`), Claude Code hooks and subagents, CI.
