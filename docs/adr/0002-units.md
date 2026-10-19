# ADR-002: Internal unit system

- Status: Accepted (maintainer, iteration 1)
- Date: 2026-10-19

## Decision

Use the consistent mm-N-MPa-tonne-s system for every quantity inside Tessera. Convert only at
the edges (file import, GUI, reports) with `tessera.model.units`.

## Why

- Printed parts are measured in mm; CAD and slicers export mm.
- mm-N-MPa is consistent (1 N/mm² = 1 MPa), so CalculiX decks need no conversion factors.
- Density in t/mm³ keeps modal analysis consistent later (PLA ≈ 1.24e-9 t/mm³).

## Consequences

Unit mistakes become an edge-only concern; any conversion inside a domain package is a review
finding.
