---
name: add-benchmark
description: Add an analytical validation benchmark to Tessera Studio (closed-form reference, tolerance file, mesh convergence). Use when a physics feature needs a validation case or the user asks for a new benchmark.
---

# /add-benchmark

1. Choose a problem with a closed-form or well-published solution; note its source.
2. Delegate to the `fea-verifier` agent (fallback: a general-purpose agent told to follow
   `.claude/agents/fea-verifier.md`). It writes `benchmarks/cases/test_<name>.py`:
   - marked `@pytest.mark.benchmark` (and `@pytest.mark.solver` if it runs CalculiX);
   - geometry built with the gmsh API, problem set up only through the public `tessera.core` API;
   - analytical value computed in the test and checked against the reference file to 1e-9;
   - reference + tolerance read from `benchmarks/references/<name>.toml` with tomllib,
     `pytest.skip("reference awaiting maintainer approval")` when the file is missing;
   - a mesh-convergence test over at least 3 sizes.
3. Present the proposed `benchmarks/references/<name>.toml` content (description, source,
   reference value, rel_tolerance with reasoning) to the maintainer. Only the maintainer writes
   it, or approves it explicitly; a hook blocks edits there otherwise.
4. Run `uv run pytest benchmarks -m benchmark` and report the error and convergence table.
