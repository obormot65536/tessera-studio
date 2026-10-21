# How we know the results are right

Tessera is checked against simple problems whose answers are known from textbooks. If Tessera
gets those right, its mesher, supports, loads and solver are working together. This page
summarises the nine checks that pass today. The full plan, with 28 more checks still to do, is in
[`benchmarks/CATALOG.md`](../../benchmarks/CATALOG.md).

## What the checks do not show

- They do not check printed parts. No printed material is in the v0.1 library yet.
- They check simple shapes and loads. They do not prove that a complicated part is right.
- Passing means the answer is within a stated tolerance at one mesh size. Your own study still
  needs the two-mesh check described in [study-setup.md](study-setup.md), Step 7
  (`tessera solve --check-mesh` runs it for you).

## The nine passing checks

Each row is one check. "Gate" is the largest error allowed. "Error at gate" is the error Tessera
gave at the mesh size the test uses. Numbers come from one run of the suite on 21 October 2026
with CalculiX 2.21.

| Check | What it tests | Reference answer | Gate | Mesh size | Error at gate |
|---|---|---|---|---|---|
| **TS-01** Cantilever | Tip deflection of a steel bar, 200 × 10 × 10 mm, clamped at one end, 100 N pulling its free end down | Beam theory (Timoshenko): 1.5267 mm | 2 % | 2.5 mm | 0.39 % |
| **TS-02** Plate with a hole | Stress at the edge of a 10 mm hole in a 50 mm wide steel strip under tension | Howland's strip solution: stress concentration 2.512 | 2 % | 1.0 mm | 0.32 % |
| **TS-03** Patch test | Every node of a 10 mm steel cube under 10 MPa tension, on an irregular mesh: stress and displacement | The exact uniform field: stress 10 MPa, end stretch 4.762 × 10⁻⁴ mm | 2 × 10⁻⁵ (stress), 10⁻⁵ (displacement) | 3.0 mm | 1.0 × 10⁻⁵ (stress), 1.1 × 10⁻⁶ (displacement) |
| **TS-04** Bar under its own weight | Stretch of a 200 mm steel bar hanging under gravity, and stress halfway down | ρgL²/(2E) = 7.334 × 10⁻⁶ mm; stress ρgL/2 = 7.70 × 10⁻³ MPa | 2 % | 5.0 mm | 0.004 % (stretch), 0.001 % (stress) |
| **TS-05** Simply supported beam | Midspan deflection of a steel beam, 200 × 10 × 10 mm, resting on two line supports at its bottom end edges (one pinned, one a roller), under 0.5 MPa on top | Beam theory (Timoshenko): 0.5989 mm | 2 % | 2.5 mm | 0.14 % |
| **TS-06** Thick-walled cylinder | Hoop stress and radial growth of a steel tube (inner radius 10 mm, outer 20 mm) at 10 MPa inside | Lamé's solution: hoop stress 16.67 MPa, radial growth 9.365 × 10⁻⁴ mm | 2 % stress, 1 % growth | 1.25 mm | 0.08 % (stress), 0.003 % (growth) |
| **TS-07** Thick-walled sphere | Hoop stress and radial growth of a hollow steel ball (inner radius 10 mm, outer 20 mm) at 10 MPa inside | Lamé's solution: hoop stress 7.143 MPa, radial growth 3.810 × 10⁻⁴ mm | 2 % stress, 1 % growth | 1.25 mm | 0.48 % (stress), 0.009 % (growth) |
| **TS-08** Circular plate | Centre deflection of a steel disc (radius 20 mm, 2 mm thick) under 1 MPa, clamped at its rim (A) or resting on a ring (B) | Plate theory with a shear correction: 0.01699 mm (A), 0.06699 mm (B) | 2 % (A), 3 % (B) | 0.667 mm | A: 0.65 %. B: 0.22 % |
| **TS-09** Torsion | Twist and maximum shear stress of a 10 mm round shaft and a 10 mm square bar, 60 mm long, under 10 000 N·mm | Saint-Venant torsion theory | 2 % twist, 3 % shear | 2.5 mm | Round: 0.004 % (twist), 0.22 % (shear). Square: 0.42 % (twist), 0.78 % (shear) |

Notes on the rows:

- **TS-02** uses a value from Timoshenko and Goodier's *Theory of Elasticity* that the
  maintainer has approved but not yet re-checked against the book. The reference is marked
  provisional in the repository. Treat the 2 % gate as provisional too.
- **TS-03** checks every node, not one point. Its gate is near the precision CalculiX writes its
  results with (six significant digits), so it shows the mesh, loads and supports add no error of
  their own.
- **TS-04** holds the top face so it can slide sideways but not move in z. A full clamp would add
  a small sideways effect that has nothing to do with gravity.
- **TS-05** uses line supports, so the beam ends can turn. It models half the beam's width.
  The stress right at a line support never settles as the mesh gets finer, so the check
  compares deflection only.
- **TS-06** uses symmetry faces, so one quarter of the tube is modelled.
- **TS-07** uses three symmetry faces, so one eighth of the sphere is modelled. The pressure acts
  on a surface curved in both directions.
- **TS-08** models a quarter of the disc. The disc is thick enough that the reference includes a
  shear term (+4.6 %); without it the error would read 3.5 %. Variant B rests the disc on a ring
  line support at its bottom outer edge, which holds it only vertically. Its gate is wider
  because a ring support has no single exact 3D form; the measured gap is about 0.4 %.
- **TS-09** applies the torque with Tessera's moment load, the same load you can use in a study.

## Convergence

Each check also runs at a coarser and a finer mesh. The error must not grow as the mesh gets
finer, unless it is already well inside the gate. This shows the answer is settling on the
reference and not drifting.

## Run the checks yourself

The checks need CalculiX (`ccx`) on the PATH, or `TESSERA_CCX` set. From the repository root:

```
uv run pytest benchmarks -m benchmark -s
```

The `-s` flag prints the convergence tables. The run takes a few minutes. The continuous
integration pipeline runs the solver checks on Linux.

## Where the references live

Reference values and their tolerances are kept in `benchmarks/references/`. Only the maintainer
changes them. A change to a reference or tolerance is made in a pull request with an explicit
approval.
