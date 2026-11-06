# Strength ratios

A strength ratio tells you how far a printed part is from its first failure. It is the
number of times you could multiply every load on the part, gravity included, before the part
is predicted to fail at its weakest spot. The part's other stress results do not say that.
A stress is only a number; a ratio is compared with the loads you applied.

This page explains the lines `tessera solve` prints, how to give a material its strengths,
and what the check does not cover.

Strength ratios only appear when the material of the part has strengths. One built-in material
has strengths, `pla-fff-gonabadi2020` (see "Built-in printed materials" below). For any other
material, you enter the strengths yourself, from tests on your own printed pieces. See "Give a
material its strengths" below.

## Read this first: what the ratio does and does not cover

**The ratio is a first-ply onset check only.** It says where the first failure is predicted
for the loads you applied, and how far the loads are from it. It does not model:

- progressive damage, where the part keeps carrying load after the first failure;
- fatigue (repeated loading);
- creep (slow change under a steady load over time);
- buckling (a slender part giving way under compression);
- print defects, beyond what your coupon strengths already contain.

A printed part's strength depends on how it was printed, and the coupon strengths you enter
are the only place that is taken into account. If a print has voids, weak layer bonding or
poor infill that your test pieces did not have, the ratio does not see it. Use your own
judgement about the margin you need. Tessera does not pick a safety factor.

The check also has no benchmark yet. The check of the criterion against exact cases (benchmark
TO-06) is still to come, so treat the numbers as not yet validated.

## What the ratio means

Here is example output from `tessera solve`. The part is a 40 x 10 x 10 mm cantilever. The
strengths are illustrative, made-up values, not data for any real material. This is a first-ply
onset check only: it shows where the first failure is predicted, and it does not model progressive
damage, fatigue, creep, buckling or print defects beyond what the coupon strengths contain.

```
lowest strength ratio 1.91 (Tsai-Wu, along the raster, tension), away from supports and sharp inside edges, at (5.862, 4.991, 10) mm, face 6
  The loads could be about 1.9 times higher before the first failure is predicted there.
lowest strength ratio 1.36 at (1.797, 8.224, 0) mm, face 5: near support 'base' (within 5.81 mm); the ratio there may keep falling as the mesh is refined; compare two mesh sizes (--check-mesh) before trusting either value
assumed: compressive strengths c1, c2, c3 equal to the tensile ones (material 'petg-x')
```

Read it like this:

- A ratio of 1.91 means the loads could be about 1.9 times higher before first failure is
  predicted at that spot. A ratio of 1 means the applied loads are right at the predicted
  first failure.
- The ratio is lowest at the weakest spot, so it is the number to look at first.
- The ratio at the weakest spot tells you about one spot only. It does not say how the rest of
  the part behaves beyond that point.

## The lines, word for word

`tessera solve` prints the strength lines after the displacement line and before the von
Mises lines. The exact wording is below, with what each line means.

**The headline ratio.** The lowest ratio away from supports and sharp inside edges:

```
lowest strength ratio 1.91 (Tsai-Wu, along the raster, tension), away from supports and sharp inside edges, at (5.862, 4.991, 10) mm, face 6
  The loads could be about 1.9 times higher before the first failure is predicted there.
```

- The bracket holds the failure mode: the criterion and the way the part is predicted to fail
  (see "Failure modes" below).
- The location is in mm and the face is the face number from `tessera faces`.
- The second line says what the number means. Ratios are rounded: two significant figures in
  this line, three in the first.

**A ratio below 1.** Failure is predicted at the applied loads:

```
  First failure is predicted there at these loads: it starts at about 0.62 times them.
```

**A flagged spot.** The lowest ratio is near a support or a sharp inside edge. Tessera shows it
but does not make it the headline, because the stress there depends on the mesh:

```
lowest strength ratio 1.36 at (1.797, 8.224, 0) mm, face 5: near support 'base' (within 5.81 mm); the ratio there may keep falling as the mesh is refined; compare two mesh sizes (--check-mesh) before trusting either value
```

- "near support 'base'" names the support. "within 5.81 mm" is the flag distance (twice the
  median element edge): every spot this close to a support or sharp inside edge is flagged.
- "near a sharp inside edge" is used instead when the spot is near a sharp inside edge.
- Run `tessera solve --check-mesh` to compare two mesh sizes. The line then ends with one of
  the verdicts below.

**The two-mesh verdicts.** With `--check-mesh`, the flagged spot ends with one of:

```
lowest strength ratio 1.62 (Tsai-Wu, interlayer shear) at (0, 2, 1) mm, face 3: near support 'base' (within 5 mm); settled, changed 0.62 % between the meshes
lowest strength ratio 1.62 at (0, 2, 1) mm, face 3: near support 'base' (within 5 mm); mesh-dependent, changed 23.5 % between the meshes
```

- `settled` means the ratio changed by 5 % or less between the two meshes. The spot can then
  lead the result, and the line gives its failure mode.
- `mesh-dependent` means it changed by more than 5 %. Do not trust that number. Use a finer
  mesh and check again.
- `not judged` means the check could not compare the value. The text says why. A common reason
  is that the second mesh was not finer around the spot, or the wall there is too thin for the
  check to see the change (the line then gives a mesh size to try).

**No clear headline.** Every checked spot is near a support or a sharp inside edge:

```
lowest strength ratio away from supports and sharp inside edges: none, every checked node is near a marked support or sharp inside edge
```

**No stress.** No checked spot carries stress, so the ratio is infinite and is not printed as a
number:

```
lowest strength ratio: none, no checked node carries stress
```

**Assumed strengths.** A compressive strength you did not enter is taken equal to the tensile
one. Each one is listed, for each material:

```
assumed: compressive strengths c1, c2, c3 equal to the tensile ones (material 'petg-x')
assumed: compressive strength c3 equal to the tensile one (material 'petg-x')
assumed: compressive strength equal to the tensile one (material 'resin')
```

Read these lines as a warning. Each assumed value is a guess, and it is only as good as the
tensile value it copies. Enter the real compressive strength if you measured it.

**Unchecked materials.** A material with no strengths has its elements left out of the check:

```
not checked for strength (no strength data): material 'pla'
```

A one-material part with no strengths prints no strength lines at all; see "When the
strength lines appear" below.

Ratios of 1000 or more are printed as whole numbers, for example `12346`.

### When the strength lines appear

- Only when the part's material has strengths: either a built-in one that carries them
  (`pla-fff-gonabadi2020`), or a `strength` in the project file. Then the strength lines
  replace the orthotropic note. The note says that von Mises stress is not a strength check for
  a layer-wise part. A material without strengths still gets the note when it is orthotropic.
- From Python, the same results are in `StudyResult.strength`. See `docs/api/core.md`.

## Failure modes

The mode says which stress is predicted to fail first. It is a label only: the ratio does not
change with it. It is the stress term with the largest stress relative to its strength at the
weakest spot. The mode tells you which way the part is predicted to break.

| Line says | What it means |
|---|---|
| along the raster, tension | The roads are pulled apart along their length. |
| along the raster, compression | The roads are crushed along their length. |
| across the roads, tension | Neighbouring roads in a layer pull apart. |
| across the roads, compression | Roads in a layer are pushed together. |
| across the layers, tension | The layers pull apart. Printed parts are often weakest in this direction. |
| across the layers, compression | The layers are crushed together. |
| in-layer shear | Roads slide over each other within a layer. |
| interlayer shear | Layers slide over each other. |
| von Mises | A ductile (metal) part reaches its yield strength. |
| maximum principal stress, tension | A brittle (resin) part is pulled apart. |
| maximum principal stress, compression | A brittle (resin) part is crushed. |

## Built-in printed materials

Two printed materials are built in. Each id names the study its numbers come from, because the
same plastic measures very differently on different printers and settings. Every number, its
source and its assumptions are listed in `docs/materials/pla-petg.md`.

| Id | What it is | Strengths |
|---|---|---|
| `pla-fff-gonabadi2020` | PLA printed on a desktop Ender-3 at 0.15 mm layers, 100% infill. A few values come from a second PLA study, and the interlayer shear plane is assumed. | Yes. Across the layers it is only 5 MPa in tension, so parts pulled across their layers show low ratios. |
| `petg-pellet-romeijn2022` | PETG from a large-format pellet printer (5 mm roads, 2 mm layers). Not desktop filament: across the roads it is likely less stiff than your prints. | No. No published study gives PETG strengths across the layers, so a PETG part gets the orthotropic note instead of a strength ratio: see "When the strength lines appear" above. |

These are values from one set of test pieces each, not guarantees for your printer, filament
or settings. If you can test your own prints, define your own material in the project file
(below). Give it its own id, for example `my-pla`, and set your part to it.

## Give a material its strengths

Strengths go in the project file, on a material in the top-level `materials` list. They are
in MPa, and every value must be above zero. The list and the stiffness values are described in
`docs/user/project-files.md`. The examples below are made-up values, not data for any real
PLA, PETG or resin:

```json
"materials": [
  {"kind": "orthotropic", "id": "my-pla", "...": "stiffness as in project-files.md",
   "strength": {"t1": 40.0, "t2": 30.0, "t3": 15.0, "s12": 20.0, "s13": 12.0, "s23": 12.0,
                "c1": 50.0}},
  {"kind": "isotropic", "id": "my-resin", "...": "stiffness as in project-files.md",
   "strength": {"tensile": 50.0, "compressive": 80.0, "brittle": true}}
]
```

A part names its material by id, as usual (`material_id`). An orthotropic material also needs a
print setup on the part, or the solve is refused.

### Orthotropic (printed in layers)

Use this for FDM printed parts, with the axes described in `project-files.md`:

| Field | Meaning | Required |
|---|---|---|
| `t1`, `t2`, `t3` | Tensile strength along the raster (1), across the roads (2), and across the layers (3). | Yes |
| `s12` | Shear strength within a layer. | Yes |
| `s13`, `s23` | Shear strengths between layers. | Yes |
| `c1`, `c2`, `c3` | Compressive strengths on the same axes. | No |

- A missing compressive strength is taken equal to the tensile one on that axis. It is then
  listed in the `assumed` line.
- The shear strengths have no default. You must enter them.
- The criterion is Tsai-Wu. Tessera uses a fixed interaction term, so you do not enter one. If
  every tensile and compressive strength is equal, and each shear strength is that value divided
  by the square root of 3, it gives the same answer as von Mises.

### Isotropic (the same in every direction)

Use this for metals and for resin:

| Field | Meaning | Required |
|---|---|---|
| `tensile` | Yield strength of a ductile material (a metal), or the tensile strength of a brittle one. | Yes |
| `brittle` | `true` for a brittle material such as a printing resin. The default is `false`. | No |
| `compressive` | Compressive strength. Only for a brittle material. | No |

- A ductile material (`"brittle": false`) is checked with von Mises against its yield strength.
  Tessera refuses a `compressive` value for it, because nothing would use it.
- A brittle material (`"brittle": true`) is checked on its largest tensile and its largest
  compressive stress. Its missing `compressive` is taken equal to `tensile`, and listed.

### Where to get the numbers

Enter the strengths you measured on your own printed test pieces, printed the way the part will
be printed. Of the built-in materials only `pla-fff-gonabadi2020` has strengths, from one
published study; your own tested values describe your prints better. Tessera does not check
that your numbers are right. A wrong strength gives a wrong ratio, and the ratio looks just
as confident.

Strengths you enter are data you own. They are not checked against a reference.

## Solve time

A part with a strength check reads each element's own stresses from the solver output. The
measurement in ADR-008 gave about 8 s of reading at one million elements. A part without
strengths is solved as before.

## See also

- [Project files](project-files.md): the materials list, print setup and file format.
- [Set up a study](study-setup.md): step 9 explains the von Mises lines and the mesh check.
- [API reference](../api/core.md): `StudyResult.strength` and the strength types.
