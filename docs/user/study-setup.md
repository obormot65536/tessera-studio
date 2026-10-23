# Set up a study for one part

This guide walks through a complete study of one 3D-printed or machined part, from the file to
the numbers. It is written for makers. Every step uses units of **mm, N, MPa, and mm/s²**.

The steps use the `tessera` command line. You can also edit the project file by hand; the last
section shows that form.

## What version 0.1 can and cannot do

Read this first. It saves you from setting up a study that Tessera cannot run.

- **One part per project.** No assemblies. No connections between parts. A part cannot be
  moved or rotated (its placement must stay at the identity). The command line lets you add a
  second part, but the solve refuses a project with more than one part.
- **Linear static analysis only.** The model assumes small displacements and a material that
  springs back when unloaded. It does not model buckling, large bending, plastic yield,
  contact, heat, vibration or time.
- **A short material library, and your own values.** Two library materials exist:
  `steel-generic` and `aluminium-6061-t6`. There is no printed plastic (PLA, PETG, ABS, nylon)
  in the library yet. You can define your own material in the project, including a layer-wise
  (orthotropic) one from your own coupon tests; see "Your own printed material" in Step 3.
  There is no command for that yet: use the Python API or edit the JSON.
- **Displacement and von Mises stress come out.** The largest displacement, and the largest
  von Mises stress with where it is: away from supports and sharp inside edges, anywhere, and
  on each face. Tessera does not compare them with a limit, does not compute a safety factor
  and does not say whether the part will hold.
- **Needs CalculiX.** The solver `ccx` must be installed and on the PATH, or its location set
  in the `TESSERA_CCX` environment variable.

## Step 1: Get the part as a file

Export the part from your CAD tool as:

- **STEP** (`.step` or `.stp`): one solid. Export in millimetres. Tessera does not rescale
  the geometry.
- **STL** (`.stl`): one closed body. Every edge must be shared by two triangles; open edges
  are refused. Tessera splits the surface into faces wherever two neighbouring triangles meet
  at an angle of more than 40°. Lengths are read as mm.

Put the geometry file in its own folder with the project file. The example below uses
`beam.step`, a 200 × 10 × 10 mm steel-like bar. You can make the same bar in any CAD tool to
follow along.

## Step 2: Create the project and add the part

```
tessera new beam.tessera.json
tessera part add beam.tessera.json beam.step --id beam --name Beam
```

Output:

```
created beam.tessera.json
added part 'beam'
```

- `new` makes an empty project. A name ending in `.tessera.json` gives a plain JSON file that
  you can open in a text editor. A name ending in `.tessera` gives one file that holds the
  project and a copy of the geometry. The commands in this guide work with both.
- `part add` takes the project, then the geometry file.
- `--id` is the short name you use in later commands. If you leave it out, the id is the file
  name.
- `--name` is a label for people. If you leave it out, the label is the id.
- `--material` sets the material in the same step. Step 3 does it separately.

The geometry path is stored relative to the project file. You can run the commands from any
folder.

## Step 3: Choose the material

List the material ids:

```
tessera materials
```

Output:

```
aluminium-6061-t6
steel-generic
```

Set the material of the part:

```
tessera part material beam.tessera.json beam steel-generic
```

Output:

```
part 'beam' is steel-generic
```

| Material id | Young's modulus E (MPa) | Poisson's ratio ν | Density (t/mm³) |
|---|---|---|---|
| `steel-generic` | 210 000 | 0.30 | 7.85e-9 |
| `aluminium-6061-t6` | 68 900 | 0.33 | 2.70e-9 |

Density is needed only if you use gravity. An unknown id is refused, and the message lists the
known ids.

If your part is a printed plastic, these values are not right for it. Use the closest library
entry only to get a rough idea, and do not treat the numbers as a print result.

### Your own printed material

A print is stiffer along its rasters than across them, and stiffer within a layer than across
layers. Tessera models this as an orthotropic material with three axes:

- axis 1: along the rasters (the extrusion direction);
- axis 2: across the rasters, in the layer;
- axis 3: the build direction, across the layers.

It needs nine numbers in MPa: `e1`, `e2`, `e3` (Young's moduli), `g12`, `g13`, `g23` (shear
moduli) and `nu12`, `nu13`, `nu23` (Poisson's ratios: `nu12` is the contraction along axis 2
when you pull along axis 1). Tessera does not supply these numbers: they come from your own
coupon tests or a paper you trust, and Tessera only checks that they are physically possible.

The part also needs a **print setup**, so Tessera knows which way the layers run:

- `build_direction`: the direction, in the part's CAD coordinates, that points up on the printer.
  `[0, 0, 1]` means the part is printed as modelled; `[1, 0, 0]` means it lies on its side with
  its +x end up.
- `raster_angle`: the raster direction in degrees, turned counter-clockwise about the build
  direction (0 by default).

A part with an orthotropic material and no print setup is refused when you solve, rather than
guessing the orientation. There is no command for this yet. With the Python API:

```python
from pathlib import Path

from tessera.core import edit, load_project, save_project
from tessera.model import OrthotropicDefinition

path = Path("beam.tessera.json")
project = load_project(path)
edit.add_material(
    project,
    OrthotropicDefinition(
        id="my-pla",
        e1=3000.0,
        e2=2400.0,
        e3=1800.0,
        nu12=0.35,
        nu13=0.30,
        nu23=0.30,
        g12=1000.0,
        g13=700.0,
        g23=650.0,
        density=1.24e-9,
    ),
)
edit.set_material(project, "beam", "my-pla")
edit.set_print_setup(project, "beam", "fdm", build_direction=(0.0, 0.0, 1.0), raster_angle=0.0)
save_project(project, path)
```

(These numbers are made up for the example, not real PLA data.) In the JSON, the material goes
in a top-level `"materials"` list with `"kind": "orthotropic"`, and the part gets
`"print_setup": {"process": "fdm", "build_direction": [0, 0, 1], "raster_angle": 0}`. Use an id
that is not a library id: `add_material` refuses one (see `docs/user/project-files.md`).

The solve then prints a note: von Mises stress still shows where the stress is, but it is not a
strength check for a printed part.

## Step 4: Find the face numbers

Supports and loads sit on faces. Each face has a number. List them with:

```
tessera faces beam.tessera.json beam
```

The last argument is the part's `id`. For the example beam, the output is:

```
   1  Plane          100.000 mm2  (0.00, 5.00, 5.00)  (-1.00, +0.00, +0.00)
   2  Plane          100.000 mm2  (200.00, 5.00, 5.00)  (+1.00, -0.00, -0.00)
   3  Plane         2000.000 mm2  (100.00, 0.00, 5.00)  (+0.00, -1.00, -0.00)
   4  Plane         2000.000 mm2  (100.00, 10.00, 5.00)  (-0.00, +1.00, +0.00)
   5  Plane         2000.000 mm2  (100.00, 5.00, 0.00)  (-0.00, +0.00, -1.00)
   6  Plane         2000.000 mm2  (100.00, 5.00, 10.00)  (+0.00, -0.00, +1.00)
```

Columns: face number, surface type, area, centre point (mm), and outward normal.

- The **normal** points out of the part. Use it to find the face you want: face 1 is the end at
  x = 0 (normal −x), face 2 is the end at x = 200 (normal +x).
- A **curved** face (a hole, a fillet) has no single normal, so it is shown without one. STL
  faces are listed as `Plane` or `Curved`.
- The numbers come from the file. If you export the CAD model again, run `tessera faces` again
  and check the numbers before you trust an old study.

## Step 5: Add supports

A support holds chosen directions of a face at zero movement. The directions are numbered:
**1 = x, 2 = y, 3 = z**.

**Fixed support** (the default, all three directions held). Use it where the part is clamped:

```
tessera support add beam.tessera.json beam 1 --id root
```

Output:

```
added support 'root'
```

**Roller support** (only the directions you list are held). Use it on a plane of symmetry, or
where a face may slide along a surface. The face can still move in the other directions. Give
the directions with `--dofs`, using the letters x, y and z:

```
tessera support add beam.tessera.json beam 5 --dofs z --id roller
tessera info beam.tessera.json
```

Output of `info` (only the lines about supports are shown here):

```
  support root: beam faces 1, fixed
  support roller: beam faces 5, holds z
```

This holds face 5 from moving in z and leaves it free to move in x and y. **Do not keep this
roller.** Face 5 is the whole underside of the beam. Holding it in z stops the beam from
bending, so the answer becomes almost zero. The example removes it:

```
tessera remove beam.tessera.json roller
```

Output:

```
removed support 'roller'
```

To see the effect, a study with the roller kept gave a largest displacement of 0.000126 mm.
The same study without it gave 1.52177 mm (Step 9).

Rules:

- The face numbers come from Step 4. A face that does not exist is refused. You can list
  several faces, separated by spaces, for example `1 3`.
- `--dofs` takes some of the letters x, y and z. Each letter may appear once. `xz` and `z` are
  valid. `zz` is refused.
- A study needs at least one support. Without one, the solve stops with
  `the study needs at least one support`.
- A fixed face is stiffer than a real clamp that grips only part of the surface. Keep a fixed
  support to the area that is really fixed. Use a roller where the real part can slide.

### Line supports (on an edge)

A face support is a clamp: the face cannot turn. Where the real part rests on a line, such as a
shelf on a ledge, a beam on two rods or a plate on a ring, hold an **edge** instead. The part
can then turn about the line, as the real one does.

An edge is named by the two faces that meet there, written `A-B`. List a part's edges with:

```
tessera edges beam.tessera.json beam --face 5
```

`--face 5` shows only the edges of face 5 (the underside); leave it out to see them all. For
the example beam, the output is:

```
      1-5  Line             10.000 mm  at (0.00, 5.00, 0.00)  from (0.00, 0.00, 0.00) to (0.00, 10.00, 0.00)
      2-5  Line             10.000 mm  at (200.00, 5.00, 0.00)  from (200.00, 0.00, 0.00) to (200.00, 10.00, 0.00)
      3-5  Line            200.000 mm  at (100.00, 0.00, 0.00)  from (0.00, 0.00, 0.00) to (200.00, 0.00, 0.00)
      4-5  Line            200.000 mm  at (100.00, 10.00, 0.00)  from (0.00, 10.00, 0.00) to (200.00, 10.00, 0.00)
```

Columns: the two faces, the curve type, the length, a point on the edge, and its two ends. A
circle or another closed edge shows `closed` instead of the ends.

A beam simply supported at both ends rests on the bottom edge of each end face: a pinned line
at one end (all directions held) and a roller line at the other (only z held, so the beam can
get longer as it bends):

```
tessera support add beam.tessera.json beam --edge 1-5 --id pin
tessera support add beam.tessera.json beam --edge 2-5 --dofs z --id roller
```

`tessera info` then shows:

```
  support pin: beam edges 1-5, fixed
  support roller: beam edges 2-5, holds z
```

Rules:

- Give face numbers or `--edge`, not both. Repeat `--edge` to hold several edges in one
  support, for example `--edge 3-5 --edge 4-5` (both long bottom edges of the beam).
- The two faces must meet. `--edge 1-2` on the beam is refused, because its two end faces do
  not touch: `faces 1 and 2 of part 'beam' do not meet at an edge (list the edges with
  `tessera edges`)`.
- If two faces meet along more than one separate curve, the edge holds all of them.
- **A line in the middle of a face is not an edge.** To hold a line across a face, add a split
  line in your CAD program and export again. The split line divides the face in two, and the
  line becomes the edge between them.
- The seam of a cylinder wall (the joint line where the CAD kernel closes the surface) is not
  an edge and cannot be held.
- **Peak stress at a line support depends on the mesh.** A line that carries a force has, in
  theory, an infinite stress right at the line. A finer mesh shows a higher peak there, without
  limit. Displacements away from the line are reliable. For this beam under 0.1 MPa on its top
  face, a 2.5 mm mesh gave 0.120414 mm and a peak of 30.8 MPa, and a 1.25 mm mesh gave
  0.120482 mm and a peak of 62.4 MPa, at the supports. `tessera solve` marks every line
  support, so that peak is listed as near the support (see "Why the peak stress changes with
  the mesh" in Step 9). Read the stress away from the support line (here, in the middle of the
  span), or check the displacement.

## Step 6: Add loads

Each load has an `--id`, the part and face numbers it acts on, and a value. Give the value
with the option shown for each kind. Negative numbers work as option values: `--force 0 0 -100`
is valid.

### Force

```
tessera load force beam.tessera.json beam 2 --force 0 0 -100 --id tip-push
```

Output:

```
added load 'tip-push'
```

- `--force` is the **total** force in N, as three numbers: FX FY FZ.
- It is spread evenly over the area of the faces. Here 100 N over 100 mm² is 1 N/mm² (1 MPa).
- The sign gives the direction. `-100` in z pushes down, if z points up.

### Pressure

```
tessera load pressure beam.tessera.json beam 3 --pressure 0.5 --id press
```

Output:

```
added load 'press'
```

- `--pressure` is in MPa (N/mm²).
- It acts **normal to the face**, so you do not give a direction.
- **Positive pushes into the part**, the way water pressure pushes on the inside of a tank wall.
  Negative pulls out of the part.
- It follows curved faces. A 0.5 MPa load on a 100 mm² face is 50 N.

### Moment

```
tessera load moment beam.tessera.json beam 2 --moment 0 0 1000 --id twist
```

Output:

```
added load 'twist'
```

- `--moment` is a torque in N·mm, as three numbers: MX MY MZ. Here 1000 N·mm is 1 N·m.
- It turns the faces about their area centre. Use the right-hand rule: a positive z moment
  turns the face counter-clockwise when you look down from +z.
- The faces must cover an area. A single line or edge cannot carry a moment.
- The moment adds no net push or pull; it only twists.

The pressure and moment loads above are shown only to demonstrate the commands. The example
study keeps just the tip force, so remove the other two:

```
tessera remove beam.tessera.json press
tessera remove beam.tessera.json twist
```

Output:

```
removed load 'press'
removed load 'twist'
```

Check what is left with `tessera info` (Step 8).

### Units at a glance

| Quantity | Unit on the command line |
|---|---|
| Length, position | mm |
| Force | N |
| Pressure, stress, Young's modulus | MPa (= N/mm²) |
| Moment | N·mm |
| Acceleration (gravity) | mm/s² |
| Density (in the library) | t/mm³ |

Tessera does not read units from the command. Convert any value you took from a datasheet
before you type it in. For example, 1 kN is 1000 N, and 1 N·m is 1000 N·mm.

## Step 7: Set gravity and the mesh size

Use `tessera study` for both settings:

```
tessera study beam.tessera.json --gravity 0 0 -9810 --mesh-size 2.5
```

Output:

```
study updated
```

- `--gravity` takes three numbers, GX GY GZ, in mm/s². Standard gravity is 9810 mm/s² (9.81
  m/s²). If z points up, gravity is `0 0 -9810`. It acts on the whole part, using the density of
  the material. Both library materials have a density.
- `--no-gravity` turns gravity off.
- `--mesh-size` is the target size of each small element, in mm. The solver works on a mesh of
  these elements.
- `--auto-mesh` lets the mesher choose the size.
- Give either `--gravity` or `--no-gravity`, not both. The same goes for `--mesh-size` and
  `--auto-mesh`. A `study` command with none of these options is refused.

For the example beam, gravity adds about 1.54 N, the weight of the bar: 7.85e-9 t/mm³ × 20 000
mm³ = 1.57e-4 t, and 1.57e-4 t × 9810 mm/s² is about 1.54 N. The study in this guide leaves
gravity off, so turn it off again:

```
tessera study beam.tessera.json --no-gravity
```

Output:

```
study updated
```

Smaller mesh sizes give more elements, a slower solve and a more exact answer. The example
below uses 2.5 mm. If you leave the size out, the mesher picks one. Set it yourself so you
know what you get.

Always run the study at two sizes, for example 2.5 mm and 1.25 mm. If the numbers change
little, the mesh is fine enough. The peak stress at a support or a sharp inside corner changes
most; `solve` marks it (see Step 9). `tessera solve --check-mesh` runs both sizes for you and
says which numbers settled (see "Let Tessera compare two meshes" in Step 9).

## Step 8: Check the project and fix it

Summarise the project:

```
tessera info beam.tessera.json
```

Output:

```
beam: 1 parts, 0 connections
  part beam: steel-generic
  support root: beam faces 1, fixed
  load tip-push: force (0, 0, -100) N on beam faces 2
  mesh size: 2.5 mm
```

`info` shows the parts, supports, loads, gravity (when it is on) and the mesh size.

To delete a support, load or part, use its id:

```
tessera remove beam.tessera.json tip-push
```

A part can be removed only when no support, load or connection uses it. Otherwise the message
lists the items that still use it:

```
error: part 'beam' is still used by root, tip-push; remove those first
```

Remove those items first, then the part.

## Step 9: Solve and read the result

```
tessera solve beam.tessera.json --workdir run-2.5
```

Output:

```
nodes: 13213, elements: 7263
max displacement: 1.52177 mm
von Mises stress away from supports and sharp inside edges: 116.508 MPa at (5.684, 5.092, 0) mm, face 5
max von Mises stress: 120.398 MPa at (1.83, 8.17, 10) mm, face 6
  near support 'root' (within 5.39 mm): the stress there may keep rising as the mesh is refined; compare two mesh sizes (--check-mesh) before trusting either value
max von Mises stress by face (5 of 6):
  face 6: 120.398 MPa at (1.83, 8.17, 10) mm, near support 'root'
  face 5: 119.969 MPa at (1.83, 8.17, 0) mm, near support 'root'
  face 3: 117.341 MPa at (1.25, 0, 0) mm, near support 'root'
  face 4: 117.324 MPa at (1.25, 10, 10) mm, near support 'root'
  face 1: 114.319 MPa at (0, 0, 0) mm, near support 'root'
solver files: run-2.5
```

`--workdir` names the folder for the solver files. If you leave it out, Tessera makes a new
folder in the system temporary directory.

How to read it:

- **nodes, elements**: the size of the mesh. More is slower and more exact.
- **max displacement (mm)**: the largest movement of any point of the part. Beam theory gives
  about 1.53 mm for this bar (bending 1.524 mm plus shear 0.003 mm). The solver gave 1.52 mm.
  A linear model only holds while the movement is small compared with the part. If the
  largest movement is a large fraction of the part's size, the answer is not reliable.
- **note: this part uses an orthotropic (layer-wise) material**: printed only when the part's
  material is orthotropic. Von Mises stress still shows where the stress is, but it is not a
  strength check: a printed part is weaker across its layers than along them, so do not compare
  these numbers with a datasheet strength. Layer-wise failure checks are planned for v0.5.
- **von Mises stress away from supports and sharp inside edges (MPa)**: the highest
  equivalent stress outside the marked places (see below), with its location and face. It is
  not always the stress to compare with the **yield strength** from the material's datasheet:
  read "Why the peak stress changes with the mesh" first. Tessera v0.1 does not store yield
  strength. If every node is in a marked place, the line says `none`; use a finer mesh.
- **max von Mises stress (MPa)**: the largest equivalent stress at any node, with its location
  and face (`inside the part` if it is not on a face). When it is in a marked place, the next
  line says which, and that the value depends on the mesh size.
- **max von Mises stress by face**: the highest stress on each face, highest first, so you can
  see where the part is loaded. Five faces are shown; the Python API has them all.
- **solver files**: the folder with the solver input (`job.inp`) and results (`job.frd`).
  Keep it when you report a problem.

### Why the peak stress changes with the mesh

A model of a part has places where the stress does not settle to one value as the mesh is
refined. It keeps rising with every finer mesh. Two kinds of place are marked in the output:

- **near a support** that holds two or three directions, such as a fixed face. The clamped
  edge is the problem. A face support that holds one direction (a roller or symmetry plane) is
  not marked, because it does not cause this. Every line support is marked, whatever it holds,
  because the force it carries is squeezed onto a line.
- **near a sharp inside edge**: an edge where the surface turns inward (into the material) by
  30° or more, like the inside corner of an L-bracket. Outside edges, rounded fillets and holes
  are not marked.

Peaks under small, point-like loads and at the edges of a loaded face can also keep rising with
the mesh. These places are not marked yet, so a peak there shows up as an unmarked value.

"Near" means within twice the typical element size. The output line says how far that is (5.39
mm for the example at 2.5 mm). The value is shown, not hidden: a real stress raiser, such as a
hole next to a clamp, can sit in a marked area.

The example beam at two mesh sizes:

| `mesh_size` (mm) | Nodes | Max displacement (mm) | Away from supports and sharp inside edges (MPa) | Max von Mises (MPa) |
|---|---|---|---|---|
| 2.5 | 13 213 | 1.52177 | 116.508 | 120.398, near support 'root' |
| 1.25 | 78 355 | 1.52238 | 116.435 | 143.79, near support 'root' |

The largest displacement and the stress away from the support barely changed (0.04 % and
0.06 %). The peak at the clamped end rose by 19 %. Do not compare marked peaks from different
mesh sizes. A value that settles is not always the stress you need, though, as the next
paragraphs show.

**The stress away from supports and sharp inside edges can read low.** It is the highest
stress outside the marked places, and it leaves out the place where the stress is often
highest. Where the highest stress in the part is expected at a support, it can read low. The
common case is a part clamped at its base and bent, like the example beam. Beam theory gives
120 MPa at the clamped end of the example (bending moment 100 N × 200 mm over the section
modulus 167 mm³). The value away from the support is 116.5 MPa, about
3 % low, because it is taken about 5.4 mm from the clamp. On a short part with a coarse mesh
the gap can be a third of the stress or more.

In such a case the real stress lies between the two numbers: the stress away from supports and
sharp inside edges, and the marked peak. Do not take either one on its own as the answer:

- Run the two-mesh check and see which numbers settle: `tessera solve --check-mesh` (below),
  or two solves by hand as in Step 7.
- Or judge from the trend: if the stress away from the support sits next to the marked area
  and rises toward it, the stress at the support is higher than that value.

Where the highest stress is away from every marked place, for example at a hole in the middle
of a plate, the two numbers are the same and the value is not marked.

To run the 1.25 mm check, change the mesh size and solve again:

```
tessera study beam.tessera.json --mesh-size 1.25
tessera solve beam.tessera.json --workdir run-1.25
```

### Let Tessera compare two meshes

`--check-mesh` solves the study twice: once at your mesh size, then with every element about
1/√2 as long, including the small ones gmsh puts at small holes and short edges (2 to 3 times
as many nodes). It takes 3 to 5 times as long as one solve. Every number it prints comes from
the finer mesh, and each one says how much it changed between the two:

```
tessera solve beam.tessera.json --workdir run-2.5 --check-mesh
```

Output:

```
nodes: 31807, elements: 18911
mesh check: solved with median element edge 2.7 mm (13213 nodes), then 2.05 mm; the numbers below are from the finer mesh, settled = changed by at most 5 %
max displacement: 1.52221 mm (settled, changed 0.029 % between the meshes)
von Mises stress away from supports and sharp inside edges: 116.484 MPa at (5.251, 3.927, 0) mm, face 5 (settled, changed 0.021 % between the meshes)
max von Mises stress: 130.572 MPa at (0, 10, 10) mm, face 1
  near support 'root' (within 4.1 mm): mesh-dependent, changed 7.79 % between the meshes; do not compare it with a strength
max von Mises stress by face (5 of 6):
  face 1: 130.572 MPa at (0, 10, 10) mm, near support 'root' (mesh-dependent, changed 12.4 % between the meshes)
  face 4: 130.572 MPa at (0, 10, 10) mm, near support 'root' (mesh-dependent, changed 10.1 % between the meshes)
  face 6: 130.572 MPa at (0, 10, 10) mm, near support 'root' (mesh-dependent, changed 7.79 % between the meshes)
  face 5: 129.983 MPa at (0, 10, 0) mm, near support 'root' (mesh-dependent, changed 7.7 % between the meshes)
  face 3: 129.891 MPa at (0, 0, 0) mm, near support 'root' (mesh-dependent, changed 9.66 % between the meshes)
solver files: run-2.5
```

- **settled**: the number changed by 5 % or less, so a finer mesh would not change it much.
- **mesh-dependent**: it changed by more than 5 %. Do not compare it with a strength. Use a
  finer mesh, or read it as described in "The stress away from supports and sharp inside edges
  can read low" above.
- **not judged**: the finer mesh did not get finer enough where this number is (its elements
  there are more than 0.85 times their size on the first mesh), so the check says nothing about
  it. Any mark on it still applies.
- A measured verdict replaces the marks: a marked peak that settled can be used, and an
  unmarked one that did not settle cannot (but see the next paragraph).
- Each number is compared with the same kind of number on the other mesh: the largest stress
  with the largest stress, face 3's peak with face 3's peak. If two places have almost the same
  peak, the two meshes may put it in different places; check that the location is the same.
- The stress away from supports and sharp inside edges is taken a little closer to the marked
  places on the finer mesh, because the marked band is twice the element size. Part of its
  change can come from that.
- The element edge is measured on the mesh, so it can differ a little from `--mesh-size`.
- The coarse solve's files go to the `coarse` folder inside `--workdir`.

**A very coarse mesh can hide a problem.** The check sees how much a number moves, and on a
mesh too coarse to show the stress rise at a clamp, the peak there barely moves. The example
beam at 5 mm, two elements through its 10 mm depth, reads `settled, changed 1.12 %` at the
clamped end; at 2.5 mm, four elements, the same peak is `mesh-dependent`. If a thin wall
carries the load and has only one or two elements across it, refine before you trust a
`settled` there.

### If the solve takes too long

Tessera stops the solver if one solve takes longer than 30 minutes (1800 s), so a solve that
hangs cannot block you forever. It also stops the programs the solver started, unless one of
them detached itself from the solver. Closing the terminal, or stopping Tessera with Ctrl+C,
stops the solver too; if Tessera is killed outright (for example from the task manager), check
that no `ccx` process is left running. To allow more or less time, give the limit in seconds:

```
tessera solve beam.tessera.json --timeout 3600
```

With `--check-mesh` each of the two solves gets this long, so the whole command can take up to
twice as long, plus the meshing. The message says which of the two solves was stopped. A part that needs more than 30 minutes usually has a mesh
finer than it needs; try a larger `--mesh-size` first (Step 7).

### What to change in the print

Tessera can tell you where the stress is high, but it does not yet know your printed material.
These suggestions are general design advice; Tessera does not check them.

- **Stress peak at a sharp corner or a point load.** Add a fillet (a rounded corner) in the CAD
  model and export again. Or spread the load over a larger face.
- **Too much movement.** For a beam in bending, depth matters most: stiffness goes with the
  cube of the depth. A small increase in depth can cut the movement a lot.
- **Printed parts.** Layers are weaker when pulled apart than when pulled along. Orient the part
  so the main load runs along the layers, not across them. Tessera v0.1 cannot check this.

### Exit codes

| Code | Meaning |
|---|---|
| 0 | The command finished |
| 1 | Something in the project or the command is wrong. The message says what |
| 2 | The solver failed or took longer than `--timeout`. The message says why; check `ccx` is installed |

## Editing the JSON directly

The project file is plain JSON when its name ends in `.tessera.json`. The commands above write
the same data. You can also write the file by hand. Save this as `beam.tessera.json` in the same
folder as `beam.step`:

```json
{
  "schema_version": 5,
  "name": "Beam cantilever",
  "assembly": {
    "parts": [
      {
        "id": "beam",
        "name": "Beam",
        "geometry": "beam.step",
        "material_id": "steel-generic"
      }
    ]
  },
  "study": {
    "analysis": "static",
    "supports": [
      {
        "id": "root",
        "faces": {"part_id": "beam", "faces": [1]}
      }
    ],
    "loads": [
      {
        "kind": "force",
        "id": "tip-push",
        "faces": {"part_id": "beam", "faces": [2]},
        "force": [0.0, 0.0, -100.0]
      }
    ],
    "mesh_size": 2.5
  }
}
```

This study clamps the beam at x = 0 and pushes its free end down with 100 N. Run
`tessera faces beam.tessera.json beam` first if your own geometry has different face numbers.

The same fields take these forms. Add them to the lists shown:

- in `study.loads`: `{"kind": "pressure", "id": "press", "faces": {"part_id": "beam", "faces": [3]}, "pressure": 0.5}`
- in `study.loads`: `{"kind": "moment", "id": "twist", "faces": {"part_id": "beam", "faces": [2]}, "moment": [0.0, 0.0, 1000.0]}`
- in `study.supports`: `{"id": "side", "faces": {"part_id": "beam", "faces": [5]}, "dofs": [3]}`
- in `study.supports`, a line support: `{"id": "pin", "edges": {"part_id": "beam", "edges": [[1, 5]]}}`
- as a field of `study`: `"gravity": [0.0, 0.0, -9810.0]`

A load without a `kind` is read as a force load.

Any field name that Tessera does not know causes an error. This catches typing mistakes early.

## Messages you may see

| Message (after `error:`) | What to do |
|---|---|
| `part 'beam' has faces 1 to 6; there is no face 99 (list them with `tessera faces`)` | Use a face number from `tessera faces` |
| `unknown material 'steel-x'; known: aluminium-6061-t6, steel-generic` | Use a material id from `tessera materials` |
| `no part 'gear' in project 'beam'` | Check the part id with `tessera info`. The message names the project, which may differ from the file name |
| `the id 'root' is already used in this project` | Give the new item another `--id` |
| `nothing called 'x' in project 'beam'` | `tessera remove` found no support, load or part with that id |
| `part 'beam' is still used by root, tip-push; remove those first` | Remove the listed items first |
| `--dofs takes the directions to hold, e.g. xyz or z; got 'zz'` | Use each letter x, y or z once |
| `a face is listed twice: [2, 2]` | List each face once |
| `faces 1 and 2 of part 'beam' do not meet at an edge (list the edges with `tessera edges`)` | Use a pair from `tessera edges` |
| `give face numbers or --edge, not both` | A support holds faces or edges; add two supports for both |
| `--edge takes two face numbers like 3-5; got '1x5'` | Write the edge as two face numbers joined by `-` |
| `an edge is named by two different faces, e.g. 3-5; got 5-5` | Use two different faces |
| `edge 1-5 is listed twice` | List each edge once (`1-5` and `5-1` are the same edge) |
| `use either --gravity or --no-gravity` | Give one of the two |
| `nothing to change: give --gravity, --no-gravity, --mesh-size or --auto-mesh` | Give a setting to change |
| `the mesh size must be positive, got 0 mm` | Use a number greater than 0 |
| `invalid ForceLoad: ... Input should be a finite number` | A value is not a number, or is `nan` or infinite. Check the numbers |
| `geometry file not found: <path>` | `part add` cannot find the file. Check the path you typed |
| `geometry of part 'beam' not found: <path>` | The solve cannot find the geometry file the project names. Check it is still there |
| `the study needs at least one support` | Add a support (Step 5) |
| `v0.1 solves a single part; assemblies arrive in v0.5` | Keep one part in the project |
| `part 'beam' has no material` | Set the material (Step 3); it is needed to solve |
| `part 'beam' uses orthotropic material 'my-pla' and needs a print setup (build direction)` | Give the part a print setup (Step 3, "Your own printed material") |
| `'steel-generic' is a built-in material; choose another id` | Give your own material a different id |
| `... is not a valid project file` | Check the JSON: a misspelt field name, a missing comma or a missing required field |
| `[Errno 2] No such file or directory: ...` | The project file you named does not exist |
| `solver failed: CalculiX executable 'ccx' not found on PATH ...` | Install CalculiX 2.21 or newer, or set `TESSERA_CCX` |
| `solver failed: CalculiX did not finish within 1800 s and was stopped; give it longer with --timeout SECONDS` | Use a larger `--mesh-size`, or allow more time with `--timeout` (see "If the solve takes too long") |
| `the solver timeout must be a positive number of seconds, got 0` | Give `--timeout` a number of seconds greater than 0 |
