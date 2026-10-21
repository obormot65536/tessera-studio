# Set up a study for one part

This guide walks through a complete study of one 3D-printed or machined part, from the file to
the numbers. It is written for makers. Every step uses units of **mm, N, MPa, and mm/s²**.

## What version 0.1 can and cannot do

Read this first. It saves you from setting up a study that Tessera cannot run.

- **One part per project.** No assemblies. No connections between parts. A part cannot be
  moved or rotated (its placement must stay at the identity).
- **Linear static analysis only.** The model assumes small displacements and a material that
  springs back when unloaded. It does not model buckling, large bending, plastic yield,
  contact, heat, vibration or time.
- **Isotropic materials only, from a short library.** Two materials exist: `steel-generic` and
  `aluminium-6061-t6`. There is no printed plastic (PLA, PETG, ABS, nylon) yet. So a result
  describes a part with the properties of metal. It does not yet show how a print behaves
  across its layers.
- **Two numbers come out.** The largest displacement and the largest von Mises stress. Tessera
  does not compare them with a limit, does not compute a safety factor and does not say
  whether the part will hold.
- **Needs CalculiX.** The solver `ccx` must be installed and on the PATH, or its location set
  in the `TESSERA_CCX` environment variable.

## Step 1: Get the part as a file

Export the part from your CAD tool as:

- **STEP** (`.step` or `.stp`): one solid. Export in millimetres. Tessera does not rescale
  the geometry.
- **STL** (`.stl`): one closed body. Every edge must be shared by two triangles; open edges
  are refused. Tessera splits the surface into faces wherever two neighbouring triangles meet
  at an angle of more than 40°. Lengths are read as mm.

Keep the geometry file in the same folder as the project file. The example below uses
`bar.step`, a 200 × 10 × 10 mm steel-like bar. You can make the same bar in any CAD tool to
follow along.

## Step 2: Create the project file

Start from the example at the end of this page, or create an empty JSON project:

```
tessera new bar.tessera.json
```

Open the file in a text editor and fill in the parts and the study. A project needs:

| Field | Meaning |
|---|---|
| `name` | Name of the project (text) |
| `assembly.parts` | The list of parts. v0.1 accepts exactly one |
| `id` | Short name for the part. You use it in face references |
| `name` (in a part) | A label for people |
| `geometry` | Path to the STEP or STL file |
| `material_id` | One of the library ids in Step 3. Needed for a solve |
| `study.supports` | At least one support (Step 5) |
| `study.loads` | Zero or more loads (Step 6) |
| `study.gravity` | Optional. Gravity vector (Step 6) |
| `study.mesh_size` | Element size in mm (Step 7) |

Any field name that Tessera does not know causes an error. This catches typing mistakes early.

## Step 3: Choose the material

Use the `material_id` exactly as written:

| `material_id` | Young's modulus E (MPa) | Poisson's ratio ν | Density (t/mm³) |
|---|---|---|---|
| `steel-generic` | 210 000 | 0.30 | 7.85e-9 |
| `aluminium-6061-t6` | 68 900 | 0.33 | 2.70e-9 |

Density is needed only if you use gravity. An unknown id is refused, and the message lists the
known ids.

If your part is a printed plastic, these values are not right for it. Use the closest library
entry only to get a rough idea, and do not treat the numbers as a print result.

## Step 4: Find the face numbers

Supports and loads sit on faces. Each face has a number. List them with:

```
tessera faces bar.tessera.json bar
```

The last argument is the part's `id`. For the example bar, the output is:

```
   1  Plane         100.000 mm2  (0.00, 5.00, 5.00)  (-1.00, +0.00, +0.00)
   2  Plane         100.000 mm2  (200.00, 5.00, 5.00)  (+1.00, -0.00, -0.00)
   3  Plane        2000.000 mm2  (100.00, 0.00, 5.00)  (+0.00, -1.00, -0.00)
   4  Plane        2000.000 mm2  (100.00, 10.00, 5.00)  (-0.00, +1.00, +0.00)
   5  Plane        2000.000 mm2  (100.00, 5.00, 0.00)  (-0.00, +0.00, -1.00)
   6  Plane        2000.000 mm2  (100.00, 5.00, 10.00)  (+0.00, -0.00, +1.00)
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

```json
{"id": "root", "faces": {"part_id": "bar", "faces": [1]}}
```

**Roller support** (only the listed directions held). Use it on a plane of symmetry, or where a
face may slide along a surface. The face can still move in the other directions:

```json
{"id": "side", "faces": {"part_id": "bar", "faces": [5]}, "dofs": [3]}
```

This holds the face at z = 0 from moving in z, and leaves it free to move in x and y.

Rules:

- The `faces` list holds face numbers from Step 4. A face that does not exist is refused.
- `dofs` holds some of 1, 2, 3. Do not list a number twice.
- A study needs at least one support. Without one, Tessera stops with
  `the study needs at least one support`.
- A fixed face is stiffer than a real clamp that grips only part of the surface. Keep a fixed
  support to the area that is really fixed. Use a roller where the real part can slide.

## Step 6: Add loads and gravity

Each entry in `study.loads` has an `id`, the `faces` it acts on, and a value. Use the `kind`
field to say which type it is.

### Force (`"kind": "force"`, or no kind)

```json
{"kind": "force", "id": "push", "faces": {"part_id": "bar", "faces": [2]}, "force": [0.0, 0.0, -100.0]}
```

- `force` is the **total** force in N, as `[Fx, Fy, Fz]`.
- It is spread evenly over the area of the faces. Here 100 N over 100 mm² is 1 N/mm² (1 MPa).
- The sign gives the direction. `-100` in z pushes down.

### Pressure (`"kind": "pressure"`)

```json
{"kind": "pressure", "id": "press", "faces": {"part_id": "bar", "faces": [3]}, "pressure": 0.5}
```

- `pressure` is in MPa (N/mm²).
- It acts **normal to the face**, so you do not give a direction.
- **Positive pushes into the part**, the way water pressure pushes on the inside of a tank wall.
  Negative pulls out of the part.
- It follows curved faces. A 0.5 MPa load on a 100 mm² face is 50 N.

### Moment (`"kind": "moment"`)

```json
{"kind": "moment", "id": "twist", "faces": {"part_id": "bar", "faces": [2]}, "moment": [0.0, 0.0, 1000.0]}
```

- `moment` is a torque in N·mm, as `[Mx, My, Mz]`. Here 1000 N·mm is 1 N·m.
- It turns the faces about their area centre. Use the right-hand rule: a positive z moment
  turns the face counter-clockwise when you look down from +z.
- The faces must cover an area. A single line or edge cannot carry a moment.
- The moment adds no net push or pull; it only twists.

### Gravity (top-level `study.gravity`)

```json
"gravity": [0.0, 0.0, -9810.0]
```

- Acceleration in **mm/s²**. Standard gravity is 9810 mm/s² (9.81 m/s²). If z points up, as
  in the examples here, gravity is `[0, 0, -9810]`.
- It acts on the whole part, using the density of the material. Both library materials have
  a density.
- For the example bar: 7.85e-9 t/mm³ × 20 000 mm³ = 1.57e-4 t. Times 9810 mm/s² that is about
  1.54 N, the weight of the bar.

### Units at a glance

| Quantity | Unit in the file |
|---|---|
| Length, position | mm |
| Force | N |
| Pressure, stress, Young's modulus | MPa (= N/mm²) |
| Moment | N·mm |
| Acceleration (gravity) | mm/s² |
| Density (in the library) | t/mm³ |

Tessera does not read units from the file. Convert any value you took from a datasheet
before you type it in. For example, 1 kN is 1000 N, and 1 N·m is 1000 N·mm.

## Step 7: Choose the mesh size

`mesh_size` is the target size of each small element, in mm. The solver works on a mesh of
these elements.

- Smaller numbers give more elements, a slower solve and a more exact answer.
- For the example bar, `mesh_size: 2.5` gave 13 213 nodes. The whole `tessera solve` command
  took about 7 seconds on the machine used to write this guide.
- If you leave `mesh_size` out, the mesher picks a size. Set it yourself so you know what you
  get.

Always run the study at two sizes, for example 2.5 mm and 1.25 mm. If the numbers change
little, the mesh is fine enough. The peak stress is the number that changes most (see Step 8).

## Step 8: Solve and read the result

```
tessera solve bar.tessera.json
```

For the example bar at `mesh_size: 2.5` the output is:

```
nodes: 13213, elements: 7263
max displacement: 1.52177 mm
max von Mises stress: 120.398 MPa
solver files: /home/you/.../tessera-1234abcd
```

How to read it:

- **nodes, elements**: the size of the mesh. More is slower and more exact.
- **max displacement (mm)**: the largest movement of any point of the part. Beam theory gives
  about 1.53 mm for this bar (bending 1.524 mm plus shear 0.003 mm). The solver gave 1.52 mm.
  A linear model only holds while the movement is small compared with the part. If the
  largest movement is a large fraction of the part's size, the answer is not reliable.
- **max von Mises stress (MPa)**: the largest equivalent stress at any node. It is a single
  number that combines the stress components. Compare it with the **yield strength** of the
  material from its datasheet. Tessera v0.1 does not store yield strength.
- **solver files**: the folder with the solver input (`job.inp`) and results (`job.frd`).
  Keep it when you report a problem.

### Why the peak stress changes with the mesh

The peak stress is the least reliable number here. In the example it sits at the clamped end,
on the edge and corner of the fixed face. At such a place the stress does not settle as the
mesh is refined. For the example bar:

| `mesh_size` (mm) | Max von Mises (MPa) | Peak position (mm) |
|---|---|---|
| 2.5 | 120.4 | (1.8, 8.2, 10.0), on the top edge of the clamped end |
| 1.25 | 143.8 | (0.0, 0.0, 10.0), a corner of the clamped end |

The largest displacement barely changed (1.5218 mm and 1.5224 mm). So the displacement is the
stable number, and the peak stress is not. Do not compare peak stresses from different mesh
sizes.

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
| 0 | The solve finished and printed the summary |
| 1 | Something in the project is wrong. The message says what |
| 2 | The solver failed. The message says why; check `ccx` is installed |

## Complete example

This file is a complete, valid v0.1 project. Save it as `bar.tessera.json` in the same folder as
`bar.step`. Run `tessera faces bar.tessera.json bar` first to check the face numbers for your
own geometry.

```json
{
  "schema_version": 4,
  "name": "Bar cantilever",
  "assembly": {
    "parts": [
      {
        "id": "bar",
        "name": "Bar",
        "geometry": "bar.step",
        "material_id": "steel-generic"
      }
    ]
  },
  "study": {
    "analysis": "static",
    "supports": [
      {
        "id": "root",
        "faces": {"part_id": "bar", "faces": [1]}
      }
    ],
    "loads": [
      {
        "kind": "force",
        "id": "tip-push",
        "faces": {"part_id": "bar", "faces": [2]},
        "force": [0.0, 0.0, -100.0]
      }
    ],
    "mesh_size": 2.5
  }
}
```

This study clamps the bar at x = 0 and pushes its free end down with 100 N.

The other forms from Steps 5 and 6, for reference. Add them to the lists shown:

- in `study.loads`: `{"kind": "pressure", "id": "press", "faces": {"part_id": "bar", "faces": [3]}, "pressure": 0.5}`
- in `study.loads`: `{"kind": "moment", "id": "twist", "faces": {"part_id": "bar", "faces": [2]}, "moment": [0.0, 0.0, 1000.0]}`
- in `study.supports`: `{"id": "side", "faces": {"part_id": "bar", "faces": [5]}, "dofs": [3]}`
- as a field of `study`: `"gravity": [0.0, 0.0, -9810.0]`

## Messages you may see

| Message (after `error:`) | What to do |
|---|---|
| `part 'bar' has no face 99` | Use a face number from `tessera faces` |
| `unknown material 'steel-x'; known: aluminium-6061-t6, steel-generic` | Fix `material_id` |
| `part 'bar' has no material` | Add `material_id` (needed to solve) |
| `the study needs at least one support` | Add a support (Step 5) |
| `v0.1 solves a single part; assemblies arrive in v0.5` | Keep one part in `assembly.parts` |
| `... is not a valid project file` | Check the JSON: a misspelt field name, a missing comma or a missing required field |
| `error:` followed by a path that does not exist | The `geometry` path is wrong. Paths are relative to the JSON file |
| `[Errno 2] No such file or directory: ...` | The project file you named does not exist |
| `solver failed: CalculiX executable 'ccx' not found on PATH ...` | Install CalculiX 2.21 or newer, or set `TESSERA_CCX` |
