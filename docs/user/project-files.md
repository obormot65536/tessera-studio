# Project files

A Tessera project holds the parts you want to study, a material for each part, and the study
itself (supports, loads and mesh size). You can save it in two forms. This page explains both,
how to make and share them, and what happens when you open one.

Commands below use the `tessera` command-line tool. Run `tessera --help` for the full list.

## The two forms

| | `name.tessera` | `name.tessera.json` |
|---|---|---|
| What it is | A zip file: `project.json` plus a copy of each geometry file | Plain text (JSON) |
| Geometry | Copied inside the file | Referenced by path, relative to the JSON file |
| Edit by hand | Not directly; unzip it first | Yes, in any text editor |
| Share | One file | The JSON and its geometry files, kept in the same folder layout |
| Opens after the original STEP/STL moved | Yes | Only if the paths still point to the files |

Any file name that does not end in `.tessera` is saved as plain project JSON. For example,
`bracket.json` is JSON too.

**Use `.tessera.json` while you work.** It is easy to read, change and keep in version control.

**Use `.tessera` to share or archive.** It holds everything the project needs in one file.

Only STEP (`.step`, `.stp`) and STL (`.stl`) geometry can go into a `.tessera` file. A part
that uses another file type cannot be packed.

## Make a new project

```
tessera new bracket.tessera.json
tessera new bracket.tessera
```

The first command creates an editable JSON project. The second creates an empty container.
Both start with no parts and no study. Add parts and study settings by editing the JSON, or
through the Python API (see `docs/api/core.md`). The name of the project is taken from the file
name; use `--name "Bracket v2"` to set another name.

## Pack a project for sharing

```
tessera pack bracket.tessera.json bracket.tessera
```

`pack` reads the JSON project, reads every part's STEP or STL file, and writes one `.tessera`
file that holds both. The target must end in `.tessera`. If the target already exists, it is
replaced.

Pack stops with an error, and writes nothing, when:

- a part's geometry file does not exist (the message says `not found`);
- a part uses a file that is not STEP or STL.

Fix the path in the JSON and run `pack` again.

## Check a project

```
tessera info bracket.tessera
```

This prints the project name, the number of parts and the number of connections. It does not
check that the geometry files exist for a `.tessera.json` file; `solve` and `faces` do.

## What happens when you open a .tessera file

Tessera copies the geometry out of the container into a private temporary folder. The copy
lives until Tessera exits, and is then deleted. Your original file is not changed.

Because the geometry lives in that temporary folder, a project opened from a `.tessera` file
**cannot be saved as `.tessera.json`** in the same session. Tessera says so and does not save.
You have two choices:

- save it as `.tessera` (the copy is written again, and nothing is lost), or
- unzip the file and edit the JSON inside it (next section).

## Edit a project someone sent you

A `.tessera` file is a normal zip file. Unzip it with any zip tool, or on the command line:

```
unzip bracket.tessera -d bracket
```

The folder now holds:

```
bracket/
  project.json
  geometry/
    bar.step
```

`project.json` is a complete project. Its geometry paths (`geometry/bar.step`) are relative to
the folder, so you can open and edit it directly:

```
tessera info bracket/project.json
```

When you are done, either keep working with `project.json` as your editable form, or pack the
folder again:

```
tessera pack bracket/project.json bracket-edited.tessera
```

## What a project file records

Every file stores a `schema_version` number. The current number is **6** (5 added line
supports; 6 added materials defined in the project and the raster angle of a print setup). It
is the version of the file layout, not of Tessera. Tessera saves every file in the current
version.

- **Older files open.** A project saved by an earlier version of Tessera still loads. Saving it
  again writes the current version.
- **Newer files are refused.** If a file says it uses a later schema than your Tessera reads,
  Tessera stops with a message such as:

  ```
  error: bracket.tessera uses project schema 7; this version of Tessera reads up to 6
  ```

  Update Tessera to open the file. Do not edit the number to make it open; the file may use
  fields this version does not know.

### Materials defined in the project

From schema 6, a project can hold its own materials in a top-level `materials` list, for
example the stiffness you measured on your own printed test pieces. Each entry has a `kind`
and an `id` that is unique within the project. Stiffness values are in MPa and density in
t/mm³ (1.24 g/cm³ is `1.24e-9`); `density` is optional.

```json
"materials": [
  {"kind": "isotropic", "id": "my-petg", "youngs_modulus": 2100.0, "poisson_ratio": 0.38},
  {"kind": "orthotropic", "id": "my-pla",
   "e1": 3000.0, "e2": 2400.0, "e3": 1800.0,
   "nu12": 0.35, "nu13": 0.30, "nu23": 0.30,
   "g12": 1000.0, "g13": 700.0, "g23": 650.0, "density": 1.24e-9}
]
```

An orthotropic material is stiffer in some directions than others, as a printed part is. Its
axes follow the print: 1 along the printed lines (the raster), 2 across the lines within a
layer, 3 up through the layers (the build direction). `nu12`, `nu13` and `nu23` are the major
Poisson ratios: `nu12` is how much the part shrinks along axis 2 when stretched along axis 1.
Tessera refuses a file whose constants cannot describe a real material (a stiffness of zero or
less, or Poisson ratios so large that stretching would release energy), and says which
material is wrong.

The examples above are made-up values, not data for real PLA or PETG. A part names a project
material by its id in `material_id`, as it would a library one. If a project material and a
library material share an id, the project's is used. A part with an orthotropic material needs
a print setup (below) to be solved.

### Print setup

A part's `print_setup` records how it is printed: `process` (`"fdm"` or `"resin"`),
`build_direction` (the part's direction that points up on the printer, default `[0, 0, 1]`)
and, from schema 6, `raster_angle`: the angle of the printed lines in degrees, counter-clockwise
when looking down the build direction (default `0`). The solve uses `build_direction` and
`raster_angle` to turn an orthotropic material's axes; an isotropic part ignores them.

Saving replaces the old file in one step. If a save fails, the previous file is kept.

## Sharing and bug reports

A `.tessera` file contains a copy of your geometry, so treat it as you would the CAD file
itself.

It does **not** contain the folder path of your original files or your user name. Tessera
does not store where the geometry came from.

When you report a bug, attach:

1. the `.tessera` file (or the JSON and its geometry, zipped together);
2. the exact command you ran and the full output it printed;
3. the solver folder that `tessera solve` printed (see `docs/user/study-setup.md`), if the
   problem happened during a solve. It holds the solver input and output files.

## Paths in a .tessera.json file

In a `.tessera.json` file, a part's `geometry` path is relative to the folder that contains the
JSON file. For example, `"geometry": "parts/bar.step"` means the `parts` folder next to the JSON.
An absolute path also works, but it will not match on another computer. A `.tessera` file
always stores the geometry under `geometry/`.
