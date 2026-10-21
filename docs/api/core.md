# API reference: `tessera.core`

`tessera.core` is the public Python API. The command-line tool and the planned desktop app call
it, and so can your own scripts. Import it like this:

```python
from pathlib import Path

from tessera.core import (
    describe_faces,
    edit,
    is_container_name,
    load_project,
    new_project,
    save_project,
    solve,
)
```

Those seven names are the whole public surface (`tessera.core.__all__`). `edit` is a module
with the functions that change a project (see below). Everything else in the package is
internal.

Units are mm, N, MPa, t and s throughout (ADR-002). Lengths are mm, forces N, pressures and
stresses MPa, moments N·mm, accelerations mm/s², and density t/mm³.

## Functions

### `new_project`

```python
def new_project(name: str) -> Project
```

Returns an empty project with the given name: no parts, no connections, and a default static
study. It does not touch the file system. It raises nothing in normal use.

### `is_container_name`

```python
def is_container_name(path: Path) -> bool
```

`True` when the file name ends in `.tessera` (case-insensitive). It only looks at the name. It
decides whether `save_project` writes a container or plain JSON.

### `save_project`

```python
def save_project(project: Project, path: Path) -> None
```

Writes the project to `path`.

- If `path` ends in `.tessera`, a zip container is written. It holds `project.json` and a copy
  of each part's geometry under `geometry/`.
- Any other name writes plain project JSON. Geometry paths are stored relative to the JSON
  file's folder.

The file is written to a temporary file in the same folder and then renamed over the target. A
failed save leaves the previous file in place. The saved file always carries the current
`schema_version` (4).

Raises:

- `FileNotFoundError`: a part's geometry file does not exist. The message names the part.
- `InputError`:
  - a container is requested and a part's geometry is not STEP (`.step`, `.stp`) or STL
    (`.stl`);
  - JSON is requested and a part uses geometry from a container that was opened in this
    process (that geometry lives in a temporary folder).
  - the in-memory project is invalid, for example a support or load names a part id that is
    not in `assembly.parts`. Nothing is written.
- `OSError`: the file cannot be written, for example because its folder does not exist.

### `load_project`

```python
def load_project(path: Path) -> Project
```

Reads a `.tessera` container or a project JSON file. The file type is found from its content,
not its name: a zip file is read as a container, anything else as JSON.

- JSON: part geometry paths are resolved against the JSON file's folder. The geometry files are
  not checked for existence here.
- Container: the geometry is copied to a private temporary folder, which Tessera removes when
  the process exits. The part geometry paths point there.

Files with a `schema_version` newer than `SCHEMA_VERSION` (4) are refused. Files without a
version, or with an older one, load.

Raises:

- `FileNotFoundError`: `path` does not exist.
- `InputError`:
  - `path` is not a file (for example, a folder);
  - the JSON cannot be parsed, or does not match the project model (unknown fields are
    errors);
  - the file uses a newer schema than this version reads (the message gives both numbers);
  - a container is damaged, encrypted or uses an unsupported compression method;
  - a container has no `project.json`, lists a geometry entry it does not contain, or lists a
    geometry name outside `geometry/`;
  - a geometry entry, or the total geometry, is larger than the size limit (1 GiB per entry,
    4 GiB in total).

### `describe_faces`

```python
def describe_faces(project: Project, part_id: str) -> list[FaceInfo]
```

Returns the faces of one part, ordered by face index. Use the `index` values in `FaceSet.faces`
to place supports and loads. Each `FaceInfo` gives the face's surface type, area, centroid and,
for planar faces, the outward unit normal.

For STEP the indices are the CAD kernel's face numbers. For STL they are numbered by Tessera
from the triangles (ADR-005), so they change if the STL file changes.

Faces are remembered by the file's content for the rest of the process, so asking again for an
unchanged file (as `tessera.core.edit` does for every support and load) does not read it again;
a file that changed on disk is read again.

Raises:

- `InputError`:
  - the project is invalid (for example two parts share an id, or a support names an unknown
    part);
  - no part has the id `part_id`;
  - the part's geometry is not STEP or STL;
  - the file cannot be read, or a STEP file does not contain exactly one solid;
  - an STL file is not a single closed solid, or has inconsistent winding or disconnected
    surfaces.
- `FileNotFoundError`: the part's geometry file does not exist. The message names the part.

### `solve`

```python
def solve(
    project: Project, workdir: Path | None = None, adapter: SolverAdapter | None = None
) -> StudyResult
```

Meshes the part, builds the static problem, runs the solver and summarises the result.

- `workdir`: the folder for solver files (`job.inp`, `job.frd` and CalculiX logs). The folder
  is created if needed. If `None`, a new folder named `tessera-*` is created in the system
  temporary directory. Tessera does not delete it.
- `adapter`: the solver. `None` uses CalculiX (`tessera.solvers.CalculixAdapter`). The adapter
  must support the `static` analysis.

v0.1 solves exactly one part with no connections and no placement (identity transform). The
study must have at least one support. The part must have a `material_id` from the built-in
library.

Raises:

- `NotImplementedError`: the project has more than one part, any connection, or a placement
  other than the identity transform. These are planned for later versions.
- `InputError`:
  - the project is invalid, for example a support or load names a part that does not exist;
  - a face number in a support or load that the part does not have;
  - the part has no material, or the material id is unknown;
  - the study has no supports;
  - gravity is set but the material has no density;
  - the geometry is unsupported, or cannot be meshed;
  - the solver adapter cannot run the analysis.
- `FileNotFoundError`: the part's geometry file does not exist. The message names the part.
- `SolverError` (from `tessera.model`): the CalculiX executable cannot be found or started, the
  run timed out, CalculiX reported an error, or it produced no usable results. `SolverError.log`
  holds the solver output.

Returns a `StudyResult`.

### `tessera.core.edit`

```python
from tessera.core import edit
```

Functions that change a project in memory: parts, materials, supports, loads and study
settings. They do not read or write files, except that `add_part` reads the geometry file.
Call `save_project` afterwards to keep the change. The command line `tessera part`, `support`,
`load`, `study` and `remove` commands are thin wrappers over these functions.

`Vec3` in the signatures below means `tuple[float, float, float]`.

Common rules:

- **Ids.** Parts, connections, supports and loads share one set of ids, so an id is used once in
  the whole project. An empty id, or an id already in use, is refused. If an id argument is
  `None`, a new one is made: `support-1`, `support-2`, and so on for supports, `load-1` and so
  on for loads.
- **Faces.** Face numbers are checked against the part's real geometry (see `describe_faces`).
  A face number the part does not have is refused, and so is an empty list or a face listed
  twice.
- **Numbers.** Forces, pressures, moments and gravity must be finite (no `nan` or infinity).
  A zero force is accepted.
- **Refused changes.** A refused change leaves the project unchanged. The functions that build
  a new item or change a setting apply the change to a copy first, check the whole project,
  and only then change `project`. The other functions run all their checks before they touch
  `project`.
- **Geometry.** The functions that look at a part's faces read its geometry file. If the file
  is missing, `FileNotFoundError` is raised, and the message names the part.

#### `material_ids`

```python
def material_ids() -> list[str]
```

The ids of the built-in materials, sorted, for `material_id`. Currently `aluminium-6061-t6` and
`steel-generic`. It raises nothing.

#### `add_part`

```python
def add_part(
    project: Project,
    geometry_path: Path,
    *,
    part_id: str | None = None,
    name: str | None = None,
    material_id: str | None = None,
) -> PartInstance
```

Adds a STEP or STL file as a new part and returns it.

- `part_id` defaults to the file name, with characters other than letters, digits, `_` and `-`
  replaced by `-`. If that id is taken, `-2`, `-3`, ... is added until it is free.
- `name` defaults to the id. `material_id` is optional; if given it must be a built-in id.
- The file is read once, so a file Tessera cannot use is refused at once.
- The geometry path is stored as an absolute path in the project. `save_project` writes it
  relative to the JSON file.

Raises:

- `FileNotFoundError`: the geometry file does not exist.
- `InputError`: the file is not a STEP or STL file (by extension); the material is unknown; the
  id is empty or already used; the geometry cannot be read (for example, it is not a single
  closed STL solid); gravity is on and the material has no density; or the project is invalid
  after the change.

#### `set_material`

```python
def set_material(project: Project, part_id: str, material_id: str) -> None
```

Sets the material of a part.

Raises:

- `InputError`: no part has `part_id`; the material is unknown; gravity is on and the material
  has no density.
- `FileNotFoundError`: the part's geometry file does not exist.

#### `add_support`

```python
def add_support(
    project: Project,
    part_id: str,
    faces: Iterable[int],
    *,
    dofs: tuple[int, ...] = (1, 2, 3),
    support_id: str | None = None,
) -> Support
```

Holds the listed faces. `dofs` lists the held directions: 1 = x, 2 = y, 3 = z. The default holds
all three (a fixed support). The directions are stored in the order given; the command line
sorts them.

Raises `InputError` when:

- the id is empty or used, or `faces` is empty, has a face the part does not have, or lists a
  face twice;
- `dofs` is empty, has a value other than 1, 2 or 3, or lists a direction twice.

#### `add_force_load`

```python
def add_force_load(
    project: Project,
    part_id: str,
    faces: Iterable[int],
    force: Vec3,
    *,
    load_id: str | None = None,
) -> ForceLoad
```

A total force in N (`FX, FY, FZ`), spread evenly over the faces' area.

Raises `InputError` for an empty or used id, a bad face list (see the common rules), a
non-finite component, or a force of zero.

#### `add_pressure_load`

```python
def add_pressure_load(
    project: Project,
    part_id: str,
    faces: Iterable[int],
    pressure: float,
    *,
    load_id: str | None = None,
) -> PressureLoad
```

A pressure in MPa, normal to the faces. Positive pushes into the part.

Raises `InputError` for an empty or used id, a bad face list, or a pressure that is zero or
not finite.

#### `add_moment_load`

```python
def add_moment_load(
    project: Project,
    part_id: str,
    faces: Iterable[int],
    moment: Vec3,
    *,
    load_id: str | None = None,
) -> MomentLoad
```

A moment in N·mm (`MX, MY, MZ`) about the faces' area centroid.

Raises `InputError` for an empty or used id, a bad face list, a non-finite component, or a
moment of zero.

#### `set_gravity`

```python
def set_gravity(project: Project, acceleration: Vec3 | None) -> None
```

Gravity in mm/s² for every part, for example `(0, 0, -9810)`. `None` turns gravity off.

Raises `InputError` for a non-finite component, an invalid project, or a part whose material
has no density (gravity needs it).

#### `set_mesh_size`

```python
def set_mesh_size(project: Project, size: float | None) -> None
```

Target element size in mm. `None` lets the mesher choose.

Raises `InputError` when `size` is not greater than 0. This includes `nan`.

#### `remove`

```python
def remove(project: Project, item_id: str) -> str
```

Removes the support, load or part with this id, and returns which kind it was: `"support"`,
`"load"` or `"part"`. Supports and loads are looked at first.

A part can be removed only if no support, load or connection uses it. Otherwise the
`InputError` lists the items that still use it, and nothing is removed.

Raises `InputError` when nothing has this id, or when the part is still used.

## Model types

These come from `tessera.model` and are also available from `tessera.core` where noted. Import
them from `tessera.model` for the full set.

### Project and study

| Type | Fields (type, default) | Notes |
|---|---|---|
| `Project` | `name: str`, `schema_version: int = 4`, `assembly: Assembly`, `study: Study` | Root object. Extra fields are errors. |
| `Assembly` | `parts: list[PartInstance]`, `connections: list[Connection]` | Part ids must be unique. |
| `PartInstance` | `id: str`, `name: str`, `geometry: Path`, `transform: Transform`, `material_id: str \| None`, `print_setup: PrintSetup \| None` | `geometry` is a STEP or STL path. |
| `Transform` | `rotation`, `translation` | Rotation is a row-major 3 × 3 matrix; translation in mm. Only the identity is solved in v0.1. |
| `PrintSetup` | `process: "fdm" \| "resin"`, `slicer_profile: str \| None`, `build_direction` | Stored only. Not used by the solver in v0.1. |
| `Connection` | `id`, `kind: ConnectionKind`, `a: FaceSet`, `b: FaceSet`, `friction`, `preload_n` | Stored only. Connections are not solved in v0.1. |
| `ConnectionKind` | `"bonded"`, `"contact"`, `"bolt"`, `"pin"` | |
| `FaceSet` | `part_id: str`, `faces: tuple[int, ...]` | At least one face. Face numbers come from `describe_faces`. |
| `Study` | `analysis: "static"`, `supports: list[Support]`, `loads: list[Load]`, `gravity: (float, float, float) \| None`, `mesh_size: float \| None` | `gravity` in mm/s². `mesh_size` in mm, greater than 0; `None` lets the mesher choose. |
| `Support` | `id: str`, `faces: FaceSet`, `dofs: tuple[1 \| 2 \| 3, ...] = (1, 2, 3)` | Holds the listed displacement components at zero: 1 = x, 2 = y, 3 = z. No repeats. |
| `ForceLoad` | `kind: "force"`, `id: str`, `faces: FaceSet`, `force: (float, float, float)` | Total force in N, spread evenly over the faces' area. |
| `PressureLoad` | `kind: "pressure"`, `id: str`, `faces: FaceSet`, `pressure: float` | MPa, normal to the faces. Positive pushes into the part. |
| `MomentLoad` | `kind: "moment"`, `id: str`, `faces: FaceSet`, `moment: (float, float, float)` | N·mm about the faces' area centroid. No net force. |

`Load` is the union `ForceLoad | PressureLoad | MomentLoad`. Files without a `kind` are read as
force loads.

Validation: a `Study` item that names a part id not in `assembly.parts`, or two parts with the
same id, is an error. Items you add later with `append` skip pydantic's checks at that moment,
so `solve`, `describe_faces` and `save_project` check the whole project again and raise
`InputError`.

### Results and geometry

| Type | Fields | Notes |
|---|---|---|
| `StudyResult` | `mesh: Mesh`, `solution: SolveResult`, `max_displacement: float` (mm), `max_von_mises: float` (MPa), `workdir: Path \| None` | Returned by `solve`. `displacement_magnitude` gives the per-node displacement length. |
| `SolveResult` | `displacements: (N, 3)` mm, `stresses: (N, 6)` MPa in the order xx, yy, zz, xy, yz, zx, `solver: str`, `log: str`, `metadata: dict[str, str]` | Per node, 0-based. Nodes that no element uses are NaN. |
| `Mesh` | `nodes: (N, 3)` mm, `elements: (M, 10)` quadratic tetrahedra (C3D10 order), `element_sets`, `node_sets`, `surfaces` | `n_nodes` and `n_elements` properties. |
| `FaceInfo` | `index: int`, `kind: str`, `area: float` (mm²), `centroid: (x, y, z)` (mm), `normal: (x, y, z) \| None` | `kind` is the CAD surface type (`"Plane"`, `"Cylinder"`, ...) for STEP, and `"Plane"` or `"Curved"` for STL. `normal` is `None` for curved faces. |

`max_von_mises` is the largest nodal von Mises stress. Peak values at corners and load points
depend on the mesh size; see `docs/user/study-setup.md`.

### Errors

| Type | Module | Meaning |
|---|---|---|
| `InputError` | `tessera.model` | A problem in what the user asked for (a subclass of `ValueError`). Its message is safe to show as-is. |
| `SolverError` | `tessera.model` | The solver failed. Has a `log` attribute. Subclass of `RuntimeError`. |

Clients show `InputError` and `SolverError` messages to the user. `FileNotFoundError` and
`NotImplementedError` are also raised by `tessera.core` for the cases listed above. Any other
exception is a bug.

### Solver interface

`solve` accepts any object with this shape (`tessera.model.SolverAdapter`, ADR-004):

```python
class SolverAdapter(Protocol):
    name: str

    def capabilities(self) -> frozenset[AnalysisKind]: ...
    def run(self, job: SolveJob, workdir: Path) -> SolveResult: ...
```

`AnalysisKind` has one member, `STATIC = "static"`.

## Example

The example starts from a project file with one part (`bar`) and no supports or loads.

```python
from pathlib import Path

from tessera.core import describe_faces, load_project, save_project, solve
from tessera.model import FaceSet, ForceLoad, Support

project = load_project(Path("bar.tessera.json"))

faces = describe_faces(project, "bar")
root = next(f.index for f in faces if f.normal == (-1.0, 0.0, 0.0))
tip = next(f.index for f in faces if f.normal == (1.0, 0.0, 0.0))

project.study.supports.append(Support(id="root", faces=FaceSet(part_id="bar", faces=(root,))))
project.study.loads.append(
    ForceLoad(id="tip", faces=FaceSet(part_id="bar", faces=(tip,)), force=(0.0, 0.0, -100.0))
)
project.study.mesh_size = 2.5

result = solve(project, workdir=Path("run-coarse"))
print(result.max_displacement, result.max_von_mises)

save_project(project, Path("bar-with-study.tessera"))
```

`normal == (-1.0, 0.0, 0.0)` compares floats exactly. For real use, compare with a tolerance,
as the benchmarks do.
