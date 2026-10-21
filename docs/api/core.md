# API reference: `tessera.core`

`tessera.core` is the public Python API. The command-line tool and the planned desktop app call
it, and so can your own scripts. Import it like this:

```python
from pathlib import Path

from tessera.core import (
    describe_faces,
    is_container_name,
    load_project,
    new_project,
    save_project,
    solve,
)
```

Those six names are the whole public surface (`tessera.core.__all__`). Everything else in the
package is internal.

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
- `pydantic.ValidationError` (a subclass of `ValueError`): the in-memory project is invalid, for
  example a support or load names a part id that is not in `assembly.parts`. This is not
  wrapped as `InputError`. Nothing is written.
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

Raises:

- `InputError`:
  - no part has the id `part_id`;
  - the part's geometry is not STEP or STL;
  - the file cannot be read, or a STEP file does not contain exactly one solid;
  - an STL file is not a single closed solid, or has inconsistent winding or disconnected
    surfaces.
- `FileNotFoundError`: the part's geometry file does not exist.

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
  - a face number in a support or load that the part does not have;
  - the part has no material, or the material id is unknown;
  - the study has no supports;
  - gravity is set but the material has no density;
  - the geometry is unsupported, or cannot be meshed;
  - the solver adapter cannot run the analysis.
- `FileNotFoundError`: the part's geometry file does not exist.
- `SolverError` (from `tessera.model`): the CalculiX executable cannot be found or started, the
  run timed out, CalculiX reported an error, or it produced no usable results. `SolverError.log`
  holds the solver output.

Returns a `StudyResult`.

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

Validation: when a project is built or loaded, a `Study` item that names a part id not in
`assembly.parts` is an error. Items you add later with `append` are not checked against the
part list, and `solve` does not check `part_id` in a support or load (v0.1 solves one part, so
it uses that part's faces). Use the ids from `assembly.parts`.

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
