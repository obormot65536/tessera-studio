# API reference: `tessera.core`

`tessera.core` is the public Python API. The command-line tool and the planned desktop app call
it, and so can your own scripts. Import it like this:

```python
from pathlib import Path

from tessera.core import (
    MESH_CHECK_REFINEMENT,
    SOLVE_TIMEOUT,
    WALL_ELEMENTS_MIN,
    describe_edges,
    describe_faces,
    edit,
    is_container_name,
    load_project,
    new_project,
    save_project,
    solve,
)
```

Those ten names are the whole public surface (`tessera.core.__all__`). `edit` is a module
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
`schema_version` (5).

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

Files with a `schema_version` newer than `SCHEMA_VERSION` (5) are refused. Files without a
version, or with an older one, load. Files up to schema 4 hold only face supports and load
unchanged.

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

### `describe_edges`

```python
def describe_edges(project: Project, part_id: str) -> list[EdgeInfo]
```

Returns the edges of one part, for line supports (ADR-006). An edge is named by the two faces
it separates, as in `EdgeSet.edges`; the list is ordered by that face pair. Each `EdgeInfo`
gives the curve type, the total length, a point on the edge, its two ends (or `None` for a
closed edge such as a circle) and how many separate curves the two faces share.

- STEP: the CAD kernel's curves between two different faces. Seams (a curve with the same face
  on both sides, such as the joint line of a cylinder wall) and zero-length curves are left
  out.
- STL: the triangle edges between two Tessera faces, so the same feature angle as the faces.

A line inside a face is not an edge. Add a split line in your CAD program so that the line
becomes the edge between two faces.

Edges are remembered together with the faces, by file content. Raises the same errors as
`describe_faces`.

### `solve`

```python
def solve(
    project: Project,
    workdir: Path | None = None,
    adapter: SolverAdapter | None = None,
    *,
    check_mesh: bool = False,
    timeout: float | None = None,
) -> StudyResult
```

Meshes the part, builds the static problem, runs the solver and summarises the result.

- `workdir`: the folder for solver files (`job.inp`, `job.frd` and CalculiX logs). The folder
  is created if needed. If `None`, a new folder named `tessera-*` is created in the system
  temporary directory. Tessera does not delete it.
- `adapter`: the solver. `None` uses CalculiX (`tessera.solvers.CalculixAdapter`). The adapter
  must support the `static` analysis.
- `check_mesh`: also run the two-mesh check (see "Mesh check" below). The study is solved twice,
  so it takes 3 to 5 times as long. The coarse solve's files go to `workdir / "coarse"`.
- `timeout`: seconds each solver run may take before CalculiX, and the processes it started,
  are stopped and `SolverTimeoutError` is raised. `None` uses `SOLVE_TIMEOUT` (1800 s). The
  limit applies to each solve on its own: with `check_mesh` the two solves may take up to
  twice as long together. Meshing is not counted. It applies to the built-in solver only;
  passing both `adapter` and `timeout` raises `ValueError` (set the limit on your adapter,
  e.g. `CalculixAdapter(timeout=...)`).

v0.1 solves exactly one part with no connections and no placement (identity transform). The
study must have at least one support. The part must have a `material_id` from the built-in
library.

Raises:

- `NotImplementedError`: the project has more than one part, any connection, or a placement
  other than the identity transform. These are planned for later versions.
- `InputError`:
  - the project is invalid, for example a support or load names a part that does not exist;
  - a face number in a support or load that the part does not have;
  - an edge support on two faces that do not share an edge in the mesh;
  - the part has no material, or the material id is unknown;
  - the study has no supports;
  - gravity is set but the material has no density;
  - the geometry is unsupported, or cannot be meshed;
  - the solver adapter cannot run the analysis;
  - `timeout` is not a positive, finite number of seconds.
- `FileNotFoundError`: the part's geometry file does not exist. The message names the part.
- `SolverError` (from `tessera.model`): the CalculiX executable cannot be found or started,
  CalculiX reported an error, or it produced no usable results. `SolverError.log` holds the
  solver output.
- `SolverTimeoutError` (from `tessera.model`, a subclass of `SolverError`): a solver run took
  longer than `timeout` and was stopped. `SolverTimeoutError.timeout` is the limit in seconds;
  `log` holds what the solver wrote before it was stopped. With `check_mesh` the message ends
  with the solve that was stopped, e.g. `(in the coarse solve of the mesh check)`.

Returns a `StudyResult`: the largest displacement, the raw peak von Mises stress, the peak away
from supports and sharp inside edges, and the peak on each face (see "Peak stress" below). With
`check_mesh`, also `mesh_check`.

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
- **Edges.** Each edge is a pair of face numbers. Both faces must exist, be different and meet
  at an edge of the part (see `describe_edges`). An empty list, or an edge listed twice (in
  either order), is refused.
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

#### `add_edge_support`

```python
def add_edge_support(
    project: Project,
    part_id: str,
    edges: Iterable[tuple[int, int]],
    *,
    dofs: tuple[int, ...] = (1, 2, 3),
    support_id: str | None = None,
) -> Support
```

A line support (ADR-006): holds the edges where the listed pairs of faces meet, for example
`[(3, 5)]`. `dofs` works as in `add_support`. Holding all three directions makes a pinned line,
holding one (for example `(3,)`) a roller line. Unlike a face support, the part can still
rotate about the line. Pairs are stored ascending: `(5, 3)` is stored as `(3, 5)`. A pair
stands for everything the two faces share; if they meet along several separate curves, all of
them are held.

Raises `InputError` when:

- the id is empty or used, `edges` is empty, a pair names the same face twice, a face the part
  does not have, or an edge listed twice;
- two faces do not meet at an edge (the message points to `tessera edges`);
- `dofs` is invalid, as for `add_support`.

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
| `Project` | `name: str`, `schema_version: int = 5`, `assembly: Assembly`, `study: Study` | Root object. Extra fields are errors. |
| `Assembly` | `parts: list[PartInstance]`, `connections: list[Connection]` | Part ids must be unique. |
| `PartInstance` | `id: str`, `name: str`, `geometry: Path`, `transform: Transform`, `material_id: str \| None`, `print_setup: PrintSetup \| None` | `geometry` is a STEP or STL path. |
| `Transform` | `rotation`, `translation` | Rotation is a row-major 3 × 3 matrix; translation in mm. Only the identity is solved in v0.1. |
| `PrintSetup` | `process: "fdm" \| "resin"`, `slicer_profile: str \| None`, `build_direction` | Stored only. Not used by the solver in v0.1. |
| `Connection` | `id`, `kind: ConnectionKind`, `a: FaceSet`, `b: FaceSet`, `friction`, `preload_n` | Stored only. Connections are not solved in v0.1. |
| `ConnectionKind` | `"bonded"`, `"contact"`, `"bolt"`, `"pin"` | |
| `FaceSet` | `part_id: str`, `faces: tuple[int, ...]` | At least one face. Face numbers come from `describe_faces`. |
| `EdgeSet` | `part_id: str`, `edges: tuple[tuple[int, int], ...]` | At least one edge, each named by two different faces. Pairs are stored ascending and must be unique. JSON: `"edges": [[3, 5], [5, 6]]`. Pairs come from `describe_edges`. |
| `Study` | `analysis: "static"`, `supports: list[Support]`, `loads: list[Load]`, `gravity: (float, float, float) \| None`, `mesh_size: float \| None` | `gravity` in mm/s². `mesh_size` in mm, greater than 0; `None` lets the mesher choose. |
| `Support` | `id: str`, `faces: FaceSet \| None`, `edges: EdgeSet \| None`, `dofs: tuple[1 \| 2 \| 3, ...] = (1, 2, 3)` | Exactly one of `faces` and `edges` is set. Holds the listed displacement components at zero: 1 = x, 2 = y, 3 = z. No repeats. The `part_id` property gives the part either way. Unused `faces`/`edges` is left out of the JSON. |
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
| `StudyResult` | `mesh: Mesh`, `solution: SolveResult`, `max_displacement: float` (mm), `max_von_mises: float` (MPa), `workdir: Path \| None`, `peak_von_mises: StressSpot`, `clear_peak_von_mises: StressSpot \| None`, `face_peaks: tuple[StressSpot, ...]`, `node_flags: NodeFlags`, `mesh_check: MeshCheck \| None` | Returned by `solve`. `displacement_magnitude` gives the per-node displacement length. See "Peak stress" below. |
| `StressSpot` | `value: float` (MPa), `node: int`, `location: (x, y, z)` (mm), `face: int \| None`, `flags: tuple[StressFlag, ...]`, `wall: WallSpan \| None` | A stress at one node. `face` is the face the node lies on (the lowest number on an edge between faces), `None` inside the part. `flagged` is `True` when `flags` is not empty. `wall` is the wall at the node (see "Thin walls" below); `None` where it was not measured. |
| `WallSpan` | `thickness: float` (mm), `elements: float`, `direction: (x, y, z)` | How thick the part is at a node: the diameter of the largest ball inside the part touching the surface there. `elements` = `thickness` / the local median corner edge length of the mesh it was measured on. `direction` is the unit direction from the node to the ball's centre. |
| `StressFlag` | `kind: FlagKind`, `support: str \| None` | Why a node is flagged. `support` is the support's id for `FlagKind.SUPPORT`. `str()` gives `near support 'root'` or `near a sharp inside edge`. |
| `FlagKind` | `"support"`, `"sharp-inside-edge"` | |
| `NodeFlags` | `masks: dict[StressFlag, (N,) bool]`, `n_nodes: int`, `distance: float` (mm), `edge_length: float` (mm) | One mask per reason that flags at least one node. `flagged` is the (N,) union; `at(node)` lists the reasons for one node. `edge_length` is h, `distance` is d = 2 h. |
| `SolveResult` | `displacements: (N, 3)` mm, `stresses: (N, 6)` MPa in the order xx, yy, zz, xy, yz, zx, `solver: str`, `log: str`, `metadata: dict[str, str]`, `point_stresses: PointStresses \| None` | Per node, 0-based. Nodes that no element uses are NaN. `point_stresses` is `None` unless the `SolveJob` asked for them (`point_stresses=True`, keyword-only, off by default: it costs about 0.6 kB per element in the solver's `.dat` file plus parse time) and the solver provides them (CalculiX does, rule `tet4`); when present, `metadata["point_rule"]` names their integration rule (e.g. `tet4`). |
| `PointStresses` | `points: (M, Q, 3)` mm, `stresses: (M, Q, 6)` MPa in the order xx, yy, zz, xy, yz, zx | Each element's own stresses at its Q integration points (ADR-008), not extrapolated or averaged; row j is element j. `n_elements` and `n_points` properties. Missing entries are NaN. Points lie inside the element, so their maximum reads low at surfaces: every reported peak stays nodal. |
| `Mesh` | `nodes: (N, 3)` mm, `elements: (M, 10)` quadratic tetrahedra (C3D10 order), `element_sets`, `node_sets`, `surfaces` | `n_nodes` and `n_elements` properties. |
| `EdgeInfo` | `faces: (int, int)`, `kind: str`, `length: float` (mm), `start: (x, y, z) \| None`, `end: (x, y, z) \| None`, `point: (x, y, z)`, `curves: int` | `kind` is the CAD curve type (`"Line"`, `"Circle"`, ..., `"Mixed"` when the curves differ) for STEP, and for STL `"Line"` when every curve is straight, otherwise `"Curved"`. `start` and `end` are `None` for a closed edge or several separate curves. `point` lies on the edge (the middle of a single straight edge). |
| `FaceInfo` | `index: int`, `kind: str`, `area: float` (mm²), `centroid: (x, y, z)` (mm), `normal: (x, y, z) \| None` | `kind` is the CAD surface type (`"Plane"`, `"Cylinder"`, ...) for STEP, and `"Plane"` or `"Curved"` for STL. `normal` is `None` for curved faces. |

#### Peak stress

In a linear model the stress has no finite value at the edges of a support that holds two or
three directions, along a line support and at sharp inside edges, so the nodal stress there keeps rising as the mesh
is refined. `solve` flags those places (ADR-007):

- **Supports:** nodes within d of a node of a face support that holds two or three directions.
  A face support that holds one direction (a symmetry plane or roller) is never flagged.
- **Sharp inside edges:** nodes within d of an edge of the outer surface that turns inward by
  30° or more, and edges shared by more than two surface faces (where two bodies touch along
  a line). Outward (convex) edges, smooth fillets and holes are not flagged; the turn is
  measured on the curved element faces, so a coarse mesh of a round hole does not count.
- **Line supports** (ADR-006): nodes within d of a line support, whatever directions it holds.
- d = 2 h, where h is the median corner-to-corner element edge length (`node_flags.edge_length`).

The fields:

- `max_von_mises`: the largest nodal von Mises stress anywhere, flagged or not. Its meaning has
  not changed.
- `peak_von_mises`: the same value as a `StressSpot`, with its node, location, face and flags.
- `clear_peak_von_mises`: the largest von Mises stress at an unflagged node, "away from
  supports and sharp inside edges"; `None` if no unflagged node has a result. It is not a
  converged stress and can read low: where the highest stress is expected at a support (a
  clamped base under bending), it is taken about d from the clamp and misses the rise toward
  it; the true stress there lies between `clear_peak_von_mises` and `peak_von_mises`. Do not
  quote it alone against a strength in that case; compare two mesh sizes instead.
- `face_peaks`: the largest von Mises stress on every face that has a result, highest first,
  each with its own flags.
- `node_flags`: the per-node flags, for use with any other stress result.

Nodes that no element uses have no result (NaN) and are skipped. Peaks under point-like loads
and at the edges of loaded faces are not flagged yet.

#### Integration-point stresses

`SolveResult.point_stresses` (ADR-008) holds each element's own stresses at its integration
points. `tessera.post` works with them without knowing the solver:

- `element_nodal_stresses(mesh, point_stresses)`: (M, 10, 6) stresses at each element's ten
  nodes, from the affine field through that element's points in physical coordinates. Not
  averaged: a node shared by two elements has a value in each. Exact for a linear stress field.
  From CalculiX's printed stresses (7 significant digits) it is exact to within about
  1.4e-6 of the element's largest stress. Raises `ValueError` if the points do not match the
  mesh or there are fewer than four per element. An element with a NaN or infinite entry gets
  NaN; so does an element whose points lie in one plane, and a `DegenerateElementsWarning`
  (its `elements` attribute) lists those.
- `point_von_mises(point_stresses)`: (M, Q) von Mises stress at every point. `von_mises` itself
  takes any (..., 6) array and returns (...).
- `flag_points(job, points)`: the flags of ADR-007 at any (P, 3) points, e.g.
  `point_stresses.points.reshape(-1, 3)`; `flag_nodes(job)` is `flag_points` at the mesh nodes.

None of the reported numbers use them: peaks, flags and the mesh check stay nodal. Integration
points lie inside the elements, so a maximum over them reads low at a surface.

#### Mesh check

`solve(project, check_mesh=True)` solves the study at the project's `mesh_size` (or gmsh's own
sizes when it is `None`), then again with every element size gmsh uses scaled by
`MESH_CHECK_REFINEMENT` (1/√2), including the small elements it places at small holes and short
curves. The returned `StudyResult` is the fine solve, so every number in it comes from the fine
mesh, and `mesh_check` compares them with the coarse solve (ADR-007):

| Type | Fields | Notes |
|---|---|---|
| `MeshCheck` | `coarse_edge_length`, `fine_edge_length` (mm), `coarse_nodes: int`, `threshold: float`, `refinement_limit: float`, `max_displacement: SpotCheck`, `peak_von_mises: SpotCheck`, `clear_peak_von_mises: SpotCheck \| None`, `face_peaks: dict[int, SpotCheck]`, `wall_elements_min: float` | Edge lengths are each mesh's median corner edge length. `clear_peak_von_mises` is `None` unless both meshes have a clear peak; `face_peaks` has the faces with a peak on both meshes. `wall_elements_min` is the thin-wall minimum used (5). |
| `SpotCheck` | `coarse: float`, `fine: float`, `change: float`, `refinement: float`, `verdict: Verdict`, `elements_across: float \| None` | `change` = \|fine − coarse\| / \|fine\|; `inf` if only the fine value is zero. `refinement` is the local element edge length, fine / coarse, where the fine value was taken (for the displacement, the ratio of the median edge lengths). `elements_across` is the fine spot's wall thickness over the coarse local edge length; `None` for the displacement or where no wall was measured. |
| `Verdict` | `SETTLED = "settled"`, `MESH_DEPENDENT = "mesh-dependent"`, `NOT_REFINED = "not-refined"`, `THIN_WALL = "thin-wall"` | `NOT_REFINED` if `refinement > refinement_limit` (or it could not be measured); otherwise, if `change <= threshold`, `THIN_WALL` for a flagged spot with `elements_across < wall_elements_min` and settled for any other; otherwise mesh-dependent. `NOT_REFINED` and `THIN_WALL` are no verdict: the spot's flags still apply. |

`threshold` is `tessera.post.MESH_CHECK_THRESHOLD` (5 %) and `refinement_limit` is
`tessera.post.MESH_CHECK_REFINEMENT_LIMIT` (0.85). A number the fine mesh did not refine around
is not judged: an unchanged value there says nothing, and its flags still apply.
`tessera.post.compare_meshes(coarse, fine, threshold, refinement_limit, wall_elements_min)`,
`spot_check` and
`local_edge_length(mesh, point)` (the median corner edge length of the elements at the node
nearest the point) do the comparison for results you solved yourself.

Numbers are compared like with like, by value: the raw peak with the raw peak, and each face's
peak with the same face's peak, wherever on it they sit. If two places have nearly the same
peak, the two meshes may put it in different places; compare the locations. The clear peak is
taken outside a flag band of 2 h, which narrows on the fine mesh, so part of its change can come
from the band rather than from convergence.

A verdict is measured, so it replaces the flags: report a flagged spot that settled, not an
unflagged one that is mesh-dependent.

#### Thin walls

On a mesh too coarse to resolve a clamp (two elements through a beam's depth in TS-01) the
singular peak there barely moves between the meshes and would read settled. So `solve` measures
the wall (`StressSpot.wall`) at the raw peak, the clear peak and every flagged spot, and the
check judges a flagged spot only where the coarse mesh had at least `WALL_ELEMENTS_MIN` (5,
approved by the maintainer) elements across it; a flagged spot that changed within the threshold
on a thinner wall, or whose wall could not be measured, is `THIN_WALL` (ADR-007 Amendment 2).
Unflagged spots are never guarded: a smooth stress is represented well with one or two quadratic
elements across. A `mesh-dependent` reading stays mesh-dependent whatever the wall. Unflagged
face peaks are not measured (`wall` is `None`).

The wall is the diameter of the largest ball inside the part that touches the surface at the
node (a shrinking ball): its centre moves in along a face's inward normal and its radius shrinks
to r = \|p − q\|² / (2 d · (q − p)), the ball through the node p and the surface point q
nearest the centre, until no surface point is inside it. Where faces meet at an edge or a
corner (normals more than 30° apart), each face gets its own ball and the smallest counts. The
surface is the mesh's outer element faces, each TRI6 split into four flat triangles; their
centroids are the samples, in one k-d tree per mesh, and q is refined to the sample triangle's
plane. A sample stops the ball only where the ball would leave the part through it, so a
neighbouring face at a convex edge (which the ball meets edge-on) and the other face at an inside
edge do not collapse it. So a plate reads its thickness, a bar's end face its section (not its
length), and a node near another face less than the wall: the ball touches that face too. At a
spot flagged by a face support, the samples on that support's faces are left out
(`samples_on(surface, nodes)`, passed as `wall_span(..., ignore=)`): the part continues into the
fixture there, so a peak next to a clamped face reads the section, not twice its distance to
the clamp. Line supports hold no face, so nothing is left out for them. A node inside the part
is measured at the nearest surface node. A support never doubles the wall, not even a
one-direction one: it cannot be told from a roller on a table.

`WallSpan.elements` uses the mesh the spot comes from (the fine one after a check);
`SpotCheck.elements_across` uses the coarse mesh, which decides the guard.

`tessera.core.WALL_ELEMENTS_MIN` is the same value as `tessera.post.WALL_ELEMENTS_MIN`, for
clients that import only `tessera.core`. `tessera.post.outer_surface(mesh)` builds the surface
once, and `tessera.post.wall_span(surface, node)` measures one node (`None` if no ball is stopped,
as on an open surface). `tessera.post.node_edge_length(mesh, node)` is `local_edge_length` for a
known node.

In the `SolveJob` that physics builds, each support is a `FixedSupport` whose node set is named
`support:<id>`. A line support has `FixedSupport.line = True`, so results code can find the
supports whose peak stress is mesh-dependent.

### Errors

| Type | Module | Meaning |
|---|---|---|
| `InputError` | `tessera.model` | A problem in what the user asked for (a subclass of `ValueError`). Its message is safe to show as-is. |
| `SolverError` | `tessera.model` | The solver failed. Has a `log` attribute. Subclass of `RuntimeError`. |
| `SolverTimeoutError` | `tessera.model` | The solver took longer than its timeout and was stopped. Has `timeout` (s) and `log`. Subclass of `SolverError`. |

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

### Solve-job materials and orientation

These types are in `tessera.model` (ADR-009). Projects cannot name an orthotropic material yet;
the CalculiX adapter writes one as engineering constants and a section's `MaterialFrame` as an
`*ORIENTATION`, and returns displacements and stresses in global (part) axes.

| Type | Fields | Notes |
|---|---|---|
| `IsotropicMaterial` | `name`, `youngs_modulus`, `poisson_ratio`, `density = None` | MPa; density in t/mm³. |
| `OrthotropicMaterial` | `name`, `e1`, `e2`, `e3`, `nu12`, `nu13`, `nu23`, `g12`, `g13`, `g23`, `density = None` | Axes 1 = raster, 2 = across the road, 3 = build. Major Poisson ratios: nu_ij = -eps_j / eps_i under uniaxial sigma_i; nu_ji = nu_ij E_j / E_i is derived. Raises `ValueError` unless every E and G is positive and finite, every Poisson ratio is finite, and the compliance is positive definite. |
| `Material` | | The union `IsotropicMaterial \| OrthotropicMaterial`; `SolveJob.materials` is a `tuple[Material, ...]`. |
| `MaterialFrame` | `axis_1`, `axis_2` | Material axes 1 and 2 in part coordinates, unit and orthogonal within 1e-9. `axis_3` = `axis_1` × `axis_2`; `rotation` is the 3 × 3 matrix whose rows are axes 1, 2, 3; `MaterialFrame.from_axes(a, b)` builds one from two non-parallel vectors (Gram-Schmidt). |
| `Section` | `element_set`, `material`, `orientation: MaterialFrame \| None = None` | `None` means the part axes. |

`tessera.model.elasticity` has the pure functions: `orthotropic_compliance(e1, ..., g23)` returns
the 6 × 6 compliance in Voigt order xx, yy, zz, xy, yz, zx (the order of `SolveResult.stresses`)
with engineering shear strains; `check_orthotropic(e1, ..., g23)` raises the `ValueError` above;
`rotate_stress(stress, frame)` takes (6,) or (N, 6) stresses in part axes and returns them in the
frame's material axes (11, 22, 33, 12, 23, 31).

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
clear = result.clear_peak_von_mises
if clear is not None:
    print(f"{clear.value:.1f} MPa away from supports and sharp inside edges, face {clear.face}")
for flag in result.peak_von_mises.flags:
    print(f"the raw peak is {flag}")

save_project(project, Path("bar-with-study.tessera"))
```

`normal == (-1.0, 0.0, 0.0)` compares floats exactly. For real use, compare with a tolerance,
as the benchmarks do.
