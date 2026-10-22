"""Project data model: an assembly of part instances joined by connections.

A single part is an assembly with one instance, so every code path handles
assemblies. All lengths are in mm (see ``tessera.model.units``).
"""

from __future__ import annotations

from enum import StrEnum
from pathlib import Path
from typing import Annotated, Any, Literal

from pydantic import (
    BaseModel,
    ConfigDict,
    Discriminator,
    Field,
    SerializerFunctionWrapHandler,
    Tag,
    field_validator,
    model_serializer,
    model_validator,
)

from tessera.model.elasticity import check_orthotropic

SCHEMA_VERSION = 6
"""2: Support.dofs (directional supports). 3: PressureLoad, Study.gravity, load ``kind``.
4: MomentLoad. 5: Support.edges (line supports, ADR-006).
6: Project.materials, PrintSetup.raster_angle (ADR-009). Older files load unchanged."""

Vec3 = tuple[float, float, float]


class _Model(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=False, allow_inf_nan=False)


class Transform(_Model):
    """Rigid placement of a part: rotation (row-major 3x3) then translation in mm."""

    rotation: tuple[Vec3, Vec3, Vec3] = ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))
    translation: Vec3 = (0.0, 0.0, 0.0)


class PrintSetup(_Model):
    """How a printed part is made. Filled from a slicer project in v0.5."""

    process: Literal["fdm", "resin"]
    slicer_profile: str | None = None
    build_direction: Vec3 = (0.0, 0.0, 1.0)
    raster_angle: float = 0.0
    """Degrees, counter-clockwise looking down the build direction (ADR-009 D2)."""


_Density = Annotated[float, Field(gt=0.0)] | None
"""t/mm^3; ``None`` when unknown (gravity then refuses the material)."""


class IsotropicDefinition(_Model):
    """A linear isotropic material defined in the project (ADR-009 D3). Moduli in MPa."""

    kind: Literal["isotropic"] = "isotropic"
    id: str
    youngs_modulus: Annotated[float, Field(gt=0.0)]
    poisson_ratio: Annotated[float, Field(gt=-1.0, lt=0.5)]
    density: _Density = None


class OrthotropicDefinition(_Model):
    """A linear orthotropic material defined in the project (ADR-009 D1, D3).

    Engineering constants in MPa on the material axes (1 = raster, 2 = across the road in the
    layer, 3 = build), with the major Poisson ratios nu_ij = -eps_j / eps_i under uniaxial
    sigma_i. Validated with ``tessera.model.elasticity.check_orthotropic``.
    """

    kind: Literal["orthotropic"] = "orthotropic"
    id: str
    e1: float
    e2: float
    e3: float
    nu12: float
    nu13: float
    nu23: float
    g12: float
    g13: float
    g23: float
    density: _Density = None

    @model_validator(mode="after")
    def _check_constants(self) -> OrthotropicDefinition:
        try:
            check_orthotropic(
                self.e1,
                self.e2,
                self.e3,
                self.nu12,
                self.nu13,
                self.nu23,
                self.g12,
                self.g13,
                self.g23,
            )
        except ValueError as exc:
            raise ValueError(f"material {self.id!r}: {exc}") from None
        return self


MaterialDefinition = Annotated[
    IsotropicDefinition | OrthotropicDefinition, Field(discriminator="kind")
]


class PartInstance(_Model):
    id: str
    name: str
    geometry: Path
    transform: Transform = Field(default_factory=Transform)
    material_id: str | None = None
    print_setup: PrintSetup | None = None
    """``None`` means the part is not printed (e.g. a steel bolt)."""


class ConnectionKind(StrEnum):
    BONDED = "bonded"
    CONTACT = "contact"
    BOLT = "bolt"
    PIN = "pin"


class FaceSet(_Model):
    """Faces of one part instance, by geometry face index."""

    part_id: str
    faces: Annotated[tuple[int, ...], Field(min_length=1)]


class EdgeSet(_Model):
    """Edges of one part instance, each named by the two faces it separates (ADR-006).

    A pair stands for everything the two faces share; if they meet along several separate
    curves, all of them are meant. Each pair is stored ascending, e.g. ``(3, 5)``.
    """

    part_id: str
    edges: Annotated[tuple[tuple[int, int], ...], Field(min_length=1)]

    @field_validator("edges")
    @classmethod
    def _ascending_unique_pairs(
        cls, edges: tuple[tuple[int, int], ...]
    ) -> tuple[tuple[int, int], ...]:
        pairs: list[tuple[int, int]] = []
        for a, b in edges:
            if a == b:
                raise ValueError(f"an edge lies between two different faces, got {a}-{b}")
            pair = (min(a, b), max(a, b))
            if pair in pairs:
                raise ValueError(f"edge {pair[0]}-{pair[1]} is listed twice")
            pairs.append(pair)
        return tuple(pairs)


class Connection(_Model):
    id: str
    kind: ConnectionKind
    a: FaceSet
    b: FaceSet
    friction: Annotated[float, Field(ge=0.0)] | None = None
    preload_n: Annotated[float, Field(ge=0.0)] | None = None


class Assembly(_Model):
    parts: list[PartInstance] = Field(default_factory=list)
    connections: list[Connection] = Field(default_factory=list)

    @model_validator(mode="after")
    def _check_references(self) -> Assembly:
        ids = [p.id for p in self.parts]
        if len(ids) != len(set(ids)):
            raise ValueError("part instance ids must be unique")
        known = set(ids)
        for c in self.connections:
            for side in (c.a, c.b):
                if side.part_id not in known:
                    raise ValueError(f"connection {c.id!r} refers to unknown part {side.part_id!r}")
        return self


class Support(_Model):
    """Support: the listed global displacement components of the faces or edges are held at zero.

    Exactly one of ``faces`` and ``edges`` is set. The default fixes all three components (a
    fixed support, or a pinned line on edges); a single component normal to a plane face makes
    a symmetry or roller support, and on edges a roller line (ADR-006).
    """

    id: str
    faces: FaceSet | None = None
    edges: EdgeSet | None = None
    dofs: Annotated[tuple[Literal[1, 2, 3], ...], Field(min_length=1)] = (1, 2, 3)
    """Global components held at zero: 1 = x, 2 = y, 3 = z."""

    @property
    def part_id(self) -> str:
        """The part whose faces or edges are held."""
        held = self.faces or self.edges
        assert held is not None  # the validator guarantees one of them
        return held.part_id

    @model_validator(mode="after")
    def _check_support(self) -> Support:
        if (self.faces is None) == (self.edges is None):
            raise ValueError(f"support {self.id!r} holds faces or edges: give exactly one")
        if len(set(self.dofs)) != len(self.dofs):
            raise ValueError(f"support {self.id!r} lists a component twice: {self.dofs}")
        return self

    @model_serializer(mode="wrap")
    def _omit_unused_set(self, handler: SerializerFunctionWrapHandler) -> dict[str, Any]:
        # A face support is written exactly as in schema 4, with no "edges": null.
        data: dict[str, Any] = handler(self)
        for key in ("faces", "edges"):
            if data.get(key) is None:
                data.pop(key, None)
        return data


class ForceLoad(_Model):
    """Total force in N, spread uniformly (per unit area) over the faces."""

    kind: Literal["force"] = "force"
    id: str
    faces: FaceSet
    force: Vec3


class PressureLoad(_Model):
    """Uniform pressure in MPa normal to the faces; positive pushes into the part."""

    kind: Literal["pressure"] = "pressure"
    id: str
    faces: FaceSet
    pressure: float


class MomentLoad(_Model):
    """Moment in N mm about the area centroid of the faces, e.g. a torque on a shaft end.

    It is spread like a distributing coupling: nodal forces proportional to each node's area
    share times its distance from the centroid, statically equivalent to the moment, with zero
    net force and no added stiffness.
    """

    kind: Literal["moment"] = "moment"
    id: str
    faces: FaceSet
    moment: Vec3


def _load_kind(value: object) -> str:
    # Schema-2 files have force loads without a ``kind``.
    if isinstance(value, dict):
        return str(value.get("kind", "force"))  # pyright: ignore[reportUnknownMemberType]
    return str(getattr(value, "kind", ""))


Load = Annotated[
    Annotated[ForceLoad, Tag("force")]
    | Annotated[PressureLoad, Tag("pressure")]
    | Annotated[MomentLoad, Tag("moment")],
    Discriminator(_load_kind),
]


class Study(_Model):
    analysis: Literal["static"] = "static"
    supports: list[Support] = Field(default_factory=list)
    loads: list[Load] = Field(default_factory=list)
    gravity: Vec3 | None = None
    """Acceleration due to gravity in mm/s^2 acting on every part, e.g. (0, 0, -9810)."""
    mesh_size: Annotated[float, Field(gt=0.0)] | None = None
    """Target element size in mm; ``None`` lets the mesher choose."""


class Project(_Model):
    schema_version: int = SCHEMA_VERSION
    name: str
    assembly: Assembly = Field(default_factory=Assembly)
    study: Study = Field(default_factory=Study)
    materials: list[MaterialDefinition] = Field(default_factory=list)
    """Materials defined in this project; at solve time they win over the library."""

    @model_validator(mode="after")
    def _check_material_ids(self) -> Project:
        seen: set[str] = set()
        for material in self.materials:
            if material.id in seen:
                raise ValueError(f"material id {material.id!r} is used twice")
            seen.add(material.id)
        return self

    @model_validator(mode="after")
    def _check_study_references(self) -> Project:
        known = {p.id for p in self.assembly.parts}
        used = [(s.id, s.part_id) for s in self.study.supports]
        used += [(load.id, load.faces.part_id) for load in self.study.loads]
        for item_id, part_id in used:
            if part_id not in known:
                raise ValueError(f"{item_id!r} refers to unknown part {part_id!r}")
        return self
