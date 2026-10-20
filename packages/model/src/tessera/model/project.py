"""Project data model: an assembly of part instances joined by connections.

A single part is an assembly with one instance, so every code path handles
assemblies. All lengths are in mm (see ``tessera.model.units``).
"""

from __future__ import annotations

from enum import StrEnum
from pathlib import Path
from typing import Annotated, Literal

from pydantic import BaseModel, ConfigDict, Discriminator, Field, Tag, model_validator

SCHEMA_VERSION = 3
"""2: Support.dofs (directional supports). 3: PressureLoad, Study.gravity, load ``kind``.
Older files load unchanged."""

Vec3 = tuple[float, float, float]


class _Model(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=False)


class Transform(_Model):
    """Rigid placement of a part: rotation (row-major 3x3) then translation in mm."""

    rotation: tuple[Vec3, Vec3, Vec3] = ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))
    translation: Vec3 = (0.0, 0.0, 0.0)


class PrintSetup(_Model):
    """How a printed part is made. Filled from a slicer project in v0.5."""

    process: Literal["fdm", "resin"]
    slicer_profile: str | None = None
    build_direction: Vec3 = (0.0, 0.0, 1.0)


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
    """Support: the listed global displacement components of the faces are held at zero.

    The default fixes all three (a fixed support); a single component normal to a plane face
    makes a symmetry or roller support.
    """

    id: str
    faces: FaceSet
    dofs: Annotated[tuple[Literal[1, 2, 3], ...], Field(min_length=1)] = (1, 2, 3)
    """Global components held at zero: 1 = x, 2 = y, 3 = z."""

    @model_validator(mode="after")
    def _check_unique_dofs(self) -> Support:
        if len(set(self.dofs)) != len(self.dofs):
            raise ValueError(f"support {self.id!r} lists a component twice: {self.dofs}")
        return self


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


def _load_kind(value: object) -> str:
    # Schema-2 files have force loads without a ``kind``.
    if isinstance(value, dict):
        return str(value.get("kind", "force"))  # pyright: ignore[reportUnknownMemberType]
    return str(getattr(value, "kind", ""))


Load = Annotated[
    Annotated[ForceLoad, Tag("force")] | Annotated[PressureLoad, Tag("pressure")],
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

    @model_validator(mode="after")
    def _check_study_references(self) -> Project:
        known = {p.id for p in self.assembly.parts}
        for item in [*self.study.supports, *self.study.loads]:
            if item.faces.part_id not in known:
                raise ValueError(f"{item.id!r} refers to unknown part {item.faces.part_id!r}")
        return self
