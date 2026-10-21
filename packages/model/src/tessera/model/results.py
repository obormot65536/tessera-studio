"""Results of a solved study, as returned by ``tessera.core.solve``."""

from dataclasses import dataclass, field
from enum import StrEnum
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

from tessera.model.mesh import Mesh
from tessera.model.solver import SolveResult

Vec3 = tuple[float, float, float]


class FlagKind(StrEnum):
    """Why the stress at a node is known to depend on the mesh (ADR-007)."""

    SUPPORT = "support"
    """Near a line support (ADR-006), or a face support that holds two or three directions."""
    SHARP_INSIDE_EDGE = "sharp-inside-edge"
    """Near an edge where the outer surface turns inward by 30 degrees or more."""


@dataclass(frozen=True)
class StressFlag:
    """One reason a node's stress is flagged."""

    kind: FlagKind
    support: str | None = None
    """The support's id for ``FlagKind.SUPPORT``, otherwise ``None``."""

    def __str__(self) -> str:
        if self.kind is FlagKind.SUPPORT:
            return f"near support {self.support!r}"
        return "near a sharp inside edge"


@dataclass(frozen=True, eq=False)
class NodeFlags:
    """Per-node flags: nodes whose stress keeps changing as the mesh is refined (ADR-007).

    The flags depend only on the mesh and the supports, not on the stress, so every result
    (von Mises now, safety factors later) uses the same ones.
    """

    masks: dict[StressFlag, NDArray[np.bool_]]
    """One (N,) mask per reason that flags at least one node."""
    n_nodes: int
    distance: float
    """Flag distance d = 2 h, mm: a node within d of a support or inside edge is flagged."""
    edge_length: float
    """h: the median corner-to-corner element edge length, mm."""

    @property
    def flagged(self) -> NDArray[np.bool_]:
        """(N,) True where any reason applies."""
        out = np.zeros(self.n_nodes, dtype=bool)
        for mask in self.masks.values():
            out |= mask
        return out

    def at(self, node: int) -> tuple[StressFlag, ...]:
        """The reasons that apply to one node, in the order of ``masks``."""
        return tuple(flag for flag, mask in self.masks.items() if mask[node])


@dataclass(frozen=True)
class StressSpot:
    """A stress value at one node, with where it is and whether it is flagged."""

    value: float
    """MPa"""
    node: int
    """0-based node index."""
    location: Vec3
    """mm"""
    face: int | None
    """Geometry face the node lies on (the lowest number on an edge between faces), or
    ``None`` inside the part."""
    flags: tuple[StressFlag, ...] = ()

    @property
    def flagged(self) -> bool:
        return bool(self.flags)


@dataclass(frozen=True)
class StudyResult:
    mesh: Mesh
    solution: SolveResult
    max_displacement: float
    """Largest nodal displacement magnitude, mm."""
    max_von_mises: float
    """Largest nodal von Mises stress, MPa (the raw peak, flagged or not)."""
    workdir: Path | None = None
    """Where the solver's input and output files were kept."""
    peak_von_mises: StressSpot = field(kw_only=True)
    """The raw peak: where ``max_von_mises`` is, and its flags."""
    clear_peak_von_mises: StressSpot | None = field(kw_only=True)
    """Largest von Mises stress away from supports and sharp inside edges (unflagged nodes);
    ``None`` if every node is flagged."""
    face_peaks: tuple[StressSpot, ...] = field(kw_only=True)
    """Largest von Mises stress on each geometry face, highest first."""
    node_flags: NodeFlags = field(kw_only=True)

    @property
    def displacement_magnitude(self) -> np.ndarray:
        return np.linalg.norm(self.solution.displacements, axis=1)
