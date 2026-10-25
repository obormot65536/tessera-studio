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
class WallSpan:
    """How thick the part is at a node, and in elements (ADR-007 Amendment 2)."""

    thickness: float
    """Diameter of the largest ball inside the part that touches the surface at the node, mm:
    the smallest over the faces meeting there. Never doubled across a support."""
    elements: float
    """``thickness`` / the local median corner edge length of this mesh at the node."""
    direction: Vec3
    """Unit direction from the node to the centre of that ball (the face's inward normal)."""


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
    wall: WallSpan | None = None
    """The wall at the node (for a node inside the part, at the nearest surface node); ``None``
    where it was not measured or could not be (an open surface)."""

    @property
    def flagged(self) -> bool:
        return bool(self.flags)


class FailureCriterion(StrEnum):
    """How a material's strength ratio is computed (ADR-010 D2)."""

    TSAI_WU = "tsai-wu"
    """Orthotropic materials; interaction coefficients F*_ij = -1/2."""
    VON_MISES = "von-mises"
    """Ductile isotropic materials, against the yield strength."""
    MAX_PRINCIPAL = "max-principal"
    """Brittle isotropic materials (resin), against the tensile and compressive strengths."""


class FailureMode(StrEnum):
    """Which way a spot is predicted to fail first (ADR-010 D2).

    For Tsai-Wu it names the stress term with the largest maximum-stress ratio there; it
    explains the strength ratio and never changes it.
    """

    ALONG_RASTER_TENSION = "along-raster-tension"
    ALONG_RASTER_COMPRESSION = "along-raster-compression"
    ACROSS_ROADS_TENSION = "across-roads-tension"
    """Across the roads within a layer (material axis 2)."""
    ACROSS_ROADS_COMPRESSION = "across-roads-compression"
    ACROSS_LAYERS_TENSION = "across-layers-tension"
    """Along the build direction (material axis 3): the layers pulling apart."""
    ACROSS_LAYERS_COMPRESSION = "across-layers-compression"
    IN_LAYER_SHEAR = "in-layer-shear"
    """tau_12."""
    INTERLAYER_SHEAR = "interlayer-shear"
    """tau_13 or tau_23: the layers sliding on each other."""
    VON_MISES = "von-mises"
    MAX_PRINCIPAL_TENSION = "max-principal-tension"
    MAX_PRINCIPAL_COMPRESSION = "max-principal-compression"


@dataclass(frozen=True)
class StrengthSpot:
    """A strength ratio at one node, with where it is, why and whether it is flagged."""

    ratio: float
    """How many times every load of the study (gravity included) could be multiplied before the
    criterion predicts first failure here; ``inf`` where the stress is zero (ADR-010 D1)."""
    node: int
    """0-based node index."""
    element: int
    """0-based index of the element whose stress gave the lowest ratio at the node."""
    location: Vec3
    """mm"""
    face: int | None
    """Geometry face the node lies on, as for ``StressSpot.face``."""
    mode: FailureMode
    flags: tuple[StressFlag, ...] = ()
    wall: WallSpan | None = None

    @property
    def flagged(self) -> bool:
        return bool(self.flags)


@dataclass(frozen=True)
class StrengthResult:
    """The strength check of a study (ADR-010 D5), from element-nodal stresses (ADR-008)."""

    criteria: dict[str, FailureCriterion]
    """Criterion used for each checked material, by material name."""
    lowest: StrengthSpot
    """The lowest ratio in the part, flagged or not."""
    clear_lowest: StrengthSpot | None
    """The lowest ratio away from supports and sharp inside edges (unflagged nodes); ``None``
    if every checked node is flagged. The headline unless the mesh check settles ``lowest``."""
    face_lowest: tuple[StrengthSpot, ...]
    """Lowest ratio on each geometry face, lowest first."""
    assumed: dict[str, tuple[str, ...]] = field(default_factory=dict)
    """Strengths taken equal to the tensile ones, by material name (e.g. ``("c1", "c3")``)."""
    unchecked: tuple[str, ...] = ()
    """Materials without strength data; their elements are not checked."""


class Verdict(StrEnum):
    """What the two-mesh check found for one reported number (ADR-007)."""

    SETTLED = "settled"
    """It changed by at most the threshold between the two meshes."""
    MESH_DEPENDENT = "mesh-dependent"
    """It changed by more than the threshold: do not compare it with a strength."""
    NOT_REFINED = "not-refined"
    """The fine mesh was not fine enough there to judge it: no verdict, the flags still apply."""
    THIN_WALL = "thin-wall"
    """A flagged spot that changed by at most the threshold, but the coarse mesh had fewer than
    ``MeshCheck.wall_elements_min`` elements across the wall there, too few to show the rise,
    or the wall could not be measured: no verdict, the flags still apply (ADR-007 Amendment 2)."""


@dataclass(frozen=True)
class SpotCheck:
    """One reported number at the coarse and the fine mesh of a two-mesh check."""

    coarse: float
    fine: float
    """The value reported in the result."""
    change: float
    """|fine - coarse| / |fine|; ``inf`` if only the fine value is zero."""
    refinement: float
    """How much finer the fine mesh is there: local element edge length, fine / coarse."""
    verdict: Verdict
    elements_across: float | None = field(default=None, kw_only=True)
    """Elements of the coarse mesh across the wall at the spot: the fine spot's wall thickness
    / the coarse local edge length; ``None`` if no wall was measured there (a flagged spot that
    changed within the threshold is then ``thin-wall``: the guard fails closed)."""


@dataclass(frozen=True)
class MeshCheck:
    """The two-mesh check of a study (ADR-007): solved again with every element size scaled by
    1/sqrt(2).

    Every number in the ``StudyResult`` carrying this check comes from the fine mesh.
    """

    coarse_edge_length: float
    """Median corner edge length h of the coarse mesh, mm."""
    fine_edge_length: float
    """The same for the fine mesh, mm."""
    coarse_nodes: int
    threshold: float
    """Largest relative change that still counts as settled."""
    refinement_limit: float
    """Largest local refinement ratio (fine / coarse edge length) at which a number is judged."""
    max_displacement: SpotCheck
    peak_von_mises: SpotCheck
    clear_peak_von_mises: SpotCheck | None
    """``None`` if either mesh has no clear peak."""
    face_peaks: dict[int, SpotCheck]
    """By geometry face, for faces with a peak on both meshes."""
    lowest_strength_ratio: SpotCheck | None = field(default=None, kw_only=True)
    """The lowest strength ratio (``StrengthResult.lowest``) on both meshes; ``None`` without a
    strength check (ADR-010 D5)."""
    wall_elements_min: float = field(default=5.0, kw_only=True)
    """Fewest coarse elements across the wall at which a flagged spot is judged settled; the
    default is ``tessera.post.WALL_ELEMENTS_MIN`` (a post test keeps them equal)."""


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
    mesh_check: MeshCheck | None = field(default=None, kw_only=True)
    """The two-mesh check, if it was asked for (``solve(..., check_mesh=True)``)."""
    orthotropic: bool = field(default=False, kw_only=True)
    """True if any material in the study is orthotropic (ADR-009 D4). Von Mises still shows where
    the stress is, but it is not a strength check for such a part."""
    strength: StrengthResult | None = field(default=None, kw_only=True)
    """Strength ratios (ADR-010), when a material of the study has strength data."""

    @property
    def displacement_magnitude(self) -> np.ndarray:
        return np.linalg.norm(self.solution.displacements, axis=1)
