"""Shared data model, units and interfaces (bottom layer)."""

from tessera.model.errors import InputError
from tessera.model.geometry_info import FaceInfo, TriangleSurface
from tessera.model.mesh import TET10_EDGES, TRI6_EDGES, Mesh
from tessera.model.project import (
    SCHEMA_VERSION,
    Assembly,
    Connection,
    ConnectionKind,
    FaceSet,
    ForceLoad,
    Load,
    PartInstance,
    PressureLoad,
    PrintSetup,
    Project,
    Study,
    Support,
    Transform,
)
from tessera.model.results import StudyResult
from tessera.model.solver import (
    AnalysisKind,
    FixedSupport,
    IsotropicMaterial,
    NodalLoads,
    Section,
    SolveJob,
    SolverAdapter,
    SolverError,
    SolveResult,
)

__all__ = [
    "SCHEMA_VERSION",
    "TET10_EDGES",
    "TRI6_EDGES",
    "AnalysisKind",
    "Assembly",
    "Connection",
    "ConnectionKind",
    "FaceInfo",
    "FaceSet",
    "FixedSupport",
    "ForceLoad",
    "InputError",
    "IsotropicMaterial",
    "Load",
    "Mesh",
    "NodalLoads",
    "PartInstance",
    "PressureLoad",
    "PrintSetup",
    "Project",
    "Section",
    "SolveJob",
    "SolveResult",
    "SolverAdapter",
    "SolverError",
    "Study",
    "StudyResult",
    "Support",
    "Transform",
    "TriangleSurface",
]
