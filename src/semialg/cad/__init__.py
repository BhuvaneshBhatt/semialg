from __future__ import annotations

import sys as _sys
import types as _types

from .bounds import (
    AlgebraicRootFunction,
    CADBound,
    CADCellBoundsCertificate,
    CertifiedRootComparison,
    DelineabilityCertificate,
    RootOrderCertificate,
    verify_cad_cell_bounds,
)
from .decomposition import (
    CADDiagnostics,
    CompleteCAD,
    decomp_collins_complete,
    decomp_from_proj_tower,
)
from .integration import (
    CADCellIntegral,
    IntrinsicCellStratum,
    IntrinsicStratification,
    full_dimensional_cell_integral,
    full_dimensional_solution_integrals,
    intrinsic_cell_integral,
    intrinsic_solution_integrals,
    stratify_intrinsic_solution,
)
from .lifting.sign_invariance import SignInvarianceCheck, verify_cad_sign_inv, verify_recorded_signs
from .lifting.stack import CADCell, sign_table
from .performance_cache import CADCacheStats, cad_cache_stats, clear_cad_caches
from .projection.collins import (
    ProjectionLevel,
    ProjectionPolynomial,
    ProjectionTower,
    build_collins_proj_set,
    collins_proj_entries,
    collins_projection_step,
)
from .projection.reduced import ProjectionValidity, ReducedProjectionTower

__all__ = [
    "CADCell",
    "CADBound",
    "AlgebraicRootFunction",
    "CertifiedRootComparison",
    "DelineabilityCertificate",
    "RootOrderCertificate",
    "CADCellBoundsCertificate",
    "verify_cad_cell_bounds",
    "CADCellIntegral",
    "IntrinsicCellStratum",
    "IntrinsicStratification",
    "stratify_intrinsic_solution",
    "full_dimensional_cell_integral",
    "full_dimensional_solution_integrals",
    "intrinsic_cell_integral",
    "intrinsic_solution_integrals",
    "CADDiagnostics",
    "CompleteCAD",
    "ProjectionLevel",
    "ProjectionPolynomial",
    "ProjectionTower",
    "ProjectionValidity",
    "ReducedProjectionTower",
    "SignInvarianceCheck",
    "build_collins_proj_set",
    "collins_projection_step",
    "collins_proj_entries",
    "decomp_collins_complete",
    "decomp_from_proj_tower",
    "sign_table",
    "verify_cad_sign_inv",
    "verify_recorded_signs",
    "CADCacheStats",
    "cad_cache_stats",
    "clear_cad_caches",
]


# Keep ``from semialg import cad`` callable even after Python has loaded the
# internal ``semialg.cad`` package and assigned it on the parent module.
# Tests and user code historically use both ``semialg.cad.cells`` and the
# top-level ``cad(...)`` helper, so the subpackage behaves as a thin callable
# proxy to the public decomposition API.
class _CallableCADModule(_types.ModuleType):
    def __call__(self, *args, **kwargs):
        from ..decomposition import cad as _public_cad

        return _public_cad(*args, **kwargs)


_sys.modules[__name__].__class__ = _CallableCADModule
