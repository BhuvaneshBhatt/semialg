from __future__ import annotations

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
    decomp_groebner_variety,
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
from .performance_cache import (
    CADCacheLimits,
    CADCacheStats,
    cad_cache_stats,
    clear_cad_caches,
    configure_cad_cache_limits,
)
from .projection.collins import (
    ProjectionLevel,
    ProjectionPolynomial,
    ProjectionTower,
    build_collins_proj_set,
    collins_proj_entries,
    collins_projection_step,
)
from .projection.groebner import (
    GroebnerVarietyProjection,
    build_groebner_variety_projection,
)
from .projection.reduced import ProjectionValidity, ReducedProjectionTower
from .selected_samples import SelectedCADSample, extract_selected_cad_samples

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
    "GroebnerVarietyProjection",
    "ReducedProjectionTower",
    "SignInvarianceCheck",
    "build_collins_proj_set",
    "build_groebner_variety_projection",
    "collins_projection_step",
    "collins_proj_entries",
    "decomp_collins_complete",
    "decomp_groebner_variety",
    "decomp_from_proj_tower",
    "sign_table",
    "verify_cad_sign_inv",
    "verify_recorded_signs",
    "CADCacheLimits",
    "CADCacheStats",
    "cad_cache_stats",
    "clear_cad_caches",
    "configure_cad_cache_limits",
    "SelectedCADSample",
    "extract_selected_cad_samples",
]
