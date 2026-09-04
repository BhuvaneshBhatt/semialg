"""Public imports for structured and cylindrical CAD-cell APIs."""

from .bounds import (
    AlgebraicRootFunction,
    CADBound,
    CADCellBoundsCertificate,
    CertifiedRootComparison,
    DelineabilityCertificate,
    RootOrderCertificate,
    verify_cad_cell_bounds,
)
from .cylindrical_solution import (
    CylindricalCoordinateConstraint,
    CylindricalDecompositionCertificate,
    CylindricalSolution,
    CylindricalSolutionCell,
    cylindrical_solution_from_structured,
    extract_cylindrical_solution,
)
from .explicit_cells import (
    extract_explicit_cylindrical_solution,
    extract_vertical_bounds_from_cad_2d,
    structured_cad_cells_to_vertical_bounds_2d,
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
from .structured_cells import (
    StructuredCADCell,
    StructuredCADCellDecomposition,
    StructuredCADLevel,
    extract_structured_cad_cells,
)

__all__ = [
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
    "CylindricalCoordinateConstraint",
    "CylindricalSolutionCell",
    "CylindricalSolution",
    "CylindricalDecompositionCertificate",
    "cylindrical_solution_from_structured",
    "extract_cylindrical_solution",
    "extract_explicit_cylindrical_solution",
    "StructuredCADLevel",
    "StructuredCADCell",
    "StructuredCADCellDecomposition",
    "extract_structured_cad_cells",
    "structured_cad_cells_to_vertical_bounds_2d",
    "extract_vertical_bounds_from_cad_2d",
]
