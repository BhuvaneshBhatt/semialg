"""Semialgebraic reasoning, CAD, and real quantifier elimination for Python."""

from __future__ import annotations

from importlib import import_module

__version__ = "0.2.0b1"

__all__ = [
    "__version__",
    "CADComponent",
    "CADFunction",
    "CADOutput",
    "CADTreeNode",
    "CADOptions",
    "CADResult",
    "CompleteQEResult",
    "ComponentResult",
    "GenericCADResult",
    "GenericCADFunction",
    "GenericSplit",
    "InstanceResult",
    "MeasureResult",
    "RegionIntegralResult",
    "StandardRegion",
    "PointRegion",
    "IntervalRegion",
    "BoxRegion",
    "SimplexRegion",
    "PolygonRegion",
    "TetrahedronRegion",
    "PolyhedronRegion",
    "ParallelogramRegion",
    "ParallelepipedRegion",
    "PrismRegion",
    "PyramidRegion",
    "BallRegion",
    "SphereRegion",
    "SphericalShellRegion",
    "CylinderRegion",
    "ConeRegion",
    "StadiumRegion",
    "CapsuleRegion",
    "ParametricRegion",
    "TransformedRegion",
    "BooleanRegion",
    "RegionUnion",
    "RegionIntersection",
    "RegionDifference",
    "RegionSymmetricDifference",
    "ParametricIntegralResult",
    "metric_jacobian_factor",
    "reduce_parametric_integral",
    "integrate_over_parametric_region",
    "integrate_over_standard_region",
    "RegionIntegralPiece",
    "ReducedRegionIntegral",
    "ImplicitFormulaPiece",
    "SymbolicBoxBounds",
    "VerticalBoundCell2D",
    "CADAdjacencyEdge",
    "CADConnectedComponent",
    "CADConnectivityGraph",
    "build_cad_adjacency_graph",
    "extract_cad_connectivity",
    "CylindricalCoordinateConstraint",
    "CylindricalSolutionCell",
    "CylindricalSolution",
    "CylindricalDecompositionCertificate",
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
    "cylindrical_solution_from_structured",
    "extract_cylindrical_solution",
    "extract_explicit_cylindrical_solution",
    "StructuredCADLevel",
    "StructuredCADCell",
    "StructuredCADCellDecomposition",
    "extract_structured_cad_cells",
    "structured_cad_cells_to_vertical_bounds_2d",
    "extract_vertical_bounds_from_cad_2d",
    "semialgebraic_level_function",
    "decompose_implicit_formula",
    "extract_symbolic_box_bounds",
    "decompose_cylindrical_formula_to_vertical_bounds_2d",
    "FunctionRangeResult",
    "OptimizationResult",
    "ParametricOptimizationResult",
    "ParametricFunctionRangeResult",
    "OptimizationCertificationPolicy",
    "polynomial_locus_dimension",
    "PowerPolicy",
    "PreprocessResult",
    "RootFunction",
    "SemialgOptions",
    "SolveDomain",
    "UnsupportedDomainError",
    "cad",
    "cad_text",
    "component_instances",
    "component_instances_text",
    "find_instance",
    "find_instance_formula",
    "find_instance_text",
    "generic_cad",
    "generic_cad_text",
    "qe_by_complete_cad",
    "reduce_formula",
    "reduce_text",
    "resolve_formula",
    "resolve_text",
    "root_of",
    "RootClassificationCell",
    "RootClassificationResult",
    "classify_real_roots",
    "sample_point",
    "sample_points",
    "semialgebraic_measure",
    "integrate_over_region",
    "reduce_region_integral",
    "function_range",
    "semialgebraic_minimize",
    "semialgebraic_maximize",
    "sign_at",
    "sign_vector",
    "algebraic_cache_stats",
    "clear_algebraic_caches",
    "ExactComputationContext",
    "computation_context",
    "current_computation_context",
    "IntervalComponent",
    "SemialgebraicSolution",
    "EquivalenceResult",
    "ImplicationResult",
    "TautologyResult",
    "SatisfiabilityResult",
    "SolutionPlotData",
    "is_satisfiable",
    "is_tautology",
    "implies",
    "equivalent",
    "solve_semialgebraic",
    "discretize_solution",
    "plot_solution",
    "discretize_region_geometry",
    "plot_region_geometry",
    "canonicalize_one_dimensional_formula",
    "normalize_domain_sensitive_constraints",
    "is_real_valued",
    "function_domain",
    "DomainNormalizationResult",
    "SolvabilityConditionsResult",
    "RootCountConditionsResult",
    "solvability_conditions",
    "root_count_conditions",
    "region_union",
    "region_intersection",
    "region_difference",
    "region_complement",
    "region_closure",
    "region_interior",
    "region_boundary",
    "region_dimension",
    "region_components",
    "BooleanSimplificationResult",
    "PiecewiseSimplificationResult",
    "simplify_boole",
    "simplify_piecewise",
    "RegionMomentResult",
    "RegionCentroidResult",
    "RegionCovarianceResult",
    "region_moment",
    "region_centroid",
    "region_covariance",
    "SimplifiedSystem",
    "AssumptionSimplificationResult",
    "SignProofResult",
    "simplify_system",
    "prove_positive",
    "prove_nonnegative",
    "prove_negative",
    "prove_nonpositive",
    "region_subset",
    "region_equal",
    "region_disjoint",
    "region_bounded",
    "region_closed",
    "region_compact",
    "simplify_under_assumptions",
    "ConditionalBranch",
    "ParameterStratifiedResult",
    "ParameterStratificationCertificate",
    "conditional_result",
    "verify_parameter_stratification",
    "ParameterStratum",
    "ParameterizedCylindricalDecomposition",
    "parameterized_cylindrical_decomposition",
    "semialgebraicize",
    "ZeroDimensionalSolveResult",
    "is_zero_dimensional",
    "solve_zero_dimensional_system",
    "SubresultantPRSResult",
    "subresultant_prs",
    "principal_subresultant_coefficients",
    "BorderBasisError",
    "BorderBasisResult",
    "compute_border_basis",
    "compute_border_basis_linear",
    "DimensionMismatchError",
    "Exists",
    "ForAll",
    "apply_quantifiers",
    "split_quantifiers",
    "require_point_dimension",
    "require_same_length",
    "zip_equal",
]

_LAZY_EXPORTS = {
    "Exists": ".quantifiers",
    "ForAll": ".quantifiers",
    "apply_quantifiers": ".quantifiers",
    "split_quantifiers": ".quantifiers",
    "DimensionMismatchError": ".errors",
    "require_point_dimension": ".dimension_validation",
    "require_same_length": ".dimension_validation",
    "zip_equal": ".dimension_validation",
    "RegionDifference": ".standard_regions",
    "RegionIntersection": ".standard_regions",
    "RegionSymmetricDifference": ".standard_regions",
    "RegionUnion": ".standard_regions",
    "algebraic_cache_stats": ".algebraic",
    "build_cad_adjacency_graph": ".connectivity",
    "cad": ".decomposition",
    "cad_text": ".decomposition",
    "canonicalize_one_dimensional_formula": ".decision",
    "classify_real_roots": ".root_classification",
    "clear_algebraic_caches": ".algebraic",
    "component_instances": ".decomposition",
    "component_instances_text": ".decomposition",
    "computation_context": ".context",
    "compute_border_basis": ".algebraic",
    "compute_border_basis_linear": ".algebraic",
    "conditional_result": ".conditional",
    "current_computation_context": ".context",
    "cylindrical_solution_from_structured": ".cad.cells",
    "decompose_cylindrical_formula_to_vertical_bounds_2d": ".implicit_geometry",
    "decompose_implicit_formula": ".implicit_geometry",
    "discretize_region_geometry": ".solution_geometry",
    "discretize_solution": ".solution_geometry",
    "equivalent": ".decision",
    "extract_cad_connectivity": ".connectivity",
    "extract_cylindrical_solution": ".cad.cells",
    "extract_explicit_cylindrical_solution": ".cad.cells",
    "extract_structured_cad_cells": ".cad.cells",
    "extract_symbolic_box_bounds": ".implicit_geometry",
    "extract_vertical_bounds_from_cad_2d": ".cad.cells",
    "find_instance": ".solve",
    "find_instance_formula": ".solve",
    "find_instance_text": ".solve",
    "function_domain": ".domain_solve",
    "function_range": ".optimization",
    "generic_cad": ".decomposition",
    "generic_cad_text": ".decomposition",
    "implies": ".decision",
    "integrate_over_parametric_region": ".parametric_integration",
    "integrate_over_region": ".region_integrate",
    "integrate_over_standard_region": ".standard_region_integrate",
    "is_real_valued": ".domain_solve",
    "is_satisfiable": ".decision",
    "is_tautology": ".decision",
    "metric_jacobian_factor": ".parametric_integration",
    "normalize_domain_sensitive_constraints": ".domain_solve",
    "parameterized_cylindrical_decomposition": ".parameter_stratification",
    "plot_region_geometry": ".solution_geometry",
    "plot_solution": ".solution_geometry",
    "principal_subresultant_coefficients": ".algebraic",
    "prove_negative": ".reasoning",
    "prove_nonnegative": ".reasoning",
    "prove_nonpositive": ".reasoning",
    "prove_positive": ".reasoning",
    "qe_by_complete_cad": ".qe",
    "reduce_formula": ".solve",
    "reduce_parametric_integral": ".parametric_integration",
    "reduce_region_integral": ".region_integrate",
    "reduce_text": ".solve",
    "region_boundary": ".regions.operations",
    "region_bounded": ".reasoning",
    "region_centroid": ".moments",
    "region_closed": ".reasoning",
    "region_closure": ".regions.operations",
    "region_compact": ".reasoning",
    "region_complement": ".regions.operations",
    "region_components": ".regions.operations",
    "region_covariance": ".moments",
    "region_difference": ".regions.operations",
    "region_dimension": ".regions.operations",
    "region_disjoint": ".reasoning",
    "region_equal": ".reasoning",
    "region_interior": ".regions.operations",
    "region_intersection": ".regions.operations",
    "region_moment": ".moments",
    "region_subset": ".reasoning",
    "region_union": ".regions.operations",
    "resolve_formula": ".solve",
    "resolve_text": ".solve",
    "root_count_conditions": ".parameters",
    "root_of": ".reconstruct",
    "sample_point": ".sampling",
    "sample_points": ".sampling",
    "semialgebraic_level_function": ".implicit_geometry",
    "semialgebraic_maximize": ".optimization",
    "semialgebraic_measure": ".measure",
    "semialgebraic_minimize": ".optimization",
    "semialgebraicize": ".preprocess",
    "sign_at": ".sampling",
    "sign_vector": ".sampling",
    "simplify_boole": ".symbolic_simplify",
    "simplify_piecewise": ".symbolic_simplify",
    "simplify_system": ".reasoning",
    "simplify_under_assumptions": ".reasoning",
    "solvability_conditions": ".parameters",
    "solve_semialgebraic": ".decision",
    "structured_cad_cells_to_vertical_bounds_2d": ".cad.cells",
    "subresultant_prs": ".algebraic",
    "verify_parameter_stratification": ".conditional",
}

_DECOMP_ATTRS = {
    "CADComponent",
    "CADFunction",
    "CADOutput",
    "CADTreeNode",
    "CADOptions",
    "CADResult",
    "ComponentResult",
    "GenericCADResult",
    "GenericCADFunction",
}
_PREPROCESS_ATTRS = {"PowerPolicy", "PreprocessResult"}
_QE_ATTRS = {"CompleteQEResult"}
_RECONSTRUCT_ATTRS = {"RootFunction"}
_GENERIC_ATTRS = {"GenericSplit"}
_ROOT_CLASS_ATTRS = {"RootClassificationCell", "RootClassificationResult"}
_ALGEBRAIC_ATTRS = {"SubresultantPRSResult", "BorderBasisError", "BorderBasisResult"}
_MEASURE_ATTRS = {"MeasureResult"}
_REGION_INTEGRATE_ATTRS = {"RegionIntegralResult", "RegionIntegralPiece", "ReducedRegionIntegral"}
_STANDARD_REGION_ATTRS = {
    "StandardRegion",
    "PointRegion",
    "IntervalRegion",
    "BoxRegion",
    "SimplexRegion",
    "PolygonRegion",
    "TetrahedronRegion",
    "PolyhedronRegion",
    "ParallelogramRegion",
    "ParallelepipedRegion",
    "PrismRegion",
    "PyramidRegion",
    "BallRegion",
    "SphereRegion",
    "SphericalShellRegion",
    "CylinderRegion",
    "ConeRegion",
    "StadiumRegion",
    "CapsuleRegion",
    "ParametricRegion",
    "TransformedRegion",
    "BooleanRegion",
}
_PARAM_INTEGRAL_ATTRS = {"ParametricIntegralResult"}
_IMPLICIT_UTILS_ATTRS = {"ImplicitFormulaPiece", "SymbolicBoxBounds", "VerticalBoundCell2D"}
_CONNECTIVITY_ATTRS = {
    "CADAdjacencyEdge",
    "CADConnectedComponent",
    "CADConnectivityGraph",
}

_CAD_CELL_ATTRS = {
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
    "StructuredCADLevel",
    "StructuredCADCell",
    "StructuredCADCellDecomposition",
}

_CONTEXT_ATTRS = {"ExactComputationContext"}
_OPTIMIZATION_ATTRS = {
    "FunctionRangeResult",
    "OptimizationResult",
    "ParametricOptimizationResult",
    "ParametricFunctionRangeResult",
    "OptimizationCertificationPolicy",
    "polynomial_locus_dimension",
}
_SAMPLING_ATTRS = set()
_DECISION_ATTRS = {
    "IntervalComponent",
    "SemialgebraicSolution",
    "SatisfiabilityResult",
    "TautologyResult",
    "ImplicationResult",
    "EquivalenceResult",
}
_SOLUTION_GEOMETRY_ATTRS = {"SolutionPlotData"}
_DOMAIN_SOLVE_ATTRS = {"DomainNormalizationResult"}
_PARAMETER_ATTRS = {"SolvabilityConditionsResult", "RootCountConditionsResult"}
_CONDITIONAL_ATTRS = {
    "ConditionalBranch",
    "ParameterStratifiedResult",
    "ParameterStratificationCertificate",
}
_PARAM_STRATA_ATTRS = {"ParameterStratum", "ParameterizedCylindricalDecomposition"}
_MOMENT_ATTRS = {"RegionMomentResult", "RegionCentroidResult", "RegionCovarianceResult"}
_REASONING_ATTRS = {"SimplifiedSystem", "AssumptionSimplificationResult", "SignProofResult"}
_SYMBOLIC_SIMPLIFY_ATTRS = {"BooleanSimplificationResult", "PiecewiseSimplificationResult"}

_SOLVE_ATTRS = {
    "InstanceResult",
    "ZeroDimensionalSolveResult",
    "SemialgOptions",
    "SolveDomain",
    "UnsupportedDomainError",
    "is_zero_dimensional",
    "solve_zero_dimensional_system",
}


def __getattr__(name: str):
    module_name = _LAZY_EXPORTS.get(name)
    if module_name is not None:
        mod = import_module(module_name, __name__)
        value = getattr(mod, name)
        globals()[name] = value
        return value
    if name in _CONTEXT_ATTRS:
        from . import context as mod
    elif name in _DECOMP_ATTRS:
        from . import decomposition as mod
    elif name in _PREPROCESS_ATTRS:
        from . import preprocess as mod
    elif name in _QE_ATTRS:
        from . import qe as mod
    elif name in _RECONSTRUCT_ATTRS:
        from . import reconstruct as mod
    elif name in _GENERIC_ATTRS:
        from . import generic as mod
    elif name in _ROOT_CLASS_ATTRS:
        from . import root_classification as mod
    elif name in _ALGEBRAIC_ATTRS:
        from . import algebraic as mod
    elif name in _MEASURE_ATTRS:
        from . import measure as mod
    elif name in _REGION_INTEGRATE_ATTRS:
        from . import region_integrate as mod
    elif name in _STANDARD_REGION_ATTRS:
        from . import standard_regions as mod
    elif name in _PARAM_INTEGRAL_ATTRS:
        from . import parametric_integration as mod
    elif name in _IMPLICIT_UTILS_ATTRS:
        from . import implicit_geometry as mod
    elif name in _CONNECTIVITY_ATTRS:
        from . import connectivity as mod
    elif name in _CAD_CELL_ATTRS:
        from .cad import cells as mod
    elif name in _DECISION_ATTRS:
        from . import decision as mod
    elif name in _SOLUTION_GEOMETRY_ATTRS:
        from . import solution_geometry as mod
    elif name in _DOMAIN_SOLVE_ATTRS:
        from . import domain_solve as mod
    elif name in _OPTIMIZATION_ATTRS:
        from . import optimization as mod
    elif name in _PARAMETER_ATTRS:
        from . import parameters as mod
    elif name in _CONDITIONAL_ATTRS:
        from . import conditional as mod
    elif name in _PARAM_STRATA_ATTRS:
        from . import parameter_stratification as mod
    elif name in _MOMENT_ATTRS:
        from . import moments as mod
    elif name in _SOLVE_ATTRS:
        from . import solve as mod
    elif name in _REASONING_ATTRS:
        from . import reasoning as mod
    elif name in _SYMBOLIC_SIMPLIFY_ATTRS:
        from . import symbolic_simplify as mod
    else:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    value = getattr(mod, name)
    globals()[name] = value
    return value
