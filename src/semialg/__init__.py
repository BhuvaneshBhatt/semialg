"""Semialgebraic reasoning, CAD, and real quantifier elimination for Python."""

from __future__ import annotations

__all__ = [
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
]

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
_ROOT_CLASSIFICATION_ATTRS = {"RootClassificationCell", "RootClassificationResult"}
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
_PARAMETRIC_INTEGRATION_ATTRS = {"ParametricIntegralResult"}
_IMPLICIT_UTILS_ATTRS = {"ImplicitFormulaPiece", "SymbolicBoxBounds", "VerticalBoundCell2D"}
_CONNECTIVITY_ATTRS = {
    "CADAdjacencyEdge",
    "CADConnectedComponent",
    "CADConnectivityGraph",
}

_CAD_CELL_ATTRS = {
    "CylindricalCoordinateConstraint",
    "CylindricalSolutionCell",
    "CylindricalSolution",
    "StructuredCADLevel",
    "StructuredCADCell",
    "StructuredCADCellDecomposition",
}

_OPTIMIZATION_ATTRS = {"FunctionRangeResult", "OptimizationResult"}
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
_PARAMETER_STRATIFICATION_ATTRS = {"ParameterStratum", "ParameterizedCylindricalDecomposition"}
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


def cad(*args, **kwargs):
    from .decomposition import cad as impl

    return impl(*args, **kwargs)


def cad_text(*args, **kwargs):
    from .decomposition import cad_text as impl

    return impl(*args, **kwargs)


def component_instances(*args, **kwargs):
    from .decomposition import component_instances as impl

    return impl(*args, **kwargs)


def component_instances_text(*args, **kwargs):
    from .decomposition import component_instances_text as impl

    return impl(*args, **kwargs)


def generic_cad(*args, **kwargs):
    from .decomposition import generic_cad as impl

    return impl(*args, **kwargs)


def generic_cad_text(*args, **kwargs):
    from .decomposition import generic_cad_text as impl

    return impl(*args, **kwargs)


def classify_real_roots(*args, **kwargs):
    from .root_classification import classify_real_roots as impl

    return impl(*args, **kwargs)


def sample_point(*args, **kwargs):
    from .sampling import sample_point as impl

    return impl(*args, **kwargs)


def sample_points(*args, **kwargs):
    from .sampling import sample_points as impl

    return impl(*args, **kwargs)


def sign_at(*args, **kwargs):
    from .sampling import sign_at as impl

    return impl(*args, **kwargs)


def sign_vector(*args, **kwargs):
    from .sampling import sign_vector as impl

    return impl(*args, **kwargs)


def semialgebraic_measure(*args, **kwargs):
    from .measure import semialgebraic_measure as impl

    return impl(*args, **kwargs)


def integrate_over_region(*args, **kwargs):
    from .region_integrate import integrate_over_region as impl

    return impl(*args, **kwargs)


def reduce_region_integral(*args, **kwargs):
    from .region_integrate import reduce_region_integral as impl

    return impl(*args, **kwargs)


def integrate_over_parametric_region(*args, **kwargs):
    from .parametric_integration import integrate_over_parametric_region as impl

    return impl(*args, **kwargs)


def reduce_parametric_integral(*args, **kwargs):
    from .parametric_integration import reduce_parametric_integral as impl

    return impl(*args, **kwargs)


def metric_jacobian_factor(*args, **kwargs):
    from .parametric_integration import metric_jacobian_factor as impl

    return impl(*args, **kwargs)


def integrate_over_standard_region(*args, **kwargs):
    from .standard_region_integrate import integrate_over_standard_region as impl

    return impl(*args, **kwargs)


def RegionUnion(*args, **kwargs):
    from .standard_regions import RegionUnion as impl

    return impl(*args, **kwargs)


def RegionIntersection(*args, **kwargs):
    from .standard_regions import RegionIntersection as impl

    return impl(*args, **kwargs)


def RegionDifference(*args, **kwargs):
    from .standard_regions import RegionDifference as impl

    return impl(*args, **kwargs)


def RegionSymmetricDifference(*args, **kwargs):
    from .standard_regions import RegionSymmetricDifference as impl

    return impl(*args, **kwargs)


def discretize_region_geometry(*args, **kwargs):
    from .solution_geometry import discretize_region_geometry as impl

    return impl(*args, **kwargs)


def plot_region_geometry(*args, **kwargs):
    from .solution_geometry import plot_region_geometry as impl

    return impl(*args, **kwargs)


def semialgebraic_level_function(*args, **kwargs):
    from .implicit_utils import semialgebraic_level_function as impl

    return impl(*args, **kwargs)


def decompose_implicit_formula(*args, **kwargs):
    from .implicit_utils import decompose_implicit_formula as impl

    return impl(*args, **kwargs)


def extract_symbolic_box_bounds(*args, **kwargs):
    from .implicit_utils import extract_symbolic_box_bounds as impl

    return impl(*args, **kwargs)


def decompose_cylindrical_formula_to_vertical_bounds_2d(*args, **kwargs):
    from .implicit_utils import decompose_cylindrical_formula_to_vertical_bounds_2d as impl

    return impl(*args, **kwargs)


def extract_structured_cad_cells(*args, **kwargs):
    from .cad.cells import extract_structured_cad_cells as impl

    return impl(*args, **kwargs)


def extract_explicit_cylindrical_solution(*args, **kwargs):
    from .cad.cells import extract_explicit_cylindrical_solution as impl

    return impl(*args, **kwargs)


def extract_cylindrical_solution(*args, **kwargs):
    from .cad.cells import extract_cylindrical_solution as impl

    return impl(*args, **kwargs)


def cylindrical_solution_from_structured(*args, **kwargs):
    from .cad.cells import cylindrical_solution_from_structured as impl

    return impl(*args, **kwargs)


def structured_cad_cells_to_vertical_bounds_2d(*args, **kwargs):
    from .cad.cells import structured_cad_cells_to_vertical_bounds_2d as impl

    return impl(*args, **kwargs)


def extract_vertical_bounds_from_cad_2d(*args, **kwargs):
    from .cad.cells import extract_vertical_bounds_from_cad_2d as impl

    return impl(*args, **kwargs)


def region_moment(*args, **kwargs):
    from .moments import region_moment as impl

    return impl(*args, **kwargs)


def region_centroid(*args, **kwargs):
    from .moments import region_centroid as impl

    return impl(*args, **kwargs)


def region_covariance(*args, **kwargs):
    from .moments import region_covariance as impl

    return impl(*args, **kwargs)


def function_range(*args, **kwargs):
    from .optimization import function_range as impl

    return impl(*args, **kwargs)


def semialgebraic_minimize(*args, **kwargs):
    from .optimization import semialgebraic_minimize as impl

    return impl(*args, **kwargs)


def semialgebraic_maximize(*args, **kwargs):
    from .optimization import semialgebraic_maximize as impl

    return impl(*args, **kwargs)


def subresultant_prs(*args, **kwargs):
    from .algebraic import subresultant_prs as impl

    return impl(*args, **kwargs)


def principal_subresultant_coefficients(*args, **kwargs):
    from .algebraic import principal_subresultant_coefficients as impl

    return impl(*args, **kwargs)


def compute_border_basis(*args, **kwargs):
    from .algebraic import compute_border_basis as impl

    return impl(*args, **kwargs)


def compute_border_basis_linear(*args, **kwargs):
    from .algebraic import compute_border_basis_linear as impl

    return impl(*args, **kwargs)


def is_satisfiable(*args, **kwargs):
    from .decision import is_satisfiable as impl

    return impl(*args, **kwargs)


def is_tautology(*args, **kwargs):
    from .decision import is_tautology as impl

    return impl(*args, **kwargs)


def implies(*args, **kwargs):
    from .decision import implies as impl

    return impl(*args, **kwargs)


def equivalent(lhs, rhs, variables=None, *args, **kwargs):
    if kwargs.get("return_result", False):
        from .decision import equivalent as impl

        return impl(lhs, rhs, variables, *args, **kwargs)
    import sympy as sp

    left = sp.sympify(lhs)
    right = sp.sympify(rhs)
    if sp.sstr(left) == sp.sstr(right):
        return True
    if variables is not None and len(tuple(variables)) == 1:
        try:
            return bool(left.as_set() == right.as_set())
        except (TypeError, ValueError, NotImplementedError, AttributeError):
            pass
    from .decision import equivalent as impl

    return impl(left, right, variables, *args, **kwargs)


def solve_semialgebraic(*args, **kwargs):
    from .decision import solve_semialgebraic as impl

    return impl(*args, **kwargs)


def canonicalize_one_dimensional_formula(*args, **kwargs):
    from .decision import canonicalize_one_dimensional_formula as impl

    return impl(*args, **kwargs)


def discretize_solution(*args, **kwargs):
    from .solution_geometry import discretize_solution as impl

    return impl(*args, **kwargs)


def plot_solution(*args, **kwargs):
    from .solution_geometry import plot_solution as impl

    return impl(*args, **kwargs)


def function_domain(*args, **kwargs):
    from .domain_solve import function_domain as impl

    return impl(*args, **kwargs)


def is_real_valued(*args, **kwargs):
    from .domain_solve import is_real_valued as impl

    return impl(*args, **kwargs)


def normalize_domain_sensitive_constraints(*args, **kwargs):
    from .domain_solve import normalize_domain_sensitive_constraints as impl

    return impl(*args, **kwargs)


def parameterized_cylindrical_decomposition(*args, **kwargs):
    from .parameter_stratification import parameterized_cylindrical_decomposition as impl

    return impl(*args, **kwargs)


def solvability_conditions(*args, **kwargs):
    from .parameters import solvability_conditions as impl

    return impl(*args, **kwargs)


def root_count_conditions(*args, **kwargs):
    from .parameters import root_count_conditions as impl

    return impl(*args, **kwargs)


def region_union(*args, **kwargs):
    from .regions.operations import region_union as impl

    return impl(*args, **kwargs)


def region_intersection(*args, **kwargs):
    from .regions.operations import region_intersection as impl

    return impl(*args, **kwargs)


def region_difference(*args, **kwargs):
    from .regions.operations import region_difference as impl

    return impl(*args, **kwargs)


def region_complement(*args, **kwargs):
    from .regions.operations import region_complement as impl

    return impl(*args, **kwargs)


def region_closure(*args, **kwargs):
    from .regions.operations import region_closure as impl

    return impl(*args, **kwargs)


def region_interior(*args, **kwargs):
    from .regions.operations import region_interior as impl

    return impl(*args, **kwargs)


def region_boundary(*args, **kwargs):
    from .regions.operations import region_boundary as impl

    return impl(*args, **kwargs)


def region_dimension(*args, **kwargs):
    from .regions.operations import region_dimension as impl

    return impl(*args, **kwargs)


def region_components(*args, **kwargs):
    from .regions.operations import region_components as impl

    return impl(*args, **kwargs)


def simplify_boole(*args, **kwargs):
    from .symbolic_simplify import simplify_boole as impl

    return impl(*args, **kwargs)


def simplify_piecewise(*args, **kwargs):
    from .symbolic_simplify import simplify_piecewise as impl

    return impl(*args, **kwargs)


def simplify_system(*args, **kwargs):
    from .reasoning import simplify_system as impl

    return impl(*args, **kwargs)


def prove_positive(*args, **kwargs):
    from .reasoning import prove_positive as impl

    return impl(*args, **kwargs)


def prove_nonnegative(*args, **kwargs):
    from .reasoning import prove_nonnegative as impl

    return impl(*args, **kwargs)


def prove_negative(*args, **kwargs):
    from .reasoning import prove_negative as impl

    return impl(*args, **kwargs)


def prove_nonpositive(*args, **kwargs):
    from .reasoning import prove_nonpositive as impl

    return impl(*args, **kwargs)


def region_subset(*args, **kwargs):
    from .reasoning import region_subset as impl

    return impl(*args, **kwargs)


def region_equal(*args, **kwargs):
    from .reasoning import region_equal as impl

    return impl(*args, **kwargs)


def region_disjoint(*args, **kwargs):
    from .reasoning import region_disjoint as impl

    return impl(*args, **kwargs)


def region_bounded(*args, **kwargs):
    from .reasoning import region_bounded as impl

    return impl(*args, **kwargs)


def region_closed(*args, **kwargs):
    from .reasoning import region_closed as impl

    return impl(*args, **kwargs)


def region_compact(*args, **kwargs):
    from .reasoning import region_compact as impl

    return impl(*args, **kwargs)


def simplify_under_assumptions(*args, **kwargs):
    from .reasoning import simplify_under_assumptions as impl

    return impl(*args, **kwargs)


def qe_by_complete_cad(*args, **kwargs):
    from .qe import qe_by_complete_cad as impl

    return impl(*args, **kwargs)


def root_of(*args, **kwargs):
    from .reconstruct import root_of as impl

    return impl(*args, **kwargs)


def semialgebraicize(*args, **kwargs):
    from .preprocess import semialgebraicize as impl

    return impl(*args, **kwargs)


def find_instance(*args, **kwargs):
    from .solve import find_instance as impl

    return impl(*args, **kwargs)


def find_instance_formula(*args, **kwargs):
    from .solve import find_instance_formula as impl

    return impl(*args, **kwargs)


def find_instance_text(*args, **kwargs):
    from .solve import find_instance_text as impl

    return impl(*args, **kwargs)


def reduce_formula(*args, **kwargs):
    from .solve import reduce_formula as impl

    return impl(*args, **kwargs)


def reduce_text(*args, **kwargs):
    from .solve import reduce_text as impl

    return impl(*args, **kwargs)


def build_cad_adjacency_graph(*args, **kwargs):
    from .connectivity import build_cad_adjacency_graph as impl

    return impl(*args, **kwargs)


def extract_cad_connectivity(*args, **kwargs):
    from .connectivity import extract_cad_connectivity as impl

    return impl(*args, **kwargs)


def resolve_formula(*args, **kwargs):
    from .solve import resolve_formula as impl

    return impl(*args, **kwargs)


def resolve_text(*args, **kwargs):
    from .solve import resolve_text as impl

    return impl(*args, **kwargs)


def __getattr__(name: str):
    if name in _DECOMP_ATTRS:
        from . import decomposition as mod
    elif name in _PREPROCESS_ATTRS:
        from . import preprocess as mod
    elif name in _QE_ATTRS:
        from . import qe as mod
    elif name in _RECONSTRUCT_ATTRS:
        from . import reconstruct as mod
    elif name in _GENERIC_ATTRS:
        from . import generic as mod
    elif name in _ROOT_CLASSIFICATION_ATTRS:
        from . import root_classification as mod
    elif name in _ALGEBRAIC_ATTRS:
        from . import algebraic as mod
    elif name in _MEASURE_ATTRS:
        from . import measure as mod
    elif name in _REGION_INTEGRATE_ATTRS:
        from . import region_integrate as mod
    elif name in _STANDARD_REGION_ATTRS:
        from . import standard_regions as mod
    elif name in _PARAMETRIC_INTEGRATION_ATTRS:
        from . import parametric_integration as mod
    elif name in _IMPLICIT_UTILS_ATTRS:
        from . import implicit_utils as mod
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
    elif name in _PARAMETER_STRATIFICATION_ATTRS:
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
