"""Symbolic function properties for SymPy expressions."""

from importlib.metadata import PackageNotFoundError
from importlib.metadata import version as _distribution_version

from .api import (
    analytic,
    bijective,
    boundary_behavior,
    boundedness,
    branch_cuts,
    branch_points,
    change_of_variables,
    complete_monotonicity,
    complex_zeros,
    continuous,
    convexity,
    determinant_nonvanishing,
    differentiability,
    discontinuities,
    domain,
    geometry,
    holder_regularity,
    homogeneity,
    injective,
    integrability,
    jacobian,
    level_boundedness,
    lipschitz,
    local_invertibility,
    meromorphic,
    mode_existence,
    monotonicity,
    oscillation,
    periodicity,
    poles,
    properties,
    range,
    sign,
    singularities,
    singularity_analysis,
    smooth,
    smoothness,
    surjective,
    symmetry,
    transformation_geometry,
    uniformly_continuous,
    variation,
    zeros,
)
from .assumption_closure import AssumptionClosure, close_assumptions
from .assumptions import (
    PropertyDecision,
    PropertyEnforcementError,
    PropertyKnowledge,
    PropertyProvenance,
    clear_entailment_cache,
    decide,
    entails,
    normalize_assumptions,
    require_decision,
)
from .cache_control import cache_info, clear_caches
from .certified_operations import (
    ArgumentPrincipleResult,
    ContinuationCertificate,
    ContourDeformationCertificate,
    EquationSolutionCertificate,
    HomotopyCertificate,
    InverseCertificate,
    LaurentCertificate,
    ResidueCertificate,
    SingularityClassificationCertificate,
    WindingCertificate,
    analytic_continue,
    argument_principle,
    classify_singularity,
    deform_contour,
    homotopic,
    inverse,
    laurent_series,
    monodromy,
    residue,
    solve,
    winding_number,
)
from .complex_certification import (
    ArgumentPrincipleCertificate,
    RectangleRoucheCertificate,
    RoucheCertificate,
    argument_principle_rectangle,
    rouche_polynomial_circle,
    rouche_polynomial_rectangle,
)
from .dlmf import DLMFFact, dlmf_audit_report, dlmf_facts
from .explanations import ResultExplanation, explain
from .geometry import (
    BoundaryBehavior,
    BoundaryLimit,
    FunctionGeometry,
    LevelBoundednessCertificate,
    ModeExistenceCertificate,
)
from .global_multivariate import GlobalMapCertificate, affine_global_mapping
from .infinite_root_families import affine_trig_root_family, structured_zero_geometry
from .integration import IntegralResult, ParametricContour, integrate
from .parameter_cells import minimize_strata, refine_parameter_cells, specialize_cells
from .parameter_discriminants import (
    ParameterTransitionLocus,
    polynomial_transition_locus,
)
from .proof_budget import ProofBudget, route_proofs
from .real_analysis.advanced_convexity import ConvexityCertificate
from .real_analysis.complete_monotonicity import CompleteMonotonicityCertificate
from .real_analysis.integrability import IntegrabilityResult
from .real_analysis.monotonic_certificates import (
    IndexedMonotonicPartition,
    InverseSolutionCertificate,
    MonotonicDecomposition,
    MonotonicInterval,
    MonotonicPartition,
    NonMonotonicityWitness,
    PartitionRange,
    certified_monotonic_inverse,
    clear_monotonic_partition_cache,
    injectivity_from_partition,
    monotonic_decomposition,
    monotonic_partition,
    range_from_monotonic_partition,
)
from .real_analysis.oscillation import OscillationCertificate
from .real_analysis.piecewise import (
    PiecewiseBoundary,
    PiecewiseBoundaryCertificate,
    PiecewiseBranch,
    PiecewiseDecomposition,
    certify_piecewise_boundaries,
    piecewise_decomposition,
)
from .real_analysis.regularity_extensions import HolderRegularityResult, LipschitzResult
from .real_analysis.smoothness import SmoothnessResult
from .real_analysis.structural_properties import (
    HomogeneityResult,
    PeriodicityResult,
    SymmetryResult,
)
from .real_analysis.variation import VariationResult
from .registry import (
    Family,
    FunctionRegistry,
    FunctionSlot,
    InstantiatedPropertyRule,
    PropertyRule,
    RuleCase,
    registry,
)
from .registry_derivatives import (
    ChainRuleCertificate,
    ChainRuleTerm,
    DerivativeRule,
    GradientCertificate,
    GradientComponent,
    HessianCertificate,
    IdentityRule,
    InstantiatedDerivativeRule,
    InstantiatedIdentityRule,
    JacobianCertificate,
    MultivariateChainRuleCertificate,
    registry_chain_rule,
    registry_derivative,
    registry_gradient,
    registry_hessian,
    registry_jacobian,
    rewrite_with_registry,
)
from .registry_transport import (
    AnalysisSubject,
    SubjectMap,
    TransportedPropertyRule,
    transport_instantiated_rule,
)
from .results import (
    Analyticity,
    Boundedness,
    ConditionalExpression,
    Continuity,
    Convexity,
    DetailedPropertyResult,
    Differentiability,
    Homogeneity,
    Meromorphicity,
    Monotonicity,
    Oscillation,
    Periodicity,
    PropertyResult,
    Sign,
    Singularity,
    SingularityKind,
    SingularityResult,
    Symmetry,
    TruthValue,
    as_detailed_result,
    explain_result,
)
from .singularities import SingularityAnalysis
from .special_theorems import GuardedTheorem, executable_special_theorems
from .statistical_special_functions import (
    MittagLefflerE,
    OwenT,
    RegularizedGammaP,
    RegularizedGammaQ,
)
from .theorem_closure import TheoremFact, close_theorems
from .transcendental import (
    AlgebraizationResult,
    AnalyticNeighborhoodCertificate,
    CertifiedComplexRootRectangle,
    CertifiedIntegerBracket,
    CertifiedRootInterval,
    ComplexRectangle,
    ComplexRootIsolationResult,
    ExponentialFourierCount,
    ExtremumEnvelope,
    ExtremumFamilyCertificate,
    ExtremumOrderingTheorem,
    ExtremumRange,
    IndexedExtremumValue,
    IndexedTailIndex,
    IndexedTailInversionCertificate,
    IndexedZero,
    IntervalPropertyCell,
    ParametricZeroFamily,
    PeriodicCondition,
    PeriodicSet,
    PropertyIntervalDecomposition,
    RootCertificate,
    RootIsolationResult,
    UnivariateReductionResult,
    ZeroGeometryCertificate,
    ZeroInterlacing,
    ZeroOrdering,
    algebraize,
    analytic_neighborhood,
    argument_principle_count,
    exponential_fourier_count,
    exponential_polynomial_terms,
    extremum_geometry,
    extremum_ordering_theorems,
    holomorphic_on_rectangle,
    indexed_tail_integer_bracket,
    indexed_tail_inversion,
    isolate_complex_roots,
    property_interval_decomposition,
    reduce_univariate_condition,
    robust_real_sign,
    simplify_extremum_range,
    symbolic_extremum_range,
    transcendental_roots,
    zero_geometry,
    zero_ordering_bounds,
)
from .transcendental import (
    ArgumentPrincipleCertificate as TranscendentalArgumentPrincipleCertificate,
)
from .transcendental.root_geometry import MonotoneRootCount, derivative_cell_root_count
from .transcendental.unified_real import (
    UnifiedRealAnalysis,
    analyze_transcendental_real,
)
from .transformation import (
    ChangeOfVariablesCertificate,
    DeterminantNonvanishingCertificate,
    LocalInvertibilityCertificate,
    TransformationGeometry,
)
from .verification import (
    verify_certificate,
    verify_convexity_certificate,
    verify_extremum_family,
    verify_extremum_ordering,
    verify_extremum_range,
    verify_indexed_tail_inversion,
    verify_integer_bracket,
    verify_monotonic_decomposition,
    verify_monotonic_partition,
    verify_root_certificate,
    verify_zero_geometry,
)

try:
    __version__ = _distribution_version("funcprops")
except PackageNotFoundError:
    __version__ = "0.3.0"

__all__ = [
    "Analyticity",
    "ArgumentPrincipleCertificate",
    "ArgumentPrincipleResult",
    "AssumptionClosure",
    "BoundaryBehavior",
    "BoundaryLimit",
    "Boundedness",
    "ChangeOfVariablesCertificate",
    "CompleteMonotonicityCertificate",
    "ConditionalExpression",
    "ContinuationCertificate",
    "Continuity",
    "ContourDeformationCertificate",
    "Convexity",
    "ConvexityCertificate",
    "DetailedPropertyResult",
    "DeterminantNonvanishingCertificate",
    "Differentiability",
    "EquationSolutionCertificate",
    "FunctionGeometry",
    "GlobalMapCertificate",
    "GuardedTheorem",
    "HolderRegularityResult",
    "Homogeneity",
    "HomogeneityResult",
    "HomotopyCertificate",
    "IntegrabilityResult",
    "IntegralResult",
    "InverseCertificate",
    "LaurentCertificate",
    "LevelBoundednessCertificate",
    "LipschitzResult",
    "LocalInvertibilityCertificate",
    "Meromorphicity",
    "ModeExistenceCertificate",
    "MonotoneRootCount",
    "Monotonicity",
    "Oscillation",
    "OscillationCertificate",
    "ParameterTransitionLocus",
    "ParametricContour",
    "Periodicity",
    "PeriodicityResult",
    "ProofBudget",
    "PropertyDecision",
    "PropertyKnowledge",
    "PropertyProvenance",
    "PropertyResult",
    "RectangleRoucheCertificate",
    "ResidueCertificate",
    "ResultExplanation",
    "RoucheCertificate",
    "Sign",
    "Singularity",
    "SingularityAnalysis",
    "SingularityClassificationCertificate",
    "SingularityKind",
    "SingularityResult",
    "SmoothnessResult",
    "Symmetry",
    "SymmetryResult",
    "TheoremFact",
    "TranscendentalArgumentPrincipleCertificate",
    "TransformationGeometry",
    "TruthValue",
    "UnifiedRealAnalysis",
    "VariationResult",
    "WindingCertificate",
    "affine_global_mapping",
    "affine_trig_root_family",
    "analytic",
    "analytic_continue",
    "analyze_transcendental_real",
    "argument_principle",
    "argument_principle_rectangle",
    "as_detailed_result",
    "bijective",
    "boundary_behavior",
    "boundedness",
    "branch_cuts",
    "branch_points",
    "cache_info",
    "change_of_variables",
    "classify_singularity",
    "clear_caches",
    "close_assumptions",
    "close_theorems",
    "complete_monotonicity",
    "complex_zeros",
    "continuous",
    "convexity",
    "decide",
    "deform_contour",
    "derivative_cell_root_count",
    "determinant_nonvanishing",
    "differentiability",
    "discontinuities",
    "domain",
    "entails",
    "executable_special_theorems",
    "explain",
    "explain_result",
    "geometry",
    "holder_regularity",
    "homogeneity",
    "homotopic",
    "injective",
    "integrability",
    "integrate",
    "inverse",
    "jacobian",
    "laurent_series",
    "level_boundedness",
    "lipschitz",
    "local_invertibility",
    "meromorphic",
    "minimize_strata",
    "mode_existence",
    "monodromy",
    "monotonicity",
    "oscillation",
    "periodicity",
    "poles",
    "polynomial_transition_locus",
    "properties",
    "range",
    "refine_parameter_cells",
    "registry",
    "residue",
    "rouche_polynomial_circle",
    "rouche_polynomial_rectangle",
    "route_proofs",
    "sign",
    "singularities",
    "singularity_analysis",
    "smooth",
    "smoothness",
    "solve",
    "specialize_cells",
    "structured_zero_geometry",
    "surjective",
    "symmetry",
    "transformation_geometry",
    "uniformly_continuous",
    "variation",
    "verify_certificate",
    "winding_number",
    "zeros",
]

# Consolidated proof/certification infrastructure.
from .asymptotic_tails import TailSignCertificate, certify_tail_sign
from .certified_balls import RealBallCertificate, certified_real_ball
from .certified_inequalities import InequalityCertificate, certify_inequality
from .certified_optimization import GlobalMinimumCertificate, certify_global_minimum
from .derivative_inequalities import (
    DerivativeChainCertificate,
    derivative_chain_positive,
)
from .diagnostics import Completeness, UnknownReason
from .global_analysis import GlobalInequalityCertificate, certify_global_positive
from .global_multivariate import ComposedMapCertificate, compose_global_maps
from .infinite_root_families import IndexedRootFamily, affine_sine_root_family
from .meromorphic_contours import (
    MeromorphicContourCertificate,
    meromorphic_rectangle_count,
)
from .multivariate_intervals import KrawczykBoxCertificate, krawczyk_box
from .parameter_cells import (
    exponential_affine_transition_cells,
    quadratic_discriminant_cells,
)
from .parameter_discriminants import (
    TranscendentalDiscriminant,
    exp_linear_discriminant,
    exp_linear_root_strata,
)
from .proof_explanations import ProofExplanation, explain_proof
from .proof_telemetry import ProofTelemetry, RouteObservation

# Compositional/certified analysis extensions.
from .property_composition import CompositionFact, compose_property
from .taylor_models import TaylorModel, taylor_model
from .theorem_registry import (
    THEOREMS,
    TheoremRecord,
    generated_adversarial_cases,
    validate_registry,
)
from .verified_roots import RealRootInterval, verify_real_root_interval
