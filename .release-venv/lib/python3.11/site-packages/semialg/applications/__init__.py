"""Domain-oriented applications built on semialg's certified core algorithms."""

from .barrier_certificates import BarrierVerificationResult, verify_barrier_certificate
from .constraint_analysis import (
    ConstraintRedundancyResult,
    FeasibleSetDiagnosticResult,
    analyze_constraint_redundancy,
    diagnose_feasible_set,
)
from .control import (
    PolynomialStabilityResult,
    polynomial_stability_analysis,
    polynomial_stability_region,
)
from .lyapunov import LyapunovVerificationResult, verify_lyapunov_function
from .model_comparison import PolynomialModelComparisonResult, compare_polynomial_models
from .optimization_benchmarks import (
    NumericOptimizationCheck,
    OptimizationBenchmark,
    exact_optimization_benchmark,
    validate_numeric_optimization,
)
from .parameter_regimes import (
    ParameterRegimeResult,
    analyze_parameter_regimes,
    analyze_root_count_regimes,
)
from .probability import PolynomialProbabilityResult, geometric_probability, polynomial_probability
from .response_surfaces import ResponseSurfaceResult, analyze_response_surface
from .robust_design import RobustParameterResult, robust_parameter_analysis, robust_parameter_region
from .safety import InvariantVerificationResult, verify_polynomial_invariant
from .sensitivity import (
    SensitivityAnalysisResult,
    SensitivityDirectionResult,
    analyze_polynomial_sensitivity,
)
from .validation import (
    ValidationResult,
    validate_formula_equivalence,
    validate_identity,
    validate_range,
)

__all__ = [
    "BarrierVerificationResult",
    "ConstraintRedundancyResult",
    "FeasibleSetDiagnosticResult",
    "InvariantVerificationResult",
    "LyapunovVerificationResult",
    "NumericOptimizationCheck",
    "OptimizationBenchmark",
    "PolynomialModelComparisonResult",
    "PolynomialProbabilityResult",
    "PolynomialStabilityResult",
    "ParameterRegimeResult",
    "ResponseSurfaceResult",
    "RobustParameterResult",
    "SensitivityAnalysisResult",
    "SensitivityDirectionResult",
    "ValidationResult",
    "analyze_constraint_redundancy",
    "analyze_parameter_regimes",
    "analyze_polynomial_sensitivity",
    "analyze_response_surface",
    "analyze_root_count_regimes",
    "compare_polynomial_models",
    "diagnose_feasible_set",
    "geometric_probability",
    "exact_optimization_benchmark",
    "polynomial_probability",
    "polynomial_stability_analysis",
    "polynomial_stability_region",
    "robust_parameter_analysis",
    "robust_parameter_region",
    "validate_formula_equivalence",
    "validate_identity",
    "validate_numeric_optimization",
    "validate_range",
    "verify_barrier_certificate",
    "verify_lyapunov_function",
    "verify_polynomial_invariant",
]
