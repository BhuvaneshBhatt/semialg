"""Public imports for exact reasoning and assumption simplification."""

from .reasoning_assumptions import AssumptionSimplificationResult, simplify_under_assumptions
from .reasoning_regions import (
    region_bounded,
    region_closed,
    region_compact,
    region_disjoint,
    region_equal,
    region_subset,
)
from .reasoning_signs import (
    SignClassificationResult,
    SignProofResult,
    function_sign,
    prove_negative,
    prove_nonnegative,
    prove_nonpositive,
    prove_nonzero,
    prove_positive,
    prove_zero,
)
from .reasoning_system import SimplifiedSystem, simplify_system

__all__ = [
    "SimplifiedSystem",
    "AssumptionSimplificationResult",
    "SignClassificationResult",
    "SignProofResult",
    "simplify_system",
    "prove_positive",
    "prove_nonnegative",
    "prove_negative",
    "prove_nonpositive",
    "prove_zero",
    "prove_nonzero",
    "function_sign",
    "region_subset",
    "region_equal",
    "region_disjoint",
    "region_bounded",
    "region_closed",
    "region_compact",
    "simplify_under_assumptions",
]
