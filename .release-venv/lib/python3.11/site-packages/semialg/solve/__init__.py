from __future__ import annotations

from .complete import (
    CompleteSolveResult,
    reduce_complete_expr,
    reduce_complete_formula,
    reduce_complete_text,
    resolve_complete_text,
)
from .domains import (
    SemialgOptions,
    SolveDomain,
    SolveRequest,
    UnsupportedDomainError,
    apply_assumptions,
    normalize_assumptions,
    normalize_domain,
)
from .find_instance import InstanceResult, find_instance, find_instance_formula, find_instance_text
from .preprocess import PreprocessResult, semialgebraicize
from .reduce import reduce_formula, reduce_text
from .resolve import resolve_formula, resolve_text
from .result import SolveResult
from .zero_dimensional import (
    ZeroDimensionalSolveResult,
    is_zero_dimensional,
    solve_zero_dimensional_system,
)

__all__ = [
    "CompleteSolveResult",
    "PreprocessResult",
    "SemialgOptions",
    "SolveDomain",
    "SolveRequest",
    "UnsupportedDomainError",
    "SolveResult",
    "ZeroDimensionalSolveResult",
    "InstanceResult",
    "find_instance",
    "is_zero_dimensional",
    "solve_zero_dimensional_system",
    "find_instance_formula",
    "find_instance_text",
    "apply_assumptions",
    "normalize_assumptions",
    "normalize_domain",
    "semialgebraicize",
    "reduce_formula",
    "reduce_text",
    "reduce_complete_formula",
    "reduce_complete_expr",
    "reduce_complete_text",
    "resolve_formula",
    "resolve_text",
    "resolve_complete_text",
]
