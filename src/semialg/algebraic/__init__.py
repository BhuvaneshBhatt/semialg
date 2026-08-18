"""Exact algebraic samples and sign utilities used by the CAD core."""

from .border_basis import (
    BorderBasisDiagnostics,
    BorderBasisError,
    BorderBasisResult,
    compute_border_basis,
    compute_border_basis_linear,
)
from .cache import (
    CACHE,
    CACHE_VERSION,
    RootIsolationCache,
    RootIsolationStats,
    root_isolation_costs,
)
from .comparison import compare_samples, sort_samples
from .intervals import RationalInterval
from .rational_univariate import (
    FilteredRationalUnivariateSolutions,
    RationalUnivariateError,
    RationalUnivariateFormulaResult,
    RationalUnivariatePoint,
    RationalUnivariateRepresentation,
    compute_rational_univariate_representation,
    evaluate_boolean_formula_at_point,
    evaluate_relation_at_point,
    filter_rur_solutions_by_constraints,
    sign_of_algebraic_expression,
    solve_and_filter_zero_dimensional_system_with_rur,
    solve_formula_with_rur,
    solve_rur_points,
    solve_rur_representation,
    solve_rur_semialgebraic_system,
    solve_zero_dimensional_system_with_rur,
)
from .roots import isolate_real_roots, refine_isol_intv, root_multiplicity
from .sample_points import choose_sector_sample
from .samples import AlgebraicRoot, RationalSample, Sample, sample_to_expr
from .signs import sign_at_sample
from .subresultants import (
    SubresultantPRSResult,
    principal_subresultant_coefficients,
    subresultant_prs,
)

compare_alg_numbers = compare_samples
sort_algebraic_numbers = sort_samples


def refine_pair_until_disj(left: AlgebraicRoot, right: AlgebraicRoot, *, steps: int = 8):
    lroot = left
    rroot = right
    for _ in range(max(steps, 0)):
        if lroot.interval.is_disjoint_from(rroot.interval):
            return lroot.interval, rroot.interval
        lroot = refine_isol_intv(lroot, steps=1)
        rroot = refine_isol_intv(rroot, steps=1)
    return lroot.interval, rroot.interval


def get_isolating_interval(root: AlgebraicRoot):
    return root.interval


__all__ = [
    "CACHE",
    "CACHE_VERSION",
    "RootIsolationCache",
    "RootIsolationStats",
    "root_isolation_costs",
    "RationalInterval",
    "RationalSample",
    "AlgebraicRoot",
    "Sample",
    "sample_to_expr",
    "isolate_real_roots",
    "refine_isol_intv",
    "root_multiplicity",
    "compare_samples",
    "sort_samples",
    "choose_sector_sample",
    "sign_at_sample",
    "SubresultantPRSResult",
    "subresultant_prs",
    "principal_subresultant_coefficients",
    "BorderBasisDiagnostics",
    "BorderBasisError",
    "BorderBasisResult",
    "compute_border_basis",
    "compute_border_basis_linear",
    "RationalUnivariateError",
    "RationalUnivariateRepresentation",
    "RationalUnivariatePoint",
    "FilteredRationalUnivariateSolutions",
    "RationalUnivariateFormulaResult",
    "compute_rational_univariate_representation",
    "solve_zero_dimensional_system_with_rur",
    "solve_rur_representation",
    "solve_rur_points",
    "sign_of_algebraic_expression",
    "evaluate_relation_at_point",
    "evaluate_boolean_formula_at_point",
    "filter_rur_solutions_by_constraints",
    "solve_rur_semialgebraic_system",
    "solve_and_filter_zero_dimensional_system_with_rur",
    "solve_formula_with_rur",
    "compare_alg_numbers",
    "sort_algebraic_numbers",
    "refine_pair_until_disj",
    "get_isolating_interval",
]
