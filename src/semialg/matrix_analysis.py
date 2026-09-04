"""Exact definiteness and rank analysis for symmetric semialgebraic matrices."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from itertools import combinations

import sympy as sp

from .conditional import ParameterStratifiedResult, conditional_result
from .decision import is_satisfiable
from .normalization import normalize_formula, normalize_problem_variables, normalize_variables
from .parameters import solvability_conditions
from .sampling import sample_point


@dataclass(frozen=True)
class MatrixDefinitenessResult:
    """Exact result for a symmetric-matrix definiteness query on a region."""

    outcome: bool | None
    requested: str
    matrix: sp.ImmutableMatrix
    variables: tuple[sp.Symbol, ...]
    domain: sp.Expr
    method: str
    inertia: tuple[int, int, int] | None = None
    minors: tuple[sp.Expr, ...] = ()
    counterexample: Mapping[sp.Symbol, sp.Expr] | None = None
    counterexample_vector: tuple[sp.Expr, ...] | None = None
    diagnostics: Mapping[str, object] = field(default_factory=dict)

    @property
    def certified(self) -> bool:
        return self.outcome is not None


@dataclass(frozen=True)
class MatrixRankResult:
    """Exact constant-rank result for a matrix on a semialgebraic region."""

    rank: int | None
    matrix: sp.ImmutableMatrix
    variables: tuple[sp.Symbol, ...]
    domain: sp.Expr
    method: str
    witness: Mapping[sp.Symbol, sp.Expr] | None = None
    diagnostics: Mapping[str, object] = field(default_factory=dict)

    @property
    def certified(self) -> bool:
        return self.rank is not None


def _exact_constant_sign(value: sp.Expr) -> int | None:
    simplified = sp.simplify(value)
    if simplified == 0:
        return 0
    if simplified.is_positive is True:
        return 1
    if simplified.is_negative is True:
        return -1
    try:
        sign = sp.sign(simplified)
    except (TypeError, ValueError, NotImplementedError):
        return None
    return int(sign) if sign in (-1, 0, 1) else None


def constant_symmetric_inertia(matrix: sp.MatrixBase) -> tuple[int, int, int] | None:
    """Return exact inertia ``(positive, negative, zero)`` by congruence LDL elimination."""

    immutable_matrix = sp.ImmutableMatrix(matrix)
    if immutable_matrix.rows != immutable_matrix.cols or immutable_matrix != immutable_matrix.T:
        return None
    if any(entry.free_symbols for entry in immutable_matrix):
        return None
    working_matrix = sp.MutableDenseMatrix(immutable_matrix)
    positive_count = negative_count = zero_count = 0
    while working_matrix.rows:
        pivot_index = None
        pivot_sign = None
        for index in range(working_matrix.rows):
            diagonal_sign = _exact_constant_sign(working_matrix[index, index])
            if diagonal_sign is None:
                return None
            if diagonal_sign != 0:
                pivot_index = index
                pivot_sign = diagonal_sign
                break
        if pivot_index is None:
            pivot_pair = None
            for row in range(working_matrix.rows):
                for column in range(row + 1, working_matrix.cols):
                    off_diagonal_sign = _exact_constant_sign(working_matrix[row, column])
                    if off_diagonal_sign is None:
                        return None
                    if off_diagonal_sign != 0:
                        pivot_pair = (row, column)
                        break
                if pivot_pair is not None:
                    break
            if pivot_pair is None:
                zero_count += working_matrix.rows
                break
            row, column = pivot_pair
            if row != 0:
                working_matrix.row_swap(0, row)
                working_matrix.col_swap(0, row)
                if column == 0:
                    column = row
                elif column == row:
                    column = 0
            if column != 1:
                working_matrix.row_swap(1, column)
                working_matrix.col_swap(1, column)
            positive_count += 1
            negative_count += 1
            if working_matrix.rows == 2:
                break
            pivot_block = working_matrix[:2, :2]
            coupling = working_matrix[2:, :2]
            trailing_block = working_matrix[2:, 2:]
            schur_complement = trailing_block - coupling * pivot_block.inv() * coupling.T
            working_matrix = sp.MutableDenseMatrix(schur_complement.applyfunc(sp.cancel))
            continue
        if pivot_index != 0:
            working_matrix.row_swap(0, pivot_index)
            working_matrix.col_swap(0, pivot_index)
        pivot = sp.simplify(working_matrix[0, 0])
        if pivot_sign == 1:
            positive_count += 1
        else:
            negative_count += 1
        if working_matrix.rows == 1:
            break
        trailing_block = working_matrix[1:, 1:]
        pivot_column = working_matrix[1:, 0]
        schur_complement = trailing_block - (pivot_column * pivot_column.T) / pivot
        working_matrix = sp.MutableDenseMatrix(schur_complement.applyfunc(sp.cancel))
    return positive_count, negative_count, zero_count


def principal_minors(matrix: sp.MatrixBase) -> tuple[sp.Expr, ...]:
    """Return all principal minors in deterministic size/index order."""

    minors: list[sp.Expr] = []
    for size in range(1, matrix.rows + 1):
        for indices in combinations(range(matrix.rows), size):
            minors.append(sp.factor(matrix.extract(indices, indices).det()))
    return tuple(minors)


def leading_principal_minors(matrix: sp.MatrixBase) -> tuple[sp.Expr, ...]:
    return tuple(sp.factor(matrix[:size, :size].det()) for size in range(1, matrix.rows + 1))


def _criterion(
    matrix: sp.ImmutableMatrix, requested: str
) -> tuple[tuple[sp.Expr, ...], tuple[sp.Expr, ...]]:
    """Return ``(minors, violations)`` for the requested definiteness property."""

    if requested == "positive_semidefinite":
        minors = principal_minors(matrix)
        return minors, tuple(minor < 0 for minor in minors)
    if requested == "positive_definite":
        minors = leading_principal_minors(matrix)
        return minors, tuple(minor <= 0 for minor in minors)
    if requested == "negative_semidefinite":
        return _criterion(sp.ImmutableMatrix(-matrix), "positive_semidefinite")
    if requested == "negative_definite":
        return _criterion(sp.ImmutableMatrix(-matrix), "positive_definite")
    raise ValueError(
        "requested must be 'positive_semidefinite', 'positive_definite', "
        "'negative_semidefinite', or 'negative_definite'"
    )


def _constant_outcome(inertia: tuple[int, int, int], requested: str) -> bool:
    positive_count, negative_count, zero_count = inertia
    if requested == "positive_semidefinite":
        return negative_count == 0
    if requested == "positive_definite":
        return negative_count == 0 and zero_count == 0
    if requested == "negative_semidefinite":
        return positive_count == 0
    return positive_count == 0 and zero_count == 0


def _quadratic_form_counterexample_vector(
    matrix: sp.ImmutableMatrix, requested: str
) -> tuple[sp.Expr, ...] | None:
    """Best-effort exact vector violating a constant definiteness property."""
    signed_matrix = matrix if requested.startswith("positive") else sp.ImmutableMatrix(-matrix)
    allow_zero = requested.endswith("definite") and not requested.endswith("semidefinite")
    try:
        eigenvectors = signed_matrix.eigenvects()
    except (TypeError, ValueError, NotImplementedError, sp.PolynomialError):
        eigenvectors = ()
    for eigenvalue, _multiplicity, vectors in eigenvectors:
        eigenvalue_sign = _exact_constant_sign(eigenvalue)
        violates = eigenvalue_sign == -1 or (allow_zero and eigenvalue_sign == 0)
        if violates and vectors:
            vector = vectors[0]
            return tuple(sp.simplify(vector[index]) for index in range(vector.rows))
    return None


def matrix_definiteness(
    matrix: sp.MatrixBase | Sequence[Sequence[object]],
    variables: Sequence[sp.Symbol | str] | None = None,
    *,
    domain=sp.true,
    requested: str = "positive_semidefinite",
    parameters: Sequence[sp.Symbol | str] | None = None,
    return_result: bool = False,
) -> bool | MatrixDefinitenessResult | ParameterStratifiedResult:
    """Decide exact symmetric-matrix definiteness on a semialgebraic domain.

    With ``parameters`` supplied, return certified Boolean parameter strata on
    which the requested definiteness property holds or fails.
    """

    immutable_matrix = sp.ImmutableMatrix(matrix)
    if immutable_matrix.rows != immutable_matrix.cols:
        raise ValueError("matrix definiteness requires a square matrix")
    if immutable_matrix != immutable_matrix.T:
        raise ValueError("matrix definiteness requires a symmetric matrix")
    normalized_domain = normalize_formula(domain)
    matrix_symbols = set().union(*(entry.free_symbols for entry in immutable_matrix))
    combined_expression = sp.Tuple(*immutable_matrix, normalized_domain)
    normalized_variables = (
        normalize_problem_variables(None, combined_expression)
        if variables is None and parameters is None
        else normalize_variables(variables, combined_expression, append_context_symbols=False)
    )
    normalized_parameters = (
        normalize_variables(parameters, combined_expression, append_context_symbols=False)
        if parameters is not None
        else ()
    )
    if set(normalized_variables) & set(normalized_parameters):
        raise ValueError("variables and parameters must be disjoint")

    if not normalized_variables and not normalized_parameters and not matrix_symbols:
        inertia = constant_symmetric_inertia(immutable_matrix)
        outcome = None if inertia is None else _constant_outcome(inertia, requested)
        result = MatrixDefinitenessResult(
            outcome,
            requested,
            immutable_matrix,
            (),
            normalized_domain,
            "constant_ldlt_inertia" if inertia is not None else "constant_inertia_unknown",
            inertia=inertia,
        )
        return result if return_result else bool(outcome) if outcome is not None else False

    minors, violations = _criterion(immutable_matrix, requested)
    violation_formula = sp.And(normalized_domain, sp.Or(*violations)) if violations else sp.false

    if normalized_parameters:
        violation_condition = solvability_conditions(
            violation_formula,
            normalized_variables,
            normalized_parameters,
        )
        property_condition = sp.simplify_logic(sp.Not(violation_condition), force=True)
        return conditional_result(
            normalized_parameters,
            ((property_condition, True), (sp.Not(property_condition), False)),
            method=f"matrix_{requested}_parameter_qe",
            diagnostics={"minors": tuple(sp.sstr(minor) for minor in minors)},
        )

    satisfiability = is_satisfiable(violation_formula, normalized_variables, return_result=True)
    outcome = not bool(satisfiability)
    counterexample = satisfiability.witness if not outcome else None
    counterexample_vector = None
    if counterexample is not None:
        specialized_matrix = sp.ImmutableMatrix(
            [
                [sp.simplify(entry.subs(counterexample)) for entry in immutable_matrix.row(row)]
                for row in range(immutable_matrix.rows)
            ]
        )
        if not any(entry.free_symbols for entry in specialized_matrix):
            counterexample_vector = _quadratic_form_counterexample_vector(
                specialized_matrix, requested
            )
    result = MatrixDefinitenessResult(
        outcome,
        requested,
        immutable_matrix,
        normalized_variables,
        normalized_domain,
        "principal_minor_exact_feasibility",
        minors=minors,
        counterexample=counterexample,
        counterexample_vector=counterexample_vector,
        diagnostics={"violation_formula": sp.sstr(violation_formula)},
    )
    return result if return_result else outcome


def matrix_psd_on(matrix, variables=None, *, domain=sp.true, parameters=None, return_result=False):
    """Decide whether a symmetric matrix is positive semidefinite on a region."""
    return matrix_definiteness(
        matrix,
        variables,
        domain=domain,
        requested="positive_semidefinite",
        parameters=parameters,
        return_result=return_result,
    )


def matrix_pd_on(matrix, variables=None, *, domain=sp.true, parameters=None, return_result=False):
    """Decide whether a symmetric matrix is positive definite on a region."""
    return matrix_definiteness(
        matrix,
        variables,
        domain=domain,
        requested="positive_definite",
        parameters=parameters,
        return_result=return_result,
    )


def _all_minors(matrix: sp.ImmutableMatrix, size: int) -> tuple[sp.Expr, ...]:
    if size == 0:
        return (sp.Integer(1),)
    return tuple(
        sp.factor(matrix.extract(rows, columns).det())
        for rows in combinations(range(matrix.rows), size)
        for columns in combinations(range(matrix.cols), size)
    )


def _rank_exact_formula(matrix: sp.ImmutableMatrix, rank: int) -> sp.Expr:
    max_rank = min(matrix.rows, matrix.cols)
    if not 0 <= rank <= max_rank:
        return sp.false
    rank_at_least = (
        sp.Or(*(sp.Ne(minor, 0) for minor in _all_minors(matrix, rank))) if rank else sp.true
    )
    rank_at_most = (
        sp.And(*(sp.Eq(minor, 0) for minor in _all_minors(matrix, rank + 1)))
        if rank < max_rank
        else sp.true
    )
    return sp.And(rank_at_least, rank_at_most)


def matrix_rank_stratification(
    matrix: sp.MatrixBase | Sequence[Sequence[object]],
    parameters: Sequence[sp.Symbol | str],
    *,
    parameter_domain=sp.true,
) -> ParameterStratifiedResult:
    """Partition parameter space by the exact rank of a symbolic matrix."""

    immutable_matrix = sp.ImmutableMatrix(matrix)
    domain = normalize_formula(parameter_domain)
    normalized_parameters = normalize_variables(
        parameters, sp.Tuple(*immutable_matrix, domain), append_context_symbols=False
    )
    branches = []
    for rank in range(min(immutable_matrix.rows, immutable_matrix.cols) + 1):
        condition = sp.And(domain, _rank_exact_formula(immutable_matrix, rank))
        branches.append((condition, rank))
    return conditional_result(
        normalized_parameters,
        branches,
        coverage_condition=domain,
        method="determinantal_rank_stratification",
    )


def matrix_rank_on(
    matrix: sp.MatrixBase | Sequence[Sequence[object]],
    variables: Sequence[sp.Symbol | str] | None = None,
    *,
    domain=sp.true,
    return_result: bool = False,
) -> int | None | MatrixRankResult:
    """Return the exact constant rank on ``domain``, or ``None`` if rank varies."""

    immutable_matrix = sp.ImmutableMatrix(matrix)
    normalized_domain = normalize_formula(domain)
    normalized_variables = normalize_problem_variables(
        variables, sp.Tuple(*immutable_matrix, normalized_domain)
    )
    if not is_satisfiable(normalized_domain, normalized_variables):
        result = MatrixRankResult(
            None, immutable_matrix, normalized_variables, normalized_domain, "empty_domain"
        )
        return result if return_result else None
    max_rank = min(immutable_matrix.rows, immutable_matrix.cols)
    for rank in range(max_rank + 1):
        violation = sp.And(normalized_domain, sp.Not(_rank_exact_formula(immutable_matrix, rank)))
        if not is_satisfiable(violation, normalized_variables):
            result = MatrixRankResult(
                rank,
                immutable_matrix,
                normalized_variables,
                normalized_domain,
                "determinantal_exact_rank",
            )
            return result if return_result else rank
    witness = sample_point(normalized_domain, normalized_variables)
    result = MatrixRankResult(
        None,
        immutable_matrix,
        normalized_variables,
        normalized_domain,
        "rank_varies",
        witness=witness,
    )
    return result if return_result else None


__all__ = [
    "MatrixDefinitenessResult",
    "MatrixRankResult",
    "constant_symmetric_inertia",
    "leading_principal_minors",
    "matrix_definiteness",
    "matrix_pd_on",
    "matrix_psd_on",
    "matrix_rank_on",
    "matrix_rank_stratification",
    "principal_minors",
]
