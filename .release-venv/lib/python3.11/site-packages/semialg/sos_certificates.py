"""Exact sum-of-squares certificate representation and verification.

``semialg`` owns exact verification. Optional optimization packages may search
numerically for Gram matrices, but a candidate is accepted only after it has
been converted to exact SymPy data and independently verified here.

The Gram/SDP formulation follows Parrilo [Parrilo2000, Parrilo2003]; see
``docs/references.md``.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass
from itertools import combinations
from math import comb
from typing import Any

import sympy as sp

from ._linear_relations import certified_sign
from ._zero_testing import certified_equal, certified_zero


@dataclass(frozen=True)
class SOSCertificate:
    """Exact Gram-matrix certificate ``p = z.T * Q * z`` with ``Q >= 0``."""

    polynomial: sp.Expr
    variables: tuple[sp.Symbol, ...]
    monomial_basis: tuple[sp.Expr, ...]
    gram_matrix: sp.ImmutableMatrix
    source: str = "exact"


@dataclass(frozen=True)
class SOSSearchPlan:
    """Cheap structural estimate deciding whether an automatic SOS search is worthwhile."""

    applicable: bool
    launch: bool
    degree: int
    variable_count: int
    gram_degree: int
    gram_dimension: int
    gram_variables: int
    reason: str
    dense_gram_dimension: int = 0
    basis_method: str = "dense"


@dataclass(frozen=True)
class SOSSearchResult:
    certificate: SOSCertificate | None
    complete: bool
    method: str
    notes: tuple[str, ...] = ()
    plan: SOSSearchPlan | None = None

    @property
    def certified(self) -> bool:
        return self.certificate is not None and self.complete


def _monomial_exponents(nvars: int, max_degree: int):
    if nvars == 0:
        yield ()
        return

    def rec(prefix: tuple[int, ...], left_vars: int, remaining: int):
        if left_vars == 1:
            yield (*prefix, remaining)
            return
        for value in range(remaining + 1):
            yield from rec((*prefix, value), left_vars - 1, remaining - value)

    for total in range(max_degree + 1):
        yield from rec((), nvars, total)


def _in_newton_polytope(target: tuple[int, ...], support: tuple[tuple[int, ...], ...]) -> bool:
    """Exact Caratheodory membership test for an integer point in ``conv(support)``."""
    if target in support:
        return True
    nvars = len(target)
    target_vec = sp.Matrix((*target, 1))
    max_size = min(nvars + 1, len(support))
    for size in range(2, max_size + 1):
        for points in combinations(support, size):
            mat = sp.Matrix(
                [[point[row] for point in points] for row in range(nvars)] + [[1] * size]
            )
            try:
                solution = sp.linsolve((mat, target_vec))
            except (ValueError, TypeError):
                continue
            if solution is sp.EmptySet:
                continue
            for values in solution:
                if any(value.free_symbols for value in values):
                    # A feasible rational point, if present, occurs at an endpoint;
                    # solve smaller affinely independent subsets instead.
                    continue
                if all(value.is_nonnegative is True for value in values):
                    return True
    return False


def sparse_sos_monomial_basis(
    polynomial: sp.Expr, variables: Sequence[sp.Symbol]
) -> tuple[sp.Expr, ...]:
    """Return the exact Newton-polytope Gram basis for an SOS search.

    If ``p = z.T Q z`` then every monomial exponent ``alpha`` occurring in
    ``z`` satisfies ``2*alpha in Newton(p)``.  Filtering the dense degree-half
    basis by this necessary condition is therefore lossless for SOS search.
    """
    vars_ = tuple(variables)
    poly = sp.Poly(sp.expand(polynomial), *vars_, domain=sp.QQ)
    if poly.is_zero:
        return (sp.Integer(1),)
    degree = int(poly.total_degree())
    if degree % 2:
        return tuple()
    support = tuple(tuple(int(v) for v in monom) for monom in poly.monoms())
    half_degree = degree // 2
    support_totals = {sum(monom) for monom in support}
    homogeneous = len(support_totals) == 1
    full_simplex = homogeneous and all(
        tuple(degree if i == j else 0 for i in range(len(vars_))) in support
        for j in range(len(vars_))
    )
    basis: list[sp.Expr] = []
    for alpha in _monomial_exponents(len(vars_), half_degree):
        if homogeneous and sum(alpha) != half_degree:
            continue
        if not full_simplex:
            doubled = tuple(2 * value for value in alpha)
            if any(
                doubled[i] < min(point[i] for point in support)
                or doubled[i] > max(point[i] for point in support)
                for i in range(len(vars_))
            ):
                continue
            if not _in_newton_polytope(doubled, support):
                continue
        term = sp.Integer(1)
        for var, exponent in zip(vars_, alpha, strict=True):
            term *= var**exponent
        basis.append(term)
    return tuple(basis)


def plan_sos_search(
    polynomial: sp.Expr,
    variables: Sequence[sp.Symbol],
    *,
    max_variables: int = 6,
    max_gram_dimension: int = 24,
    max_gram_variables: int = 300,
) -> SOSSearchPlan:
    """Estimate SOS search cost using the exact Newton-polytope basis."""
    vars_ = tuple(variables)
    expr = sp.expand(sp.sympify(polynomial))
    try:
        poly = sp.Poly(expr, *vars_, domain=sp.QQ) if vars_ else sp.Poly(expr, domain=sp.QQ)
    except (sp.PolynomialError, ValueError, TypeError):
        return SOSSearchPlan(False, False, -1, len(vars_), 0, 0, 0, "not_rational_polynomial")

    degree = int(poly.total_degree()) if not poly.is_zero else 0
    nvars = len(vars_)
    gram_degree = (degree + 1) // 2
    dense_dim = comb(nvars + gram_degree, gram_degree) if nvars else 1
    sparse_basis = sparse_sos_monomial_basis(expr, vars_) if degree % 2 == 0 else tuple()
    gram_dimension = len(sparse_basis) if sparse_basis else dense_dim
    gram_variables = gram_dimension * (gram_dimension + 1) // 2
    basis_method = "newton_polytope" if sparse_basis and gram_dimension < dense_dim else "dense"

    def result(applicable: bool, launch: bool, reason: str) -> SOSSearchPlan:
        return SOSSearchPlan(
            applicable,
            launch,
            degree,
            nvars,
            gram_degree,
            gram_dimension,
            gram_variables,
            reason,
            dense_gram_dimension=dense_dim,
            basis_method=basis_method,
        )

    if degree == 0:
        return result(True, False, "constant_exact_preferred")
    if degree % 2:
        return result(False, False, "odd_total_degree")
    if nvars <= 1:
        return result(True, False, "univariate_exact_preferred")
    if degree <= 2:
        return result(True, False, "quadratic_exact_preferred")
    if nvars > max_variables:
        return result(True, False, "variable_count_limit")
    if gram_dimension > max_gram_dimension:
        return result(True, False, "gram_dimension_limit")
    if gram_variables > max_gram_variables:
        return result(True, False, "gram_variable_limit")
    return result(True, True, "within_auto_budget")


def _exact_nonnegative(value: sp.Expr) -> bool:
    sign = certified_sign(sp.simplify(value))
    return sign in {0, 1}


@dataclass(frozen=True)
class PSDVerification:
    verified: bool
    method: str
    pivots: tuple[sp.Expr, ...] = ()
    reason: str = ""


def verify_psd_exact(matrix: sp.MatrixBase) -> PSDVerification:
    """Verify exact PSD by symmetric Schur-complement/LDL congruence elimination.

    A positive diagonal pivot preserves inertia under exact congruence.  A zero
    diagonal entry of a PSD matrix must have an identically zero row/column;
    this handles singular PSD matrices without dividing by zero.
    """
    work = sp.MutableDenseMatrix(matrix)
    if work.rows != work.cols or work != work.T:
        return PSDVerification(False, "ldl_congruence", reason="not_symmetric_square")
    pivots: list[sp.Expr] = []
    while work.rows:
        chosen = None
        zero_diagonals: list[int] = []
        for i in range(work.rows):
            diag = sp.simplify(work[i, i])
            sign = certified_sign(diag)
            if sign == -1:
                return PSDVerification(False, "ldl_congruence", tuple(pivots), "negative_pivot")
            if sign == 1:
                chosen = i
                break
            if sign == 0:
                zero_diagonals.append(i)
            else:
                return PSDVerification(False, "ldl_congruence", tuple(pivots), "unknown_pivot_sign")
        if chosen is None:
            # Every diagonal is zero. PSD then forces the complete matrix to be zero.
            if any(
                certified_zero(work[i, j]) is not True
                for i in range(work.rows)
                for j in range(work.cols)
            ):
                return PSDVerification(
                    False, "ldl_congruence", tuple(pivots), "zero_diagonal_nonzero_row"
                )
            pivots.extend(sp.Integer(0) for _ in range(work.rows))
            break
        if chosen != 0:
            work.row_swap(0, chosen)
            work.col_swap(0, chosen)
        pivot = sp.simplify(work[0, 0])
        pivots.append(pivot)
        if work.rows == 1:
            break
        column = work[1:, 0]
        tail = work[1:, 1:]
        scale = sp.cancel(1 / pivot)
        work = sp.MutableDenseMatrix(tail - (column * column.T) * scale).applyfunc(sp.cancel)
    return PSDVerification(True, "ldl_congruence", tuple(pivots))


def verify_sos_certificate(
    polynomial: sp.Expr,
    certificate: SOSCertificate,
    variables: Sequence[sp.Symbol] | None = None,
) -> bool:
    """Return whether an SOS certificate proves ``polynomial >= 0`` exactly.

    Malformed certificate payloads are rejected rather than leaking low-level
    SymPy conversion errors from a verification boundary.
    """

    if not isinstance(certificate, SOSCertificate):
        return False
    try:
        expr = sp.expand(sp.sympify(polynomial))
        cert_polynomial = sp.sympify(certificate.polynomial)
        cert_variables = tuple(certificate.variables)
        vars_ = tuple(variables) if variables is not None else cert_variables
        basis = tuple(map(sp.sympify, certificate.monomial_basis))
        gram = sp.ImmutableMatrix(certificate.gram_matrix)
    except (TypeError, ValueError, sp.SympifyError):
        return False
    if any(not isinstance(variable, sp.Symbol) for variable in (*cert_variables, *vars_)):
        return False
    if certified_equal(cert_polynomial, expr) is not True:
        return False
    if cert_variables != vars_:
        return False
    if gram.rows != len(basis) or gram.cols != len(basis):
        return False
    if any((term.free_symbols - set(vars_)) for term in basis):
        return False
    z = sp.Matrix(basis)
    residual = sp.expand(expr - (z.T * sp.Matrix(gram) * z)[0])
    if certified_zero(residual) is not True:
        return False
    return verify_psd_exact(gram).verified


def _exact_external_gram(result: Any) -> tuple[Any, str] | None:
    """Extract exact Gram data from the external certificate contract."""

    certificate = getattr(result, "certificate", None)
    if certificate is None:
        return None
    exact_residual_zero = getattr(certificate, "exact_residual_zero", None)
    psd_checks = getattr(certificate, "exact_psd_checks", None)
    gram_matrices = getattr(certificate, "gram_matrices", None)
    # Certificate metadata is only an applicability guard. semialg verifies
    # the reconstructed Gram identity and PSD condition independently.
    if exact_residual_zero is not True or gram_matrices is None:
        return None
    if psd_checks is not None and not all(check is True for check in psd_checks):
        return None
    try:
        blocks = [sp.ImmutableMatrix(block) for block in gram_matrices]
    except (TypeError, ValueError, sp.SympifyError):
        return None
    if not blocks or any(entry.has(sp.Float) for block in blocks for entry in block):
        return None
    gram = sp.diag(*blocks) if len(blocks) > 1 else blocks[0]
    return gram, "certificate"


def _candidate_from_external_result(
    polynomial: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    result: Any,
    *,
    source: str,
) -> SOSCertificate | None:
    if result is None or getattr(result, "success", True) is False:
        return None
    basis = getattr(result, "monomial_basis", None)
    extracted = _exact_external_gram(result)
    if basis is None or extracted is None:
        return None
    exact_gram, contract = extracted
    try:
        cert = SOSCertificate(
            sp.expand(polynomial),
            variables,
            tuple(map(sp.sympify, basis)),
            sp.ImmutableMatrix(exact_gram),
            f"{source}:{contract}",
        )
    except (TypeError, ValueError, sp.SympifyError):
        return None
    return cert if verify_sos_certificate(polynomial, cert, variables) else None


def search_sos_certificate(
    polynomial: sp.Expr,
    variables: Sequence[sp.Symbol],
    *,
    backend: str | Callable[[sp.Expr, list[sp.Symbol]], object] = "auto",
    use_planner: bool | None = None,
) -> SOSSearchResult:
    """Optionally ask ``symbopt`` for an SOS candidate, then verify it exactly.

    No search backend is required by semialg. Missing packages, numerical-only
    candidates, and failed searches are all represented as incomplete results.
    """

    vars_ = tuple(variables)
    custom_search = backend if callable(backend) else None
    backend_name = "custom" if custom_search is not None else backend.lower()
    external_names = {"cvxpy", "clarabel", "scs", "mosek"}
    if custom_search is None and backend_name not in {"auto", "symbopt", "none", *external_names}:
        raise ValueError("unsupported SOS backend")
    if backend_name == "none":
        return SOSSearchResult(None, False, "sos_skipped", ("search_disabled",))

    if use_planner is None:
        use_planner = backend_name == "auto"
    plan = plan_sos_search(polynomial, vars_) if use_planner else None
    if plan is not None and not plan.launch:
        notes = (
            f"planner:{plan.reason}",
            f"degree={plan.degree}",
            f"variables={plan.variable_count}",
            f"gram_dimension={plan.gram_dimension}",
            f"gram_variables={plan.gram_variables}",
            f"basis={plan.basis_method}",
            f"dense_gram_dimension={plan.dense_gram_dimension}",
        )
        return SOSSearchResult(None, False, "sos_planner_skipped", notes, plan)

    source = "external"
    search = custom_search
    if search is None and backend_name in external_names:
        from .sdp_backends import external_sdp_search

        source = backend_name

        def search(expr, syms):
            return external_sdp_search(expr, syms, backend=backend_name)

    if search is None:
        source = "symbopt"
        try:
            import symbopt  # type: ignore[import-not-found]
        except ImportError:
            return SOSSearchResult(None, False, "sos_symbopt", ("symbopt_unavailable",), plan)

        search = getattr(symbopt, "sos_decompose", None)
        if search is None:
            try:
                from symbopt.sos import sos_decompose as search  # type: ignore[import-not-found]
            except (ImportError, AttributeError):
                return SOSSearchResult(
                    None, False, "sos_symbopt", ("sos_decompose_unavailable",), plan
                )

    try:
        external = search(sp.expand(polynomial), list(vars_))
    except (ArithmeticError, ValueError, TypeError, RuntimeError, AttributeError) as exc:
        return SOSSearchResult(
            None,
            False,
            f"sos_{source}",
            (f"search_failed:{type(exc).__name__}",),
            plan,
        )
    certificate = _candidate_from_external_result(polynomial, vars_, external, source=source)
    if certificate is None:
        notes = ["candidate_not_exactly_verified"]
        status = getattr(external, "status", None)
        if status:
            notes.append(f"backend_status:{status}")
        return SOSSearchResult(None, False, f"sos_{source}", tuple(notes), plan)
    return SOSSearchResult(certificate, True, f"sos_{source}_verified", plan=plan)


__all__ = [
    "SOSCertificate",
    "SOSSearchPlan",
    "SOSSearchResult",
    "PSDVerification",
    "plan_sos_search",
    "sparse_sos_monomial_basis",
    "search_sos_certificate",
    "verify_psd_exact",
    "verify_sos_certificate",
]
