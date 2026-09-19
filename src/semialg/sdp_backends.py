"""Optional numerical SDP proposal backends for exact SOS certification."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from fractions import Fraction

import sympy as sp


@dataclass(frozen=True)
class ExternalSOSCandidate:
    success: bool
    monomial_basis: tuple[sp.Expr, ...]
    gram_matrix: sp.ImmutableMatrix | None
    backend: str
    status: str


def _coefficient_equations(
    polynomial: sp.Expr, variables: tuple[sp.Symbol, ...], basis: tuple[sp.Expr, ...]
):
    pairs = [(i, j) for i in range(len(basis)) for j in range(i, len(basis))]
    symbols = sp.symbols(f"q0:{len(pairs)}")
    gram_expr = sp.Integer(0)
    for symbol, (i, j) in zip(symbols, pairs, strict=True):
        factor = 1 if i == j else 2
        gram_expr += factor * symbol * basis[i] * basis[j]
    residual = sp.Poly(
        sp.expand(gram_expr - polynomial), *variables, domain=sp.QQ.frac_field(*symbols)
    )
    equations = [sp.expand(coeff) for coeff in residual.coeffs()]
    # Poly over a fraction field may omit zero coefficients; coefficient equations suffice.
    return pairs, symbols, equations


def reconstruct_exact_gram(
    polynomial: sp.Expr,
    variables: Sequence[sp.Symbol],
    basis: Sequence[sp.Expr],
    numeric_gram,
    *,
    max_denominator: int = 100000,
) -> sp.ImmutableMatrix | None:
    """Project a numerical Gram candidate onto the exact coefficient affine space."""
    vars_ = tuple(variables)
    basis_ = tuple(map(sp.sympify, basis))
    pairs = [(i, j) for i in range(len(basis_)) for j in range(i, len(basis_))]
    q = sp.symbols(f"q0:{len(pairs)}")
    expr = sp.Integer(0)
    for symbol, (i, j) in zip(q, pairs, strict=True):
        expr += (1 if i == j else 2) * symbol * basis_[i] * basis_[j]
    poly = sp.Poly(sp.expand(expr - polynomial), *vars_)
    equations = [sp.Eq(coeff, 0) for coeff in poly.coeffs()]
    solution = sp.linsolve([eq.lhs for eq in equations], q)
    if solution is sp.EmptySet:
        return None
    vector = next(iter(solution))
    free = sorted(set().union(*(value.free_symbols for value in vector)), key=sp.default_sort_key)
    numeric_values = {}
    for symbol in free:
        # Free symbols are some q_i. Use the corresponding numerical entry as a rational hint.
        try:
            index = q.index(symbol)
            i, j = pairs[index]
            value = float(numeric_gram[i][j])
        except (ValueError, TypeError, IndexError):
            value = 0.0
        numeric_values[symbol] = sp.Rational(Fraction(value).limit_denominator(max_denominator))
    exact = [sp.cancel(value.subs(numeric_values)) for value in vector]
    matrix = sp.zeros(len(basis_))
    for value, (i, j) in zip(exact, pairs, strict=True):
        matrix[i, j] = value
        matrix[j, i] = value
    return sp.ImmutableMatrix(matrix)


def external_sdp_search(
    polynomial: sp.Expr,
    variables: Sequence[sp.Symbol],
    *,
    backend: str = "cvxpy",
) -> ExternalSOSCandidate:
    """Ask CVXPY-compatible solvers for a numerical Gram proposal, then recover it exactly."""
    from .sos_certificates import sparse_sos_monomial_basis

    vars_ = tuple(variables)
    basis = sparse_sos_monomial_basis(polynomial, vars_)
    if not basis:
        return ExternalSOSCandidate(False, basis, None, backend, "empty_basis")
    try:
        import cvxpy as cp  # type: ignore[import-not-found]
        import numpy as np
    except ImportError:
        return ExternalSOSCandidate(False, basis, None, backend, "cvxpy_unavailable")

    size = len(basis)
    q = cp.Variable((size, size), symmetric=True)
    constraints = [q >> 0]
    target = sp.Poly(sp.expand(polynomial), *vars_, domain=sp.QQ)
    exponents = set(target.monoms())
    products: dict[tuple[int, ...], list[tuple[int, int, int]]] = {}
    basis_exp = [sp.Poly(term, *vars_).monoms()[0] for term in basis]
    for i in range(size):
        for j in range(i, size):
            exponent = tuple(a + b for a, b in zip(basis_exp[i], basis_exp[j], strict=True))
            products.setdefault(exponent, []).append((i, j, 1 if i == j else 2))
            exponents.add(exponent)
    coeffs = dict(target.terms())
    for exponent in exponents:
        lhs = sum(mult * q[i, j] for i, j, mult in products.get(exponent, ()))
        constraints.append(lhs == float(coeffs.get(exponent, 0)))
    problem = cp.Problem(cp.Minimize(cp.trace(q)), constraints)
    solver = None if backend == "cvxpy" else backend.upper()
    try:
        problem.solve(solver=solver, verbose=False)
    except (cp.error.SolverError, ValueError):
        return ExternalSOSCandidate(False, basis, None, backend, "solver_unavailable_or_failed")
    if q.value is None or problem.status not in {cp.OPTIMAL, cp.OPTIMAL_INACCURATE}:
        return ExternalSOSCandidate(False, basis, None, backend, str(problem.status))
    numeric = np.asarray(q.value)
    exact = reconstruct_exact_gram(polynomial, vars_, basis, numeric)
    if exact is None:
        return ExternalSOSCandidate(False, basis, None, backend, "exact_recovery_failed")
    return ExternalSOSCandidate(True, basis, exact, backend, str(problem.status))


__all__ = ["ExternalSOSCandidate", "external_sdp_search", "reconstruct_exact_gram"]
