from __future__ import annotations

from collections.abc import Sequence
from itertools import product

import sympy as sp

from .representation import RationalUnivariateError


def _as_rational_polynomial(expr: sp.Expr, variables: Sequence[sp.Symbol]) -> sp.Poly:
    poly = sp.Poly(sp.expand(expr), *variables, domain=sp.QQ)
    if poly.total_degree() < 0:
        raise RationalUnivariateError("zero polynomial is not a valid defining equation")
    return poly


def _leading_exponent_grevlex(poly: sp.Poly) -> tuple[int, ...]:
    terms = poly.terms(order="grevlex")
    if not terms:
        raise RationalUnivariateError("zero polynomial has no leading monomial")
    return tuple(int(v) for v in terms[0][0])


def _componentwise_leq(left: Sequence[int], right: Sequence[int]) -> bool:
    return all(a <= b for a, b in zip(left, right, strict=True))


def _is_not_divisible_by_any_leading_monomial(
    candidate: Sequence[int], leading_exponents: Sequence[Sequence[int]]
) -> bool:
    return not any(_componentwise_leq(leading, candidate) for leading in leading_exponents)


def _standard_exponents(
    leading_exponents: Sequence[Sequence[int]], variable_count: int
) -> tuple[tuple[int, ...], ...]:
    """Return exponent vectors of monomials outside the leading ideal.

    The standard monomial basis is finite when each variable direction has a pure-power leading monomial bound. This condition is sufficient for the quotient basis used
    here and rejects positive-dimensional inputs early.
    """

    bounds: list[int] = []
    for index in range(variable_count):
        pure_power = None
        for exponent in leading_exponents:
            if all(value == 0 for pos, value in enumerate(exponent) if pos != index):
                pure_power = int(exponent[index])
                break
        if pure_power is None or pure_power <= 0:
            raise RationalUnivariateError("system does not expose a finite standard-monomial basis")
        bounds.append(pure_power)

    candidates = product(*(range(bound) for bound in bounds))
    basis = [
        tuple(candidate)
        for candidate in candidates
        if _is_not_divisible_by_any_leading_monomial(candidate, leading_exponents)
    ]
    if not basis or basis[0] != tuple(0 for _ in range(variable_count)):
        basis.sort(key=lambda exp: (sum(exp), exp))
    return tuple(basis)


def _monomial_from_exponent(variables: Sequence[sp.Symbol], exponent: Sequence[int]) -> sp.Expr:
    monomial = sp.Integer(1)
    for variable, power in zip(variables, exponent, strict=True):
        monomial *= variable ** int(power)
    return monomial


def _normal_form(groebner_basis: sp.polys.polytools.GroebnerBasis, expr: sp.Expr) -> sp.Expr:
    try:
        _, remainder = groebner_basis.reduce(sp.expand(expr))
    except (
        sp.PolynomialError,
        ValueError,
        TypeError,
    ) as exc:  # pragma: no cover - defensive SymPy boundary
        raise RationalUnivariateError(f"Groebner reduction failed: {exc}") from exc
    return sp.expand(remainder)


def _coefficient_vector(
    expr: sp.Expr,
    variables: Sequence[sp.Symbol],
    basis_exponents: Sequence[Sequence[int]],
) -> sp.Matrix:
    poly = sp.Poly(sp.expand(expr), *variables, domain=sp.QQ)
    coefficient_rules = {tuple(mon): coeff for mon, coeff in poly.terms()}
    return sp.Matrix(
        [coefficient_rules.get(tuple(exponent), sp.Integer(0)) for exponent in basis_exponents]
    )


def _multiplication_tensor(
    groebner_basis: sp.polys.polytools.GroebnerBasis,
    variables: Sequence[sp.Symbol],
    basis_exponents: Sequence[Sequence[int]],
) -> list[list[sp.Matrix]]:
    basis_monomials = [_monomial_from_exponent(variables, exponent) for exponent in basis_exponents]
    tensor: list[list[sp.Matrix]] = []
    for left in basis_monomials:
        row: list[sp.Matrix] = []
        for right in basis_monomials:
            remainder = _normal_form(groebner_basis, left * right)
            row.append(_coefficient_vector(remainder, variables, basis_exponents))
        tensor.append(row)
    return tensor


def _multiplication_matrix(
    coordinates: sp.Matrix, tensor: Sequence[Sequence[sp.Matrix]]
) -> sp.Matrix:
    dimension = len(tensor)
    columns: list[sp.Matrix] = []
    for basis_index in range(dimension):
        column = sp.zeros(dimension, 1)
        for coeff_index in range(dimension):
            column += coordinates[coeff_index] * tensor[coeff_index][basis_index]
        columns.append(column)
    return sp.Matrix.hstack(*columns) if columns else sp.zeros(0, 0)


def _monic_polynomial(poly: sp.Poly) -> sp.Poly:
    if poly.is_zero:
        raise RationalUnivariateError("zero polynomial cannot be normalized to monic form")
    return sp.Poly(poly.as_expr() / poly.LC(), *poly.gens, domain=sp.QQ)


def _squarefree_part(poly: sp.Poly) -> sp.Poly:
    derivative = poly.diff()
    gcd = sp.gcd(poly, derivative)
    return _monic_polynomial(
        sp.Poly(sp.cancel(poly.as_expr() / gcd.as_expr()), *poly.gens, domain=sp.QQ)
    )


def _select_separating_linear_form(
    groebner_basis: sp.polys.polytools.GroebnerBasis,
    variables: Sequence[sp.Symbol],
    parameter: sp.Symbol,
    basis_exponents: Sequence[Sequence[int]],
    tensor: Sequence[Sequence[sp.Matrix]],
    *,
    max_attempts: int = 64,
) -> tuple[sp.Expr, sp.Poly, sp.Poly, sp.Matrix, list[sp.Matrix], int]:
    dimension = len(basis_exponents)
    trace_of_basis_mult = sp.Matrix(
        [sum(tensor[i][j][j] for j in range(dimension)) for i in range(dimension)]
    )
    trace_pairing = sp.zeros(dimension, dimension)
    for left in range(dimension):
        for right in range(dimension):
            trace_pairing[left, right] = (tensor[left][right].T * trace_of_basis_mult)[0]
    expected_distinct_roots = trace_pairing.rank()

    if max_attempts < 1:
        raise RationalUnivariateError("max_attempts must be positive")

    for attempt in range(1, max_attempts + 1):
        linear_form = sum((attempt**idx) * variable for idx, variable in enumerate(variables))
        remainder = _normal_form(groebner_basis, linear_form)
        coordinate_vector = _coefficient_vector(remainder, variables, basis_exponents)
        multiplication = _multiplication_matrix(coordinate_vector, tensor)
        characteristic = sp.Poly(
            multiplication.charpoly(parameter).as_expr(), parameter, domain=sp.QQ
        )
        squarefree = _squarefree_part(characteristic)
        if squarefree.degree() == expected_distinct_roots:
            derivative = characteristic.diff()
            gcd = sp.gcd(characteristic, derivative)
            denominator = sp.Poly(
                sp.cancel(derivative.as_expr() / gcd.as_expr()), parameter, domain=sp.QQ
            )
            powers = [sp.eye(dimension).col(0)]
            for _ in range(1, squarefree.degree()):
                powers.append(sp.simplify(multiplication * powers[-1]))
            return (
                linear_form,
                squarefree,
                denominator,
                trace_of_basis_mult,
                powers,
                expected_distinct_roots,
            )
    raise RationalUnivariateError("could not find a separating linear form")
