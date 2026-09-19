from __future__ import annotations

from collections.abc import Sequence
from itertools import product

import sympy as sp
from sympy.polys.orderings import grevlex

from ._domains import require_exact_rur_domain
from .representation import RationalUnivariateError


def _as_exact_polynomial(
    expr: sp.Expr,
    variables: Sequence[sp.Symbol],
    *,
    domain=None,
) -> sp.Poly:
    """Return a polynomial over an exact rational/algebraic coefficient field."""

    try:
        poly = (
            sp.Poly(sp.expand(expr), *variables, domain=domain)
            if domain is not None
            else sp.Poly(sp.expand(expr), *variables, extension=True)
        )
    except (
        sp.PolynomialError,
        sp.polys.polyerrors.CoercionFailed,
        ValueError,
        TypeError,
        NotImplementedError,
    ) as exc:
        raise RationalUnivariateError(
            "RUR requires polynomial equations over an exact algebraic coefficient field"
        ) from exc
    if poly.total_degree() < 0:
        raise RationalUnivariateError("zero polynomial is not a valid defining equation")
    return poly.set_domain(require_exact_rur_domain(poly.domain))


def _leading_exponent_grevlex(poly: sp.Poly) -> tuple[int, ...]:
    terms = poly.terms(order=grevlex)
    if not terms:
        raise RationalUnivariateError("zero polynomial has no leading monomial")
    return tuple(int(v) for v in terms[0][0])


def _componentwise_leq(left: Sequence[int], right: Sequence[int]) -> bool:
    return all(a <= b for a, b in zip(left, right, strict=True))


def _is_not_divisible_by_any_leading_monomial(
    candidate: Sequence[int], leading_exponents: Sequence[Sequence[int]]
) -> bool:
    return not any(_componentwise_leq(leading, candidate) for leading in leading_exponents)


def _standard_exponent_count(
    leading_exponents: Sequence[Sequence[int]], variable_count: int
) -> int:
    """Count standard monomials without materializing their exponent vectors.

    The recursion prunes an entire suffix once a leading monomial already
    divides the assigned prefix, and counts an entire suffix at once once every
    leading monomial is permanently blocked by an assigned exponent.  This
    avoids the memory cliff of constructing the full standard-monomial basis
    and usually avoids scanning the complete pure-power bounding box as well.
    """

    bounds: list[int] = []
    for index in range(variable_count):
        pure_powers = [
            int(exponent[index])
            for exponent in leading_exponents
            if int(exponent[index]) > 0
            and all(int(value) == 0 for pos, value in enumerate(exponent) if pos != index)
        ]
        if not pure_powers:
            raise RationalUnivariateError("system does not expose a finite standard-monomial basis")
        bounds.append(min(pure_powers))

    generators = tuple(
        tuple(int(value) for value in exponent)
        for exponent in leading_exponents
        if any(int(value) for value in exponent)
    )
    # Remove generators divisible by another leading generator.
    minimal = tuple(
        generator
        for generator in generators
        if not any(
            other != generator and _componentwise_leq(other, generator) for other in generators
        )
    )
    suffix_sizes = [1] * (variable_count + 1)
    for index in range(variable_count - 1, -1, -1):
        suffix_sizes[index] = suffix_sizes[index + 1] * bounds[index]

    prefix: list[int] = []

    def count_from(index: int) -> int:
        # If an already fully supported generator divides the assigned prefix,
        # every completion lies in the leading ideal.
        for generator in minimal:
            if all(generator[pos] == 0 for pos in range(index, variable_count)) and all(
                generator[pos] <= prefix[pos] for pos in range(index)
            ):
                return 0

        # If every generator is already blocked by some assigned coordinate, no
        # later exponent can make any generator divide a completion.
        if minimal and all(
            any(prefix[pos] < generator[pos] for pos in range(index)) for generator in minimal
        ):
            return suffix_sizes[index]

        if index == variable_count:
            return 1
        total = 0
        for exponent in range(bounds[index]):
            prefix.append(exponent)
            total += count_from(index + 1)
            prefix.pop()
        return total

    return count_from(0)


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
    domain=sp.QQ,
) -> sp.Matrix:
    poly = sp.Poly(sp.expand(expr), *variables, domain=domain)
    coefficient_rules = {tuple(mon): coeff for mon, coeff in poly.terms()}
    return sp.Matrix(
        [coefficient_rules.get(tuple(exponent), sp.Integer(0)) for exponent in basis_exponents]
    )


def _multiplication_tensor(
    groebner_basis: sp.polys.polytools.GroebnerBasis,
    variables: Sequence[sp.Symbol],
    basis_exponents: Sequence[Sequence[int]],
    domain,
) -> list[list[sp.Matrix]]:
    basis_monomials = [_monomial_from_exponent(variables, exponent) for exponent in basis_exponents]
    tensor: list[list[sp.Matrix]] = []
    for left in basis_monomials:
        row: list[sp.Matrix] = []
        for right in basis_monomials:
            remainder = _normal_form(groebner_basis, left * right)
            row.append(_coefficient_vector(remainder, variables, basis_exponents, domain))
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
    return poly.monic()


def _squarefree_part(poly: sp.Poly) -> sp.Poly:
    derivative = poly.diff()
    gcd = sp.gcd(poly, derivative)
    return _monic_polynomial(poly.quo(gcd))


def _select_separating_linear_form(
    groebner_basis: sp.polys.polytools.GroebnerBasis,
    variables: Sequence[sp.Symbol],
    parameter: sp.Symbol,
    basis_exponents: Sequence[Sequence[int]],
    tensor: Sequence[Sequence[sp.Matrix]],
    domain,
    *,
    max_attempts: int = 64,
) -> tuple[sp.Expr, sp.Poly, sp.Poly, sp.Matrix, list[sp.Matrix], int]:
    """Choose a separating linear form whose quotient-algebra characteristic polynomial distinguishes all geometric solutions."""
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
        coordinate_vector = _coefficient_vector(remainder, variables, basis_exponents, domain)
        multiplication = _multiplication_matrix(coordinate_vector, tensor)
        characteristic = sp.Poly(
            multiplication.charpoly(parameter).as_expr(), parameter, domain=domain
        )
        squarefree = _squarefree_part(characteristic)
        if squarefree.degree() == expected_distinct_roots:
            derivative = characteristic.diff()
            gcd = sp.gcd(characteristic, derivative)
            denominator = derivative.quo(gcd)
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
