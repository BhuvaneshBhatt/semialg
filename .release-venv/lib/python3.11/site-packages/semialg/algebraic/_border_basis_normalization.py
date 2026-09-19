"""Normalization and monomial-order helpers for exact border bases.

These helpers are kept out of the realization module so representation and
normalization policy remain independently testable without wrapping hot loops.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence

import sympy as sp

from ._border_basis_types import BorderBasisError
from .rational_univariate.quotient import (
    _coefficient_vector,
    _monomial_from_exponent,
    _normal_form,
)


def _sort_exponents(exponents: Iterable[Sequence[int]]) -> tuple[tuple[int, ...], ...]:
    """Sort monomial exponents by degree, preferring earlier variables.

    SymPy's grevlex standard-monomial enumeration can leave symmetric quotient
    bases in either ``x`` or ``y``. For a public border-basis object it is more
    predictable to prefer monomials involving earlier user-supplied variables,
    so ``x`` precedes ``y`` for variables ``[x, y]`` at the same total degree.
    """

    normalized = {tuple(int(value) for value in exp) for exp in exponents}
    return tuple(sorted(normalized, key=lambda exp: (sum(exp), tuple(-value for value in exp))))


def _degree_exponents(variable_count: int, degree: int) -> tuple[tuple[int, ...], ...]:
    if variable_count == 1:
        return ((degree,),)
    out: list[tuple[int, ...]] = []
    for head in range(degree, -1, -1):
        for tail in _degree_exponents(variable_count - 1, degree - head):
            out.append((head,) + tail)
    return tuple(out)


def _is_divisor_closed_after_add(order: set[tuple[int, ...]], exponent: tuple[int, ...]) -> bool:
    for index, power in enumerate(exponent):
        if power <= 0:
            continue
        divisor = list(exponent)
        divisor[index] -= 1
        if tuple(divisor) not in order:
            return False
    return True


def _preferred_order_ideal_from_quotient(
    groebner_basis: sp.polys.polytools.GroebnerBasis,
    variables: Sequence[sp.Symbol],
    standard_order: Sequence[Sequence[int]],
) -> tuple[tuple[int, ...], ...]:
    """Choose a user-variable-friendly order ideal for the quotient basis.

    Standard monomials of a Groebner basis are a valid order ideal but are tied
    to the supporting monomial order. For public border-basis output we prefer
    an equivalent order ideal whose normal forms are linearly independent and
    whose monomials use earlier variables when possible. This keeps examples
    such as ``[x**2 - 1, y - x]`` in the intuitive basis ``(1, x)``.
    """

    dimension = len(tuple(standard_order))
    if dimension <= 1:
        return tuple(tuple(exp) for exp in standard_order)
    selected: list[tuple[int, ...]] = []
    selected_set: set[tuple[int, ...]] = set()
    columns: list[sp.Matrix] = []
    degree = 0
    max_degree = max(sum(exp) for exp in standard_order) + dimension + 1
    while len(selected) < dimension and degree <= max_degree:
        for exponent in _degree_exponents(len(variables), degree):
            if exponent in selected_set or not _is_divisor_closed_after_add(selected_set, exponent):
                continue
            monomial = _monomial_from_exponent(variables, exponent)
            remainder = _normal_form(groebner_basis, monomial)
            vector = _coefficient_vector(remainder, variables, standard_order)
            trial = columns + [vector]
            if sp.Matrix.hstack(*trial).rank() == len(trial):
                selected.append(exponent)
                selected_set.add(exponent)
                columns.append(vector)
                if len(selected) == dimension:
                    return tuple(selected)
        degree += 1
    return tuple(tuple(exp) for exp in standard_order)


def _border_exponents(
    order_ideal: Sequence[Sequence[int]],
    variable_count: int,
) -> tuple[tuple[int, ...], ...]:
    order_set = {tuple(exp) for exp in order_ideal}
    border: set[tuple[int, ...]] = set()
    for exponent in order_set:
        for index in range(variable_count):
            candidate = list(exponent)
            candidate[index] += 1
            candidate_tuple = tuple(candidate)
            if candidate_tuple not in order_set:
                border.add(candidate_tuple)
    return _sort_exponents(border)


def _exponent_from_monomial(monomial: sp.Expr, variables: Sequence[sp.Symbol]) -> tuple[int, ...]:
    poly = sp.Poly(monomial, *variables, domain=sp.QQ)
    terms = poly.terms()
    if len(terms) != 1 or terms[0][1] != 1:
        raise BorderBasisError(f"order ideal entry is not a monomial: {monomial!s}")
    return tuple(int(value) for value in terms[0][0])


def _normalize_order_ideal(
    order_ideal: Iterable[Sequence[int] | sp.Expr] | None,
    variables: Sequence[sp.Symbol],
    default: Sequence[Sequence[int]],
) -> tuple[tuple[int, ...], ...]:
    if order_ideal is None:
        normalized = _sort_exponents(default)
    else:
        exponents: list[tuple[int, ...]] = []
        for item in order_ideal:
            if isinstance(item, (tuple, list)) and all(isinstance(value, int) for value in item):
                exponent = tuple(int(value) for value in item)
            else:
                exponent = _exponent_from_monomial(sp.sympify(item), variables)
            if len(exponent) != len(variables):
                raise BorderBasisError("order ideal exponent has the wrong dimension")
            exponents.append(exponent)
        normalized = _sort_exponents(exponents)
    normalized_set = set(normalized)
    zero = tuple(0 for _ in variables)
    if zero not in normalized_set:
        raise BorderBasisError("order ideal must contain 1")
    for exponent in normalized:
        for index, power in enumerate(exponent):
            for lower_power in range(power):
                divisor = list(exponent)
                divisor[index] = lower_power
                if tuple(divisor) not in normalized_set:
                    raise BorderBasisError("order ideal is not closed under divisibility")
    return normalized


def _map_permuted_exponents_to_original(
    exponents: Sequence[Sequence[int]],
    permuted_variables: Sequence[sp.Symbol],
    original_variables: Sequence[sp.Symbol],
) -> tuple[tuple[int, ...], ...]:
    positions = {var: idx for idx, var in enumerate(original_variables)}
    mapped: list[tuple[int, ...]] = []
    for exponent in exponents:
        out = [0] * len(original_variables)
        for idx, power in enumerate(exponent):
            out[positions[permuted_variables[idx]]] = int(power)
        mapped.append(tuple(out))
    return tuple(mapped)


def _total_degree(poly: sp.Expr, variables: Sequence[sp.Symbol]) -> int:
    return int(sp.Poly(poly, *variables, domain=sp.QQ).total_degree())


def _monomial_exponents_upto(variable_count: int, max_degree: int) -> tuple[tuple[int, ...], ...]:
    exponents: list[tuple[int, ...]] = []
    for degree in range(max_degree + 1):
        exponents.extend(_degree_exponents(variable_count, degree))
    return tuple(exponents)


def _macaulay_column_order(variable_count: int, max_degree: int) -> tuple[tuple[int, ...], ...]:
    """Return monomials ordered so row reduction pivots prefer large terms."""

    exponents = _monomial_exponents_upto(variable_count, max_degree)
    # Row reduction treats earlier columns as pivot candidates. Put higher
    # degree monomials first, and among equal degrees pivot monomials involving
    # later user variables before earlier ones. This leaves quotient bases in
    # earlier variables when possible, matching the public Groebner-derived
    # border-basis preference, e.g. ``(1, x)`` rather than ``(1, y)`` for
    # ``[x**2 - 1, y - x]`` with variables ``[x, y]``.
    return tuple(sorted(exponents, key=lambda exp: (sum(exp), tuple(reversed(exp))), reverse=True))


def _poly_row(
    poly: sp.Expr, variables: Sequence[sp.Symbol], columns: Sequence[Sequence[int]]
) -> list[sp.Expr]:
    p = sp.Poly(sp.expand(poly), *variables, domain=sp.QQ)
    return [p.coeff_monomial(_monomial_from_exponent(variables, exp)) for exp in columns]


def _macaulay_rows(
    polynomials: Sequence[sp.Expr],
    variables: Sequence[sp.Symbol],
    max_degree: int,
    columns: Sequence[Sequence[int]],
) -> list[list[sp.Expr]]:
    """Build exact Macaulay rows from multiples whose total degree is bounded."""

    rows: list[list[sp.Expr]] = []
    variable_count = len(variables)
    for poly in polynomials:
        poly_degree = _total_degree(poly, variables)
        if poly_degree < 0 or poly_degree > max_degree:
            continue
        for multiplier_exp in _monomial_exponents_upto(variable_count, max_degree - poly_degree):
            multiplier = _monomial_from_exponent(variables, multiplier_exp)
            row = _poly_row(sp.expand(multiplier * poly), variables, columns)
            if any(value != 0 for value in row):
                rows.append(row)
    return rows


def _is_order_ideal(exponents: Sequence[Sequence[int]]) -> bool:
    order_set = {tuple(exp) for exp in exponents}
    if not order_set:
        return False
    zero = tuple(0 for _ in next(iter(order_set)))
    if zero not in order_set:
        return False
    for exponent in order_set:
        for index, power in enumerate(exponent):
            if power <= 0:
                continue
            divisor = list(exponent)
            divisor[index] -= 1
            if tuple(divisor) not in order_set:
                return False
    return True
