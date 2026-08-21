from __future__ import annotations

from collections.abc import Iterable, Sequence

import sympy as sp

from ..cache import CACHE, expr_key
from .quotient import (
    _as_rational_polynomial,
    _coefficient_vector,
    _leading_exponent_grevlex,
    _multiplication_matrix,
    _multiplication_tensor,
    _normal_form,
    _select_separating_linear_form,
    _standard_exponents,
)
from .representation import RationalUnivariateError, RationalUnivariateRepresentation


def compute_rational_univariate_representation(
    polynomials: Iterable[sp.Expr],
    variables: Sequence[sp.Symbol],
    parameter: sp.Symbol | None = None,
    *,
    max_separating_attempts: int = 64,
) -> RationalUnivariateRepresentation:
    """Compute a RUR for a rational zero-dimensional polynomial system.

    Construct a rational univariate representation from the quotient algebra.

    Inputs must be rational polynomial equalities and the standard monomial
    basis must be finite.
    """

    variable_tuple = tuple(variables)
    raw_polynomials = tuple(sp.sympify(poly) for poly in polynomials)
    if not variable_tuple:
        raise RationalUnivariateError("at least one variable is required")
    if parameter is None:
        used_names = {str(symbol) for symbol in variable_tuple}
        name = "_rur_t"
        while name in used_names:
            name = "_" + name
        parameter = sp.Symbol(name)
    if parameter in variable_tuple:
        raise RationalUnivariateError("parameter must be distinct from system variables")

    polys = [_as_rational_polynomial(poly, variable_tuple).as_expr() for poly in raw_polynomials]
    cache_key = (
        tuple(expr_key(poly) for poly in polys),
        tuple(sp.srepr(v) for v in variable_tuple),
        sp.srepr(parameter),
        int(max_separating_attempts),
    )
    cached = CACHE.rur.get(cache_key)
    if cached is not None:
        CACHE.stats.rur_hits += 1
        return cached  # type: ignore[return-value]
    CACHE.stats.rur_misses += 1
    if len(polys) < len(variable_tuple):
        raise RationalUnivariateError(
            "at least as many equations as variables are required for rational univariate solving"
        )

    groebner_basis = sp.groebner(polys, *variable_tuple, order="grevlex", domain=sp.QQ)
    if groebner_basis.polys == [sp.Poly(1, *variable_tuple, domain=sp.QQ)]:
        result = RationalUnivariateRepresentation(
            variables=variable_tuple,
            parameter=parameter,
            defining_polynomial=sp.Poly(1, parameter, domain=sp.QQ),
            coordinate_denominator=sp.Poly(1, parameter, domain=sp.QQ),
            coordinate_numerators=tuple(
                sp.Poly(0, parameter, domain=sp.QQ) for _ in variable_tuple
            ),
            separating_linear_form=sp.Integer(0),
            standard_exponents=tuple(),
            quotient_dimension=0,
            geometric_solution_count=0,
        )
        CACHE.rur.put(cache_key, result)
        return result

    leading_exponents = [_leading_exponent_grevlex(poly) for poly in groebner_basis.polys]
    basis_exponents = _standard_exponents(leading_exponents, len(variable_tuple))
    tensor = _multiplication_tensor(groebner_basis, variable_tuple, basis_exponents)
    linear_form, defining_poly, denominator_poly, trace_vector, powers, geometric_count = (
        _select_separating_linear_form(
            groebner_basis,
            variable_tuple,
            parameter,
            basis_exponents,
            tensor,
            max_attempts=max_separating_attempts,
        )
    )

    degree = defining_poly.degree()
    coeffs_ascending = list(reversed(defining_poly.all_coeffs()))
    horner_vectors: list[sp.Matrix] = []
    for power_index in range(degree):
        vector = sp.zeros(len(basis_exponents), 1)
        for shift in range(power_index + 1):
            coeff = coeffs_ascending[degree - shift]
            vector += coeff * powers[power_index - shift]
        horner_vectors.append(sp.simplify(vector))

    coordinate_numerators: list[sp.Poly] = []
    for variable in variable_tuple:
        remainder = _normal_form(groebner_basis, variable)
        variable_vector = _coefficient_vector(remainder, variable_tuple, basis_exponents)
        variable_mult = _multiplication_matrix(variable_vector, tensor)
        trace_variable_products = variable_mult.T * trace_vector
        numerator = sp.Integer(0)
        for power_index, horner_vector in enumerate(horner_vectors):
            numerator += (horner_vector.T * trace_variable_products)[0] * parameter ** (
                degree - power_index - 1
            )
        coordinate_numerators.append(sp.Poly(sp.expand(numerator), parameter, domain=sp.QQ))

    result = RationalUnivariateRepresentation(
        variables=variable_tuple,
        parameter=parameter,
        defining_polynomial=defining_poly,
        coordinate_denominator=denominator_poly,
        coordinate_numerators=tuple(coordinate_numerators),
        separating_linear_form=linear_form,
        standard_exponents=tuple(tuple(e) for e in basis_exponents),
        quotient_dimension=len(basis_exponents),
        geometric_solution_count=geometric_count,
    )
    CACHE.rur.put(cache_key, result)
    return result
