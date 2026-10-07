"""Compatibility façade for finite quotient algebra now owned by algroots.

New Semialg code should use :class:`algroots.quotient.QuotientAlgebra` directly when it
needs equality-only quotient structure.  These private names remain temporarily
for internal/source compatibility while carrying no independent algorithms.
"""

from __future__ import annotations

from collections.abc import Sequence

import sympy as sp
from algroots.errors import QuotientAlgebraError
from algroots.quotient import (
    QuotientAlgebra,
    _componentwise_leq,
    _leading_exponent_grevlex,
    _normal_form_for_basis,
    _require_exact_domain,
)
from algroots.quotient import (
    _coefficient_vector as _algroots_coefficient_vector,
)
from algroots.quotient import (
    monomial_from_exponent as _monomial_from_exponent,
)
from algroots.quotient import (
    standard_exponent_count as _algroots_standard_exponent_count,
)
from algroots.quotient import (
    standard_exponents as _algroots_standard_exponents,
)

from .representation import RationalUnivariateError


def _translate(callable_, *args, **kwargs):
    try:
        return callable_(*args, **kwargs)
    except QuotientAlgebraError as exc:
        raise RationalUnivariateError(str(exc)) from exc


def _as_exact_polynomial(
    expr: sp.Expr,
    variables: Sequence[sp.Symbol],
    *,
    domain=None,
) -> sp.Poly:
    try:
        poly = (
            sp.Poly(sp.expand(expr), *variables, domain=domain)
            if domain is not None
            else sp.Poly(sp.expand(expr), *variables, extension=True)
        )
        if poly.total_degree() < 0:
            raise RationalUnivariateError("zero polynomial is not a valid defining equation")
        return poly.set_domain(_require_exact_domain(poly.domain))
    except RationalUnivariateError:
        raise
    except (
        sp.PolynomialError,
        sp.polys.polyerrors.CoercionFailed,
        ValueError,
        TypeError,
        NotImplementedError,
        QuotientAlgebraError,
    ) as exc:
        raise RationalUnivariateError(
            "RUR requires polynomial equations over an exact algebraic coefficient field"
        ) from exc


def _avoids_leading_monomials(candidate, leading_exponents) -> bool:
    return not any(_componentwise_leq(leading, candidate) for leading in leading_exponents)


def _standard_exponent_count(leading_exponents, variable_count: int) -> int:
    return _translate(_algroots_standard_exponent_count, leading_exponents, variable_count)


def _standard_exponents(leading_exponents, variable_count: int):
    return _translate(_algroots_standard_exponents, leading_exponents, variable_count)


def _normal_form(groebner_basis, expr: sp.Expr) -> sp.Expr:
    return _translate(_normal_form_for_basis, groebner_basis, expr)


def _coefficient_vector(expr, variables, basis_exponents, domain=sp.QQ):
    return _translate(
        _algroots_coefficient_vector,
        expr,
        variables,
        basis_exponents,
        _require_exact_domain(domain),
    )


def _variable_multiplication_matrices(groebner_basis, variables, basis_exponents, domain):
    quotient = _translate(QuotientAlgebra.from_groebner_basis, groebner_basis, variables)
    if tuple(tuple(e) for e in basis_exponents) != quotient.standard_exponents:
        raise RationalUnivariateError(
            "requested basis exponents do not match the quotient staircase"
        )
    if _require_exact_domain(domain) != quotient.domain:
        raise RationalUnivariateError(
            "requested coefficient domain does not match the quotient domain"
        )
    return quotient.variable_multiplication_matrices


__all__ = [
    "_as_exact_polynomial",
    "_leading_exponent_grevlex",
    "_componentwise_leq",
    "_avoids_leading_monomials",
    "_standard_exponent_count",
    "_standard_exponents",
    "_monomial_from_exponent",
    "_normal_form",
    "_coefficient_vector",
    "_variable_multiplication_matrices",
]
