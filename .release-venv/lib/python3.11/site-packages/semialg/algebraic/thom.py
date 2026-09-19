from __future__ import annotations

from dataclasses import dataclass
from functools import cmp_to_key

import sympy as sp

from .._zero_testing import certified_zero
from ..exact_arithmetic import compare_exact_reals
from ..structural_keys import symbol_identity_key
from .exact_sign import exact_algebraic_sign
from .intervals import RationalInterval
from .roots import rational_intvl_around


def _univariate_poly(poly: sp.Poly | sp.Expr, variable: sp.Symbol | None) -> sp.Poly:
    if isinstance(poly, sp.Poly):
        if poly.is_multivariate:
            if variable is None:
                raise ValueError("a variable is required for multivariate polynomial input")
            return sp.Poly(poly.as_expr(), variable, extension=True)
        return sp.Poly(poly.as_expr(), poly.gens[0], extension=True)
    expr = sp.sympify(poly)
    if variable is None:
        symbols = sorted(expr.free_symbols, key=symbol_identity_key)
        if len(symbols) != 1:
            raise ValueError("a variable is required unless the expression is univariate")
        variable = symbols[0]
    return sp.Poly(expr, variable, extension=True)


def _exact_algebraic_sign(expr: sp.Expr) -> int:
    sign = exact_algebraic_sign(expr)
    if sign is None:
        raise ValueError(f"could not certify exact algebraic sign of {sp.sstr(expr)}")
    return sign


@dataclass(frozen=True)
class ThomEncoding:
    """Thom encoding of one real root of a univariate polynomial.

    ``derivative_signs[i - 1]`` is the exact sign of the ``i``-th derivative
    at ``root``. Together with the defining equation, this sign condition
    uniquely identifies the represented real root.
    """

    polynomial: sp.Poly
    root: sp.Expr
    derivative_signs: tuple[int, ...]
    isolating_interval: RationalInterval

    @property
    def variable(self) -> sp.Symbol:
        return self.polynomial.gens[0]

    @property
    def degree(self) -> int:
        return int(self.polynomial.degree())

    @property
    def multiplicity(self) -> int:
        for order, sign in enumerate(self.derivative_signs, start=1):
            if sign != 0:
                return order
        return self.degree

    def sign_of(self, poly: sp.Poly | sp.Expr) -> int:
        """Return the exact sign of another univariate polynomial at this root."""

        target = _univariate_poly(poly, self.variable)
        if target.degree() < 0:
            return 0
        remainder = target.rem(self.polynomial)
        if remainder.is_zero:
            return 0
        return _exact_algebraic_sign(remainder.as_expr().subs(self.variable, self.root))

    def compare(self, other: ThomEncoding) -> int:
        """Compare two represented real roots exactly."""

        return compare_exact_reals(self.root, other.root)


def thom_encoding(
    poly: sp.Poly | sp.Expr,
    root: sp.Expr,
    variable: sp.Symbol | None = None,
) -> ThomEncoding:
    """Construct the exact Thom encoding of ``root`` for ``poly``."""

    polynomial = _univariate_poly(poly, variable)
    root = sp.sympify(root)
    if certified_zero(polynomial.as_expr().subs(polynomial.gens[0], root)) is not True:
        raise ValueError("root does not satisfy the defining polynomial")
    signs = []
    derivative = polynomial.as_expr()
    var = polynomial.gens[0]
    for _order in range(1, polynomial.degree() + 1):
        derivative = sp.diff(derivative, var)
        signs.append(_exact_algebraic_sign(derivative.subs(var, root)))
    return ThomEncoding(polynomial, root, tuple(signs), rational_intvl_around(root))


def thom_encodings(
    poly: sp.Poly | sp.Expr,
    variable: sp.Symbol | None = None,
) -> tuple[ThomEncoding, ...]:
    """Return Thom encodings of all distinct real roots in increasing order."""

    polynomial = _univariate_poly(poly, variable)
    try:
        roots = tuple(sp.real_roots(polynomial.as_expr()))
    except (NotImplementedError, sp.PolynomialError, ValueError) as exc:
        raise ValueError("could not isolate all exact real roots") from exc
    unique: list[sp.Expr] = []
    for root in roots:
        if not any(compare_exact_reals(root, old) == 0 for old in unique):
            unique.append(root)
    unique.sort(key=cmp_to_key(compare_exact_reals))
    return tuple(thom_encoding(polynomial, root) for root in unique)


__all__ = ["ThomEncoding", "thom_encoding", "thom_encodings"]
