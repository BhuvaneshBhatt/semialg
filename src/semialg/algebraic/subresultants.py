from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import sympy as sp


@dataclass(frozen=True)
class SubresultantPRSResult:
    """Exact subresultant polynomial remainder sequence for two univariate polynomials.

    The sequence is normalized as SymPy ``Poly`` objects over the requested
    domain. When SymPy's native subresultant implementation is available, this
    wrapper uses it directly; otherwise it falls back to a fraction-free
    pseudo-remainder sequence. The fallback is intended as a portability guard,
    not as a replacement for SymPy's optimized dense-polynomial routines.
    """

    variable: sp.Symbol
    polynomials: tuple[sp.Poly, ...]
    principal_coefficients: tuple[sp.Expr, ...]
    resultant: sp.Expr
    source: str

    def as_exprs(self) -> tuple[sp.Expr, ...]:
        """Return the PRS entries as SymPy expressions."""

        return tuple(poly.as_expr() for poly in self.polynomials)


def _infer_variable(
    first: sp.Expr | sp.Poly,
    second: sp.Expr | sp.Poly,
    variable: sp.Symbol | None,
) -> sp.Symbol:
    if variable is not None:
        return variable
    if isinstance(first, sp.Poly) and len(first.gens) == 1:
        return first.gens[0]
    if isinstance(second, sp.Poly) and len(second.gens) == 1:
        return second.gens[0]
    symbols = sorted(
        set(getattr(sp.sympify(first), "free_symbols", set()))
        | set(getattr(sp.sympify(second), "free_symbols", set())),
        key=lambda sym: sym.sort_key(),
    )
    if len(symbols) != 1:
        raise ValueError("subresultant_prs needs a univariate variable when it cannot be inferred")
    return symbols[0]


def _as_poly(expr: sp.Expr | sp.Poly, variable: sp.Symbol, domain=None) -> sp.Poly:
    if isinstance(expr, sp.Poly):
        poly = expr
        if poly.gens != (variable,):
            poly = sp.Poly(poly.as_expr(), variable, domain=domain or poly.domain)
        elif domain is not None:
            poly = sp.Poly(poly.as_expr(), variable, domain=domain)
        return poly
    if domain is None:
        return sp.Poly(sp.expand(expr), variable)
    return sp.Poly(sp.expand(expr), variable, domain=domain)


def _primitive_part(poly: sp.Poly) -> sp.Poly:
    if poly.is_zero:
        return poly
    _, primitive = poly.primitive()
    if primitive.LC().could_extract_minus_sign():
        primitive = -primitive
    return primitive


def _native_subresultants(first: sp.Poly, second: sp.Poly) -> tuple[sp.Poly, ...] | None:
    """Call SymPy's native subresultants if available."""

    try:
        subresultants = sp.subresultants(first, second)
    except AttributeError:
        try:
            from sympy.polys.polytools import subresultants
        except Exception:
            return None
        try:
            subresultants = subresultants(first, second)
        except Exception:
            return None
    except Exception:
        return None
    try:
        return tuple(
            sp.Poly(
                item.as_expr() if isinstance(item, sp.Poly) else item,
                first.gens[0],
                domain=first.domain,
            )
            for item in subresultants
        )
    except Exception:
        return tuple(
            sp.Poly(item.as_expr() if isinstance(item, sp.Poly) else item, first.gens[0])
            for item in subresultants
        )


def _pseudo_remainder_fallback(first: sp.Poly, second: sp.Poly) -> tuple[sp.Poly, ...]:
    """Fraction-free PRS fallback used only when native subresultants are absent."""

    if first.is_zero:
        return (second,)
    if second.is_zero:
        return (first,)
    a, b = first, second
    if a.degree() < b.degree():
        a, b = b, a
    sequence = [a, b]
    while not sequence[-1].is_zero:
        prev, curr = sequence[-2], sequence[-1]
        if curr.degree() <= 0:
            break
        rem = prev.prem(curr)
        if rem.is_zero:
            break
        sequence.append(_primitive_part(rem))
    return tuple(sequence)


def _principal_coefficients(polys: Sequence[sp.Poly]) -> tuple[sp.Expr, ...]:
    coefficients: list[sp.Expr] = []
    for poly in polys:
        if poly.is_zero:
            coefficients.append(sp.Integer(0))
            continue
        # For a subresultant S_k, the principal subresultant coefficient is the
        # coefficient of x^k. Native sequences do not need to expose k
        # separately: deg(S_k) <= k, so the leading coefficient is the PSC for
        # non-defective entries and is the useful coefficient for PRS decisions.
        coefficients.append(sp.factor(poly.LC()))
    return tuple(coefficients)


def subresultant_prs(
    first: sp.Expr | sp.Poly,
    second: sp.Expr | sp.Poly,
    variable: sp.Symbol | None = None,
    *,
    domain=None,
    primitive: bool = False,
) -> SubresultantPRSResult:
    """Return an exact subresultant PRS for two univariate polynomials.

    Parameters
    ----------
    first, second:
        Univariate polynomials as SymPy expressions or ``Poly`` objects.
    variable:
        Main variable. It is inferred when both inputs are univariate in the
        same symbol.
    domain:
        Optional SymPy polynomial domain, e.g. ``sp.QQ`` or ``sp.ZZ``.
    primitive:
        If true, primitive-normalize each returned PRS entry. The default keeps
        SymPy's native normalization because that preserves resultant/PSC scale
        conventions.
    """

    var = _infer_variable(first, second, variable)
    left = _as_poly(first, var, domain=domain)
    right = _as_poly(second, var, domain=domain or left.domain)
    if left.is_zero and right.is_zero:
        raise ValueError("subresultant_prs is undefined for two zero polynomials")
    if left.degree() < right.degree():
        left, right = right, left

    native = _native_subresultants(left, right)
    if native is not None:
        polys = native
        source = "sympy.subresultants"
    else:
        polys = _pseudo_remainder_fallback(left, right)
        source = "fraction_free_pseudo_remainder_fallback"
    if primitive:
        polys = tuple(_primitive_part(poly) for poly in polys)
    resultant = sp.factor(sp.resultant(left.as_expr(), right.as_expr(), var))
    return SubresultantPRSResult(
        variable=var,
        polynomials=tuple(polys),
        principal_coefficients=_principal_coefficients(polys),
        resultant=resultant,
        source=source,
    )


def principal_subresultant_coefficients(
    first: sp.Expr | sp.Poly,
    second: sp.Expr | sp.Poly,
    variable: sp.Symbol | None = None,
    *,
    domain=None,
) -> tuple[sp.Expr, ...]:
    """Return leading/principal coefficients from ``subresultant_prs``."""

    return subresultant_prs(first, second, variable, domain=domain).principal_coefficients


__all__ = [
    "SubresultantPRSResult",
    "subresultant_prs",
    "principal_subresultant_coefficients",
]
