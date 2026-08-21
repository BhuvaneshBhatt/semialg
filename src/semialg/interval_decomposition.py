"""Exact one-dimensional decomposition utilities shared by public APIs."""

from __future__ import annotations

from collections.abc import Mapping
from functools import cmp_to_key

import sympy as sp
from sympy.logic.boolalg import And, BooleanFalse, BooleanTrue, Not, Or

from .exact_arithmetic import compare_exact_reals, exact_truth
from .relations import split_relation


def truth_at(condition: sp.Expr, values: Mapping[sp.Symbol, sp.Expr]) -> bool:
    """Evaluate a specialized formula using exact real arithmetic."""

    return exact_truth(sp.simplify(condition.subs(values)))


def finite_real_roots(poly_expr: sp.Expr, variable: sp.Symbol) -> tuple[sp.Expr, ...]:
    """Return distinct exact real roots of a univariate polynomial."""

    poly = sp.Poly(poly_expr, variable)
    if poly.is_zero:
        return ()
    try:
        roots = sp.real_roots(poly.as_expr())
    except (TypeError, ValueError, NotImplementedError, sp.PolynomialError) as exc:
        raise NotImplementedError("exact real-root isolation failed") from exc
    result: list[sp.Expr] = []
    seen: set[str] = set()
    for root in roots:
        root_expr = sp.re(root) if not isinstance(root, sp.Expr) else root
        key = sp.sstr(root_expr)
        if key not in seen:
            result.append(root_expr)
            seen.add(key)
    return tuple(result)


def relational_polynomials(condition: sp.Expr) -> tuple[sp.Expr, ...]:
    """Collect zero-normalized polynomials from a Boolean relation formula."""

    if condition in (sp.true, sp.false) or isinstance(
        condition,
        (BooleanTrue, BooleanFalse),
    ):
        return ()
    if getattr(condition, "is_Relational", False):
        expr, _ = split_relation(condition)
        return (expr,)
    if isinstance(condition, (And, Or)):
        result: list[sp.Expr] = []
        for arg in condition.args:
            result.extend(relational_polynomials(arg))
        return tuple(result)
    if isinstance(condition, Not):
        return relational_polynomials(condition.args[0])
    raise TypeError(f"unsupported formula expression: {condition!r}")


def sample_between(left: sp.Expr, right: sp.Expr) -> sp.Expr:
    """Choose a simple exact sample strictly between ordered endpoints."""

    if left == -sp.oo and right == sp.oo:
        return sp.Integer(0)
    if left == -sp.oo:
        return sp.simplify(right - 1)
    if right == sp.oo:
        return sp.simplify(left + 1)
    return sp.simplify((left + right) / 2)


def one_dimensional_intervals(
    condition: sp.Expr,
    variable: sp.Symbol,
    bound: tuple[sp.Expr, sp.Expr] | None,
    *,
    extra_symbol_error: str | None = None,
) -> tuple[tuple[sp.Expr, sp.Expr], ...]:
    """Decompose a univariate semialgebraic formula into true open intervals."""

    if condition is sp.false or isinstance(condition, BooleanFalse):
        return ()
    lower, upper = bound if bound is not None else (-sp.oo, sp.oo)
    cuts: list[sp.Expr] = []
    if lower != -sp.oo:
        cuts.append(lower)
    if upper != sp.oo:
        cuts.append(upper)
    for poly in relational_polynomials(condition):
        if variable not in poly.free_symbols and sp.simplify(poly) != 0:
            continue
        if not poly.free_symbols <= {variable}:
            message = extra_symbol_error or (
                "1D formula contains symbols outside the decomposition variable"
            )
            raise ValueError(message)
        cuts.extend(finite_real_roots(poly, variable))

    finite = [value for value in cuts if value not in (-sp.oo, sp.oo)]
    ordered: list[sp.Expr] = []
    seen: set[str] = set()
    for root in sorted(finite, key=cmp_to_key(compare_exact_reals)):
        if lower != -sp.oo and compare_exact_reals(root, lower) < 0:
            continue
        if upper != sp.oo and compare_exact_reals(root, upper) > 0:
            continue
        key = sp.sstr(root)
        if key not in seen:
            ordered.append(root)
            seen.add(key)

    boundaries = [lower, *ordered, upper]
    result: list[tuple[sp.Expr, sp.Expr]] = []
    for left, right in zip(boundaries, boundaries[1:], strict=False):
        if left == right:
            continue
        sample = sample_between(left, right)
        if truth_at(condition, {variable: sample}):
            result.append((left, right))
    return tuple(result)
