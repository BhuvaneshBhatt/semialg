"""Shared normalization helpers for public semialgebraic APIs."""

from __future__ import annotations

from collections.abc import Mapping, Sequence

import sympy as sp

from .exact_arithmetic import compare_exact_reals
from .formula import to_sympy
from .symbol_resolution import normalize_variables as resolve_variables
from .symbol_resolution import resolve_symbol


def normalize_formula(condition: object) -> sp.Expr:
    """Return one SymPy formula, conjoining explicit formula collections."""

    if isinstance(condition, (list, tuple, set, frozenset)):
        pieces = [normalize_formula(piece) for piece in condition]
        return sp.And(*pieces) if pieces else sp.true
    if condition is True:
        return sp.true
    if condition is False:
        return sp.false
    if isinstance(condition, (sp.Basic, sp.logic.boolalg.Boolean)):
        return condition  # type: ignore[return-value]
    return to_sympy(condition)  # type: ignore[arg-type]


def normalize_variables(
    variables: Sequence[sp.Symbol | str] | None,
    *context: object,
    append_context_symbols: bool = False,
    exclude: Sequence[sp.Symbol] = (),
) -> tuple[sp.Symbol, ...]:
    """Resolve variables against symbols already present in the problem."""

    return resolve_variables(
        variables,
        context=context,
        append_context_symbols=append_context_symbols,
        exclude=exclude,
    )


def normalize_problem_variables(
    variables: Sequence[sp.Symbol | str] | None,
    *context: object,
    exclude: Sequence[sp.Symbol] = (),
) -> tuple[sp.Symbol, ...]:
    """Resolve explicit variables and append remaining symbols from the problem."""

    return normalize_variables(
        variables,
        *context,
        append_context_symbols=True,
        exclude=exclude,
    )


def normalize_sampling_variables(
    variables: Sequence[sp.Symbol | str] | None,
    *context: object,
) -> tuple[sp.Symbol, ...]:
    """Resolve sampling variables; infer context symbols only when none are explicit."""

    return normalize_variables(
        variables,
        *context,
        append_context_symbols=variables is None,
    )


def normalize_symbol_sequence(
    variables: Sequence[sp.Symbol | str],
) -> tuple[sp.Symbol, ...]:
    """Normalize an explicit ordered symbol sequence without adding context symbols."""

    return normalize_variables(variables)


def conjuncts(expr: sp.Expr) -> tuple[sp.Expr, ...]:
    """Return flattened top-level conjunction atoms in deterministic order."""

    if expr is sp.true or expr == sp.true:
        return ()
    if isinstance(expr, sp.And):
        items: list[sp.Expr] = []
        for arg in expr.args:
            items.extend(conjuncts(arg))
        return tuple(items)
    return (expr,)


def normalize_bounds(
    bounds: Sequence[tuple[sp.Symbol | str, object, object]]
    | Mapping[sp.Symbol | str, tuple[object, object]]
    | None,
    variables: Sequence[sp.Symbol],
) -> dict[sp.Symbol, tuple[sp.Expr, sp.Expr]]:
    """Normalize optional variable bounds while preserving symbol identity."""

    if bounds is None:
        return {}
    if isinstance(bounds, Mapping):
        items = []
        for key, value in bounds.items():
            if (
                not isinstance(value, Sequence)
                or isinstance(value, (str, bytes))
                or len(value) != 2
            ):
                raise ValueError(f"bound for {key!r} must be a (lower, upper) pair")
            items.append((key, value[0], value[1]))
    else:
        items = list(bounds)
        if any(len(item) != 3 for item in items):
            raise ValueError("sequence bounds must contain (variable, lower, upper) triples")
    known = tuple(variables)
    by_name = {var.name: var for var in known}
    result: dict[sp.Symbol, tuple[sp.Expr, sp.Expr]] = {}
    for raw_var, lower, upper in items:
        if isinstance(raw_var, str):
            var = by_name.get(raw_var)
            if var is None:
                resolve_symbol(raw_var, known_symbols=known)
                raise ValueError(f"bound variable {raw_var!r} is not in the variable list")
        else:
            var = raw_var
            if var not in known:
                raise ValueError(f"bound variable {var!r} is not in the variable list")
        if var in result:
            raise ValueError(f"duplicate bound for variable {var!r}")
        lower_expr = sp.sympify(lower)
        upper_expr = sp.sympify(upper)
        try:
            reversed_bounds = compare_exact_reals(lower_expr, upper_expr) > 0
        except (TypeError, ValueError, NotImplementedError):
            reversed_bounds = sp.simplify(lower_expr - upper_expr).is_positive is True
        if reversed_bounds:
            raise ValueError(f"lower bound exceeds upper bound for {var!r}")
        result[var] = (lower_expr, upper_expr)
    return result
