"""Shared normalization helpers for public semialgebraic APIs."""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence

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


def normalize_constraints(condition: object | None) -> sp.Expr:
    """Normalize optional constraints, including general iterables, to one formula."""

    if condition is None:
        return sp.true
    if isinstance(condition, Iterable) and not isinstance(
        condition, (sp.Basic, sp.logic.boolalg.Boolean, str, bytes)
    ):
        pieces = tuple(normalize_formula(piece) for piece in condition)
        return sp.And(*pieces) if pieces else sp.true
    return normalize_formula(condition)


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


def normalize_parameters(
    parameters: Sequence[sp.Symbol | str] | None,
    *context: object,
    exclude: Sequence[sp.Symbol] = (),
) -> tuple[sp.Symbol, ...]:
    """Resolve an ordered parameter list against symbols in ``context``.

    This is the parameter counterpart of :func:`normalize_variables`: explicit
    order is preserved, duplicate symbols are removed, same-name symbols with
    different assumptions are rejected as ambiguous, and no contextual symbols
    are appended implicitly.
    """

    return normalize_variables(
        parameters,
        *context,
        append_context_symbols=False,
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


def normalize_point(
    point: Mapping[sp.Symbol | str, object] | Sequence[object],
    variables: Sequence[sp.Symbol],
    *,
    context: Sequence[object] = (),
) -> dict[sp.Symbol, sp.Expr]:
    """Normalize a point while preserving contextual symbol identity.

    Mapping keys may be exact Symbol objects or unambiguous string names.
    Missing, extra, and ambiguous coordinates are rejected.
    """
    vars_ = tuple(variables)
    if isinstance(point, Mapping):
        out: dict[sp.Symbol, sp.Expr] = {}
        for raw_key, raw_value in point.items():
            symbol = resolve_symbol(raw_key, context=context, known_symbols=vars_)
            if symbol not in vars_:
                raise ValueError(f"point coordinate {symbol!r} is not in the variable list")
            if symbol in out:
                raise ValueError(f"duplicate point coordinate for {symbol!r}")
            out[symbol] = sp.sympify(raw_value)
        missing = tuple(v for v in vars_ if v not in out)
        if missing:
            raise ValueError(f"point is missing coordinates for {missing!r}")
        return out
    values = tuple(sp.sympify(value) for value in point)
    if len(values) != len(vars_):
        raise ValueError("point dimension does not match variables")
    return dict(zip(vars_, values, strict=True))


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


def disjuncts(expr: sp.Expr) -> tuple[sp.Expr, ...]:
    """Return top-level disjunction branches without Boolean simplification."""

    if expr is sp.false or expr == sp.false:
        return ()
    if isinstance(expr, sp.Or):
        return tuple(expr.args)
    return (expr,)


def require_conjunction(
    expr: sp.Expr, *, message: str = "formula must be a conjunction"
) -> tuple[sp.Expr, ...]:
    """Return flattened conjunction atoms or reject Boolean branching.

    This is the shared no-allocation-beyond-flattening helper for algorithmic
    paths that intentionally accept only conjunctions.
    """
    if expr is sp.false or expr == sp.false:
        return (sp.false,)
    if isinstance(expr, (sp.Or, sp.Not)):
        raise NotImplementedError(message)
    atoms = conjuncts(expr)
    if all(getattr(atom, "is_Relational", False) for atom in atoms):
        return atoms
    if not atoms:
        return ()
    raise TypeError(f"unsupported formula expression: {expr!r}")


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
    result: dict[sp.Symbol, tuple[sp.Expr, sp.Expr]] = {}
    for raw_var, lower, upper in items:
        if isinstance(raw_var, str):
            var = resolve_symbol(raw_var, known_symbols=known)
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
