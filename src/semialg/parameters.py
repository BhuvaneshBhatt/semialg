from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field

import sympy as sp
from sympy.logic.boolalg import Boolean

from .formula import parse_formula
from .qe import qe_by_complete_cad
from .root_classification import RootClassificationResult, classify_real_roots

FormulaLike = sp.Expr | Boolean | bool


@dataclass(frozen=True)
class SolvabilityConditionsResult:
    """Parameter conditions under which a semialgebraic system is solvable."""

    formula: sp.Expr
    constraints: sp.Expr
    variables: tuple[sp.Symbol, ...]
    parameters: tuple[sp.Symbol, ...]
    method: str = "complete_cad_qe"
    diagnostics: Mapping[str, object] = field(default_factory=dict)

    def __bool__(self) -> bool:
        return self.formula is not sp.false and self.formula != sp.false


@dataclass(frozen=True)
class RootCountConditionsResult:
    """Parameter-space conditions grouped by real-root count."""

    polynomial: sp.Expr
    variable: sp.Symbol
    parameters: tuple[sp.Symbol, ...]
    conditions_by_count: Mapping[sp.Expr, sp.Expr]
    classification: RootClassificationResult
    method: str
    diagnostics: Mapping[str, object] = field(default_factory=dict)

    def condition_for_count(self, count: int | sp.Expr) -> sp.Expr:
        return self.conditions_by_count.get(sp.sympify(count), sp.false)


def _as_real_symbol(var: sp.Symbol | str) -> sp.Symbol:
    return sp.Symbol(var, real=True) if isinstance(var, str) else var


def _normalize_formula(formula: FormulaLike | Iterable[FormulaLike]) -> sp.Expr:
    if isinstance(formula, (list, tuple, set, frozenset)):
        pieces = [sp.sympify(piece) for piece in formula]
        return sp.And(*pieces) if pieces else sp.true
    if formula is True:
        return sp.true
    if formula is False:
        return sp.false
    if isinstance(formula, (sp.Basic, Boolean)):
        return formula  # type: ignore[return-value]
    return sp.sympify(formula)


def _normalize_symbols(symbols: Sequence[sp.Symbol | str] | None) -> tuple[sp.Symbol, ...]:
    out: list[sp.Symbol] = []
    seen: set[sp.Symbol] = set()
    for item in symbols or ():
        sym = _as_real_symbol(item)
        if sym not in seen:
            out.append(sym)
            seen.add(sym)
    return tuple(out)


def _ordered_free_parameters(
    expr: sp.Expr, variables: Sequence[sp.Symbol], parameters: Sequence[sp.Symbol] | None
) -> tuple[sp.Symbol, ...]:
    if parameters is not None:
        return _normalize_symbols(parameters)
    vset = set(variables)
    return tuple(sorted(expr.free_symbols - vset, key=lambda sym: sym.name))


def _simplify_condition(expr: sp.Expr) -> sp.Expr:
    if expr is sp.true or expr == sp.true:
        return sp.true
    if expr is sp.false or expr == sp.false:
        return sp.false
    try:
        return sp.simplify_logic(sp.simplify(expr), form="dnf")
    except Exception:
        return sp.simplify(expr)


def _root_count_positive_condition(
    poly: sp.Expr, var: sp.Symbol, params: tuple[sp.Symbol, ...]
) -> sp.Expr:
    grouped = root_count_conditions(poly, var, params)
    pieces: list[sp.Expr] = []
    for count, condition in grouped.items():
        if count is sp.oo or count == sp.oo:
            pieces.append(condition)
            continue
        try:
            if int(count) > 0:
                pieces.append(condition)
        except Exception:
            continue
    return _simplify_condition(sp.Or(*pieces)) if pieces else sp.false


def _single_relational_existential_condition(
    expr: sp.Expr, variables: tuple[sp.Symbol, ...], params: tuple[sp.Symbol, ...]
) -> sp.Expr | None:
    """Fast exact conditions for one-variable polynomial atoms.

    These cover common parameter-condition queries without invoking full CAD,
    for example ``exists x. x**2 + a*x + b == 0`` and
    ``exists x. x**2 + a < 0``. Compound and higher-dimensional systems still
    fall back to the complete CAD/QE backend.
    """

    if len(variables) != 1 or not getattr(expr, "is_Relational", False):
        return None
    var = variables[0]
    lhs = sp.expand(expr.lhs - expr.rhs)  # type: ignore[attr-defined]
    try:
        poly = sp.Poly(lhs, var)
    except Exception:
        return None
    if any(sym not in set(params) | {var} for sym in lhs.free_symbols):
        return None

    rel = expr.rel_op  # type: ignore[attr-defined]
    degree = poly.degree()
    if rel == "==":
        return _root_count_positive_condition(poly.as_expr(), var, params)
    if degree > 2:
        return None

    coeffs = [sp.factor(c) for c in poly.all_coeffs()]
    if degree == 0:
        value = coeffs[0]
        if rel == "<":
            return value < 0
        if rel == "<=":
            return value <= 0
        if rel == ">":
            return value > 0
        if rel == ">=":
            return value >= 0
        if rel == "!=":
            return sp.Ne(value, 0)
        return None
    if degree == 1:
        lead, const = coeffs
        if rel in {"<", "<=", ">", ">=", "!="}:
            const_case = _single_relational_existential_condition(
                sp.Rel(const, 0, rel), variables, params
            )
            return _simplify_condition(sp.Or(sp.Ne(lead, 0), sp.And(sp.Eq(lead, 0), const_case)))
        return None

    lead, middle, const = coeffs
    disc = sp.factor(middle**2 - 4 * lead * const)
    if rel == "<":
        return _simplify_condition(
            sp.Or(
                sp.Lt(lead, 0),
                sp.And(sp.Gt(lead, 0), sp.Gt(disc, 0)),
                sp.And(sp.Eq(lead, 0), sp.Ne(middle, 0)),
                sp.And(sp.Eq(lead, 0), sp.Eq(middle, 0), sp.Lt(const, 0)),
            )
        )
    if rel == "<=":
        return _simplify_condition(
            sp.Or(
                sp.Lt(lead, 0),
                sp.And(sp.Gt(lead, 0), sp.Ge(disc, 0)),
                sp.And(sp.Eq(lead, 0), sp.Ne(middle, 0)),
                sp.And(sp.Eq(lead, 0), sp.Eq(middle, 0), sp.Le(const, 0)),
            )
        )
    if rel == ">":
        return _simplify_condition(
            sp.Or(
                sp.Gt(lead, 0),
                sp.And(sp.Lt(lead, 0), sp.Gt(disc, 0)),
                sp.And(sp.Eq(lead, 0), sp.Ne(middle, 0)),
                sp.And(sp.Eq(lead, 0), sp.Eq(middle, 0), sp.Gt(const, 0)),
            )
        )
    if rel == ">=":
        return _simplify_condition(
            sp.Or(
                sp.Gt(lead, 0),
                sp.And(sp.Lt(lead, 0), sp.Ge(disc, 0)),
                sp.And(sp.Eq(lead, 0), sp.Ne(middle, 0)),
                sp.And(sp.Eq(lead, 0), sp.Eq(middle, 0), sp.Ge(const, 0)),
            )
        )
    if rel == "!=":
        return _simplify_condition(sp.Or(sp.Ne(lead, 0), sp.Ne(middle, 0), sp.Ne(const, 0)))
    return None


def solvability_conditions(
    constraints: FormulaLike | Iterable[FormulaLike],
    variables: Sequence[sp.Symbol | str],
    parameters: Sequence[sp.Symbol | str] | None = None,
    *,
    domain: str = "reals",
    return_result: bool = False,
) -> sp.Expr | SolvabilityConditionsResult:
    """Return parameter conditions for real solvability of a constraint system."""

    if domain.lower() not in {"real", "reals", "r", "rr"}:
        raise NotImplementedError("solvability_conditions currently supports only the real domain")
    expr = _normalize_formula(constraints)
    vars_ = _normalize_symbols(variables)
    params = _ordered_free_parameters(
        expr, vars_, _normalize_symbols(parameters) if parameters is not None else None
    )

    if expr is sp.true or expr == sp.true:
        condition = sp.true
        result = SolvabilityConditionsResult(condition, expr, vars_, params, "trivial")
        return result if return_result else condition
    if expr is sp.false or expr == sp.false:
        condition = sp.false
        result = SolvabilityConditionsResult(condition, expr, vars_, params, "trivial")
        return result if return_result else condition

    fast_condition = _single_relational_existential_condition(expr, vars_, params)
    if fast_condition is not None:
        result = SolvabilityConditionsResult(
            fast_condition, expr, vars_, params, "univariate_polynomial_atom", {"fast_path": True}
        )
        return result if return_result else fast_condition

    quantifiers = tuple(("exists", var) for var in vars_)
    all_vars = tuple(
        dict.fromkeys(
            tuple(params)
            + tuple(vars_)
            + tuple(sorted(expr.free_symbols - set(params) - set(vars_), key=lambda s: s.name))
        )
    )
    qe_result = qe_by_complete_cad(
        all_vars, quantifiers, parse_formula(expr), free_variables=params
    )
    condition = _simplify_condition(qe_result.formula)
    result = SolvabilityConditionsResult(
        condition,
        expr,
        vars_,
        params,
        getattr(qe_result, "backend", "complete_cad_qe"),
        {
            "quantified_variables": tuple(sp.sstr(v) for v in vars_),
            "free_parameters": tuple(sp.sstr(p) for p in params),
            "is_sentence": bool(getattr(qe_result, "is_sentence", False)),
        },
    )
    return result if return_result else condition


def root_count_conditions(
    polynomial: sp.Poly | sp.Expr,
    variable: sp.Symbol | str,
    parameters: Sequence[sp.Symbol | str] | None = None,
    *,
    return_result: bool = False,
) -> Mapping[sp.Expr, sp.Expr] | RootCountConditionsResult:
    """Return parameter conditions grouped by distinct real-root count."""

    var = _as_real_symbol(variable)
    expr = polynomial.as_expr() if isinstance(polynomial, sp.Poly) else sp.sympify(polynomial)
    params = (
        _normalize_symbols(parameters)
        if parameters is not None
        else tuple(sorted(expr.free_symbols - {var}, key=lambda sym: sym.name))
    )
    classification = classify_real_roots(expr, var, parameters=params)
    grouped: dict[sp.Expr, list[sp.Expr]] = {}
    for cell in classification.cells:
        grouped.setdefault(sp.sympify(cell.root_count), []).append(cell.condition)
    conditions_by_count: dict[sp.Expr, sp.Expr] = {}
    for count, pieces in grouped.items():
        if not pieces:
            condition = sp.false
        elif len(pieces) == 1:
            condition = pieces[0]
        else:
            condition = sp.Or(*pieces)
        conditions_by_count[count] = _simplify_condition(condition)
    result = RootCountConditionsResult(
        sp.expand(expr),
        var,
        params,
        conditions_by_count,
        classification,
        classification.method,
        {"cell_count": len(classification.cells)},
    )
    return result if return_result else conditions_by_count


__all__ = [
    "RootCountConditionsResult",
    "SolvabilityConditionsResult",
    "root_count_conditions",
    "solvability_conditions",
]
