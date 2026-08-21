from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field

import sympy as sp
from sympy.logic.boolalg import Boolean

from .symbol_resolution import normalize_variables as _resolve_variables

FormulaLike = sp.Expr | Boolean | bool


@dataclass(frozen=True)
class DomainNormalizationResult:
    """Result of rewriting domain-sensitive real constraints.

    ``formula`` is a semialgebraic Boolean formula suitable for the existing
    solver. ``domain_constraints`` contains constraints that came only from
    expression domains, while ``rewrites`` records the high-level transforms
    used for diagnostics.
    """

    formula: sp.Expr
    variables: tuple[sp.Symbol, ...]
    domain_constraints: tuple[sp.Expr, ...] = ()
    rewrites: tuple[str, ...] = ()
    auxiliary_symbols: tuple[sp.Symbol, ...] = ()
    diagnostics: dict[str, object] = field(default_factory=dict)


def _as_formula(constraints: FormulaLike | Iterable[FormulaLike]) -> sp.Expr:
    if isinstance(constraints, (list, tuple, set, frozenset)):
        pieces = [sp.sympify(item) for item in constraints]
        return sp.And(*pieces) if pieces else sp.true
    if constraints is True:
        return sp.true
    if constraints is False:
        return sp.false
    return sp.sympify(constraints)


def _as_symbols(
    variables: Sequence[sp.Symbol | str] | None, expr: sp.Expr
) -> tuple[sp.Symbol, ...]:
    return _resolve_variables(variables, context=(expr,), append_context_symbols=True)


def _boolean_map(expr: sp.Expr, fn) -> sp.Expr:
    if expr is sp.true or expr == sp.true or expr is sp.false or expr == sp.false:
        return expr
    if isinstance(expr, sp.And):
        return sp.And(*(_boolean_map(arg, fn) for arg in expr.args))
    if isinstance(expr, sp.Or):
        return sp.Or(*(_boolean_map(arg, fn) for arg in expr.args))
    if isinstance(expr, sp.Not):
        return sp.Not(_boolean_map(expr.args[0], fn))
    if isinstance(expr, sp.Implies):
        return sp.Implies(_boolean_map(expr.args[0], fn), _boolean_map(expr.args[1], fn))
    return fn(expr)


def _domain_constraints_for_expr(expr: sp.Expr) -> tuple[sp.Expr, ...]:
    constraints: list[sp.Expr] = []
    for sub in sp.preorder_traversal(expr):
        if isinstance(sub, sp.Pow):
            base, exponent = sub.as_base_exp()
            if exponent.is_Rational:
                q = int(exponent.q)
                if q % 2 == 0:
                    constraints.append(base >= 0)
        elif sub.func == sp.log and sub.args:
            constraints.append(sub.args[0] > 0)
    # Denominators must be nonzero.  Apply ``together`` only to scalar
    # expression parts; passing Boolean formulas to it triggers deprecated
    # Boolean arithmetic inside SymPy.
    scalar_parts: list[sp.Expr] = []
    if getattr(expr, "is_Relational", False):
        scalar_parts.extend((expr.lhs, expr.rhs))
    elif isinstance(expr, sp.logic.boolalg.Boolean):
        for atom in sp.preorder_traversal(expr):
            if getattr(atom, "is_Relational", False):
                scalar_parts.extend((atom.lhs, atom.rhs))
    else:
        scalar_parts.append(expr)
    for part in scalar_parts:
        try:
            _, den = sp.fraction(sp.together(part))
            if den != 1:
                constraints.append(sp.Ne(den, 0))
        except (TypeError, ValueError, NotImplementedError, sp.SympifyError):
            continue
    return tuple(dict.fromkeys(constraints))


def function_domain(expr: sp.Expr, variables: Sequence[sp.Symbol | str] | None = None) -> sp.Expr:
    """Return a conservative real-domain formula for supported expressions.

    The initial scope is semialgebraic-friendly: rational denominators,
    square/even roots, and logarithm positivity conditions are recognized.
    """

    sym_expr = sp.sympify(expr)
    constraints = _domain_constraints_for_expr(sym_expr)
    return sp.And(*constraints) if constraints else sp.true


def is_real_valued(
    expr: sp.Expr,
    variables: Sequence[sp.Symbol | str] | None = None,
    *,
    assumptions: FormulaLike = True,
) -> bool:
    """Return whether supported domain conditions follow from assumptions."""

    domain = function_domain(expr, variables)
    if domain is sp.true or domain == sp.true:
        return True
    assumption_formula = _as_formula(assumptions)
    vars_ = _as_symbols(variables, sp.And(domain, assumption_formula))
    try:
        from .decision import implies

        return implies(assumption_formula, domain, vars_)
    except Exception:
        try:
            return bool(sp.simplify(sp.Implies(assumption_formula, domain)))
        except Exception:
            return False


def _relation_from_head(head, lhs: sp.Expr, rhs: sp.Expr) -> sp.Expr:
    if head is sp.StrictLessThan:
        return lhs < rhs
    if head is sp.LessThan:
        return lhs <= rhs
    if head is sp.StrictGreaterThan:
        return lhs > rhs
    if head is sp.GreaterThan:
        return lhs >= rhs
    if head is sp.Equality:
        return sp.Eq(lhs, rhs)
    if head is sp.Unequality:
        return sp.Ne(lhs, rhs)
    return head(lhs, rhs)


def _rewrite_rational_relation(rel: sp.Rel) -> sp.Expr | None:
    diff = sp.together(rel.lhs - rel.rhs)
    try:
        num, den = sp.fraction(diff)
    except Exception:
        return None
    if den == 1:
        return None
    num = sp.expand(num)
    den = sp.expand(den)
    if isinstance(rel, sp.StrictGreaterThan):
        return sp.Or(sp.And(num > 0, den > 0), sp.And(num < 0, den < 0))
    if isinstance(rel, sp.GreaterThan):
        return sp.Or(sp.And(num >= 0, den > 0), sp.And(num <= 0, den < 0))
    if isinstance(rel, sp.StrictLessThan):
        return sp.Or(sp.And(num < 0, den > 0), sp.And(num > 0, den < 0))
    if isinstance(rel, sp.LessThan):
        return sp.Or(sp.And(num <= 0, den > 0), sp.And(num >= 0, den < 0))
    if isinstance(rel, sp.Equality):
        return sp.And(sp.Eq(num, 0), sp.Ne(den, 0))
    if isinstance(rel, sp.Unequality):
        return sp.And(sp.Ne(num, 0), sp.Ne(den, 0))
    return None


def _rewrite_abs_relation(rel: sp.Rel) -> sp.Expr | None:
    lhs, rhs = rel.lhs, rel.rhs
    if lhs.func != sp.Abs:
        return None
    arg = lhs.args[0]
    if isinstance(rel, sp.StrictLessThan):
        return sp.And(rhs > 0, arg < rhs, arg > -rhs)
    if isinstance(rel, sp.LessThan):
        return sp.And(rhs >= 0, arg <= rhs, arg >= -rhs)
    if isinstance(rel, sp.StrictGreaterThan):
        return sp.Or(rhs < 0, arg > rhs, arg < -rhs)
    if isinstance(rel, sp.GreaterThan):
        return sp.Or(rhs <= 0, arg >= rhs, arg <= -rhs)
    if isinstance(rel, sp.Equality):
        return sp.And(rhs >= 0, sp.Or(sp.Eq(arg, rhs), sp.Eq(arg, -rhs)))
    if isinstance(rel, sp.Unequality):
        return sp.Or(rhs < 0, sp.And(rhs >= 0, sp.Ne(arg, rhs), sp.Ne(arg, -rhs)))
    return None


def _is_sqrt_expr(expr: sp.Expr) -> bool:
    if not isinstance(expr, sp.Pow):
        return False
    _, exponent = expr.as_base_exp()
    return bool(exponent == sp.Rational(1, 2))


def _rewrite_sqrt_relation(rel: sp.Rel) -> sp.Expr | None:
    lhs, rhs = rel.lhs, rel.rhs
    if not _is_sqrt_expr(lhs):
        return None
    radicand = lhs.base
    domain = radicand >= 0
    if isinstance(rel, sp.StrictLessThan):
        return sp.And(domain, rhs > 0, radicand < rhs**2)
    if isinstance(rel, sp.LessThan):
        return sp.And(domain, rhs >= 0, radicand <= rhs**2)
    if isinstance(rel, sp.StrictGreaterThan):
        return sp.And(domain, sp.Or(rhs < 0, radicand > rhs**2))
    if isinstance(rel, sp.GreaterThan):
        return sp.And(domain, sp.Or(rhs <= 0, radicand >= rhs**2))
    if isinstance(rel, sp.Equality):
        return sp.And(domain, rhs >= 0, sp.Eq(radicand, rhs**2))
    if isinstance(rel, sp.Unequality):
        return sp.And(domain, sp.Or(rhs < 0, rhs >= 0, sp.Ne(radicand, rhs**2)))
    return None


def _rewrite_minmax_relation(rel: sp.Rel) -> sp.Expr | None:
    lhs, rhs = rel.lhs, rel.rhs
    if lhs.func not in {sp.Max, sp.Min}:
        return None
    args = lhs.args
    if lhs.func == sp.Max:
        if isinstance(rel, sp.StrictLessThan):
            return sp.And(*(arg < rhs for arg in args))
        if isinstance(rel, sp.LessThan):
            return sp.And(*(arg <= rhs for arg in args))
        if isinstance(rel, sp.StrictGreaterThan):
            return sp.Or(*(arg > rhs for arg in args))
        if isinstance(rel, sp.GreaterThan):
            return sp.Or(*(arg >= rhs for arg in args))
        if isinstance(rel, sp.Equality):
            return sp.Or(
                *(sp.And(sp.Eq(arg, rhs), *(arg >= other for other in args)) for arg in args)
            )
    if lhs.func == sp.Min:
        if isinstance(rel, sp.StrictLessThan):
            return sp.Or(*(arg < rhs for arg in args))
        if isinstance(rel, sp.LessThan):
            return sp.Or(*(arg <= rhs for arg in args))
        if isinstance(rel, sp.StrictGreaterThan):
            return sp.And(*(arg > rhs for arg in args))
        if isinstance(rel, sp.GreaterThan):
            return sp.And(*(arg >= rhs for arg in args))
        if isinstance(rel, sp.Equality):
            return sp.Or(
                *(sp.And(sp.Eq(arg, rhs), *(arg <= other for other in args)) for arg in args)
            )
    return None


def _rewrite_piecewise_relation(rel: sp.Rel) -> sp.Expr | None:
    lhs, rhs = rel.lhs, rel.rhs
    if not isinstance(lhs, sp.Piecewise):
        return None
    branches: list[sp.Expr] = []
    previous = sp.false
    for expr, cond in lhs.args:
        active = sp.And(sp.Not(previous), cond)
        branch_rel = _relation_from_head(rel.func, expr, rhs)
        branches.append(sp.And(active, _rewrite_relation(branch_rel)))
        previous = sp.Or(previous, cond)
    return sp.Or(*branches) if branches else sp.false


def _rewrite_relation(atom: sp.Expr) -> sp.Expr:
    if not isinstance(atom, sp.core.relational.Relational):
        return atom
    # Reorient simple RHS semialgebraic expressions to the left.
    if (
        atom.lhs.func not in {sp.Abs, sp.Max, sp.Min}
        and not _is_sqrt_expr(atom.lhs)
        and not isinstance(atom.lhs, sp.Piecewise)
    ):
        if (
            atom.rhs.func in {sp.Abs, sp.Max, sp.Min}
            or _is_sqrt_expr(atom.rhs)
            or isinstance(atom.rhs, sp.Piecewise)
        ):
            if isinstance(atom, sp.StrictLessThan):
                atom = atom.rhs > atom.lhs
            elif isinstance(atom, sp.LessThan):
                atom = atom.rhs >= atom.lhs
            elif isinstance(atom, sp.StrictGreaterThan):
                atom = atom.rhs < atom.lhs
            elif isinstance(atom, sp.GreaterThan):
                atom = atom.rhs <= atom.lhs
            elif isinstance(atom, (sp.Equality, sp.Unequality)):
                atom = atom.func(atom.rhs, atom.lhs)
    for rewriter in (
        _rewrite_piecewise_relation,
        _rewrite_abs_relation,
        _rewrite_sqrt_relation,
        _rewrite_minmax_relation,
        _rewrite_rational_relation,
    ):
        rewritten = rewriter(atom)
        if rewritten is not None:
            return _rewrite_formula(rewritten)
    domain = sp.And(
        *(_domain_constraints_for_expr(atom.lhs) + _domain_constraints_for_expr(atom.rhs))
    )
    return sp.And(domain, atom) if domain is not sp.true and domain != sp.true else atom


def _rewrite_formula(expr: sp.Expr) -> sp.Expr:
    return _boolean_map(expr, _rewrite_relation)


def normalize_domain_sensitive_constraints(
    constraints: FormulaLike | Iterable[FormulaLike],
    variables: Sequence[sp.Symbol | str] | None = None,
) -> DomainNormalizationResult:
    """Rewrite supported rational/domain-sensitive constraints.

    This is intentionally conservative: unsupported atoms are preserved, with
    any recognized domain constraints conjoined. The resulting formula remains
    suitable for the existing real semialgebraic solver whenever the input falls
    in the supported semialgebraic subset.
    """

    expr = _as_formula(constraints)
    vars_ = _as_symbols(variables, expr)
    domain_constraints = tuple(dict.fromkeys(_domain_constraints_for_expr(expr)))
    special_heads = (sp.Abs, sp.Min, sp.Max, sp.Piecewise)
    if not domain_constraints and not any(expr.has(head) for head in special_heads):
        return DomainNormalizationResult(
            formula=expr,
            variables=vars_,
            diagnostics={
                "original_formula": sp.sstr(expr),
                "normalized_formula": sp.sstr(expr),
            },
        )
    rewritten = _rewrite_formula(expr)
    try:
        simplified = sp.simplify_logic(rewritten, form="dnf")
    except Exception:
        simplified = sp.simplify(rewritten)
    rewrites: list[str] = []
    if simplified != expr:
        rewrites.append("domain_sensitive_constraints")
    return DomainNormalizationResult(
        formula=simplified,
        variables=vars_,
        domain_constraints=domain_constraints,
        rewrites=tuple(rewrites),
        diagnostics={"original_formula": sp.sstr(expr), "normalized_formula": sp.sstr(simplified)},
    )


__all__ = [
    "DomainNormalizationResult",
    "function_domain",
    "is_real_valued",
    "normalize_domain_sensitive_constraints",
]
