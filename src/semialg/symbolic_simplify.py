from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field

import sympy as sp
from sympy.logic.boolalg import Boolean

from .decision import implies, is_satisfiable, is_tautology
from .normalization import normalize_formula as _normalize_formula
from .normalization import normalize_variables as _normalize_variables

FormulaLike = sp.Expr | Boolean | bool

_RECOVERABLE_ERRORS = (
    ArithmeticError,
    TypeError,
    ValueError,
    NotImplementedError,
    RuntimeError,
    sp.PolynomialError,
)


@dataclass(frozen=True)
class BooleanSimplificationResult:
    """Structured result returned by ``simplify_boole(..., return_result=True)``."""

    formula: sp.Expr
    original: sp.Expr
    variables: tuple[sp.Symbol, ...]
    assumptions: sp.Expr = sp.true
    method: str = "semantic_boolean_simplification"
    removed_redundant: tuple[sp.Expr, ...] = ()
    diagnostics: dict[str, object] = field(default_factory=dict)

    def __bool__(self) -> bool:
        return self.formula is not sp.false and self.formula != sp.false


@dataclass(frozen=True)
class PiecewiseSimplificationResult:
    """Structured result returned by ``simplify_piecewise(..., return_result=True)``."""

    expression: sp.Expr
    original: sp.Expr
    variables: tuple[sp.Symbol, ...]
    assumptions: sp.Expr = sp.true
    removed_unreachable: tuple[tuple[sp.Expr, sp.Expr], ...] = ()
    merged_branches: int = 0
    simplified_branch_values: int = 0
    diagnostics: dict[str, object] = field(default_factory=dict)

    def __bool__(self) -> bool:
        return self.expression is not sp.nan


def _safe_simplify_logic(expr: sp.Expr) -> sp.Expr:
    """Simplify logic without routing Boolean formulas through scalar simplifiers."""

    if isinstance(expr, Boolean):
        try:
            simplified = sp.simplify_logic(expr, form="dnf")
            return simplified if simplified is not None else expr
        except _RECOVERABLE_ERRORS:
            return expr
    try:
        return sp.simplify(expr)
    except _RECOVERABLE_ERRORS:
        return expr


def _semantic_unsat(expr: sp.Expr, variables: Sequence[sp.Symbol], strategy: str | None) -> bool:
    if expr is sp.false or expr == sp.false:
        return True
    if expr is sp.true or expr == sp.true:
        return False
    try:
        return not is_satisfiable(expr, variables, strategy=strategy)
    except _RECOVERABLE_ERRORS:
        return False


def _semantic_tautology(
    expr: sp.Expr, variables: Sequence[sp.Symbol], strategy: str | None
) -> bool:
    if expr is sp.true or expr == sp.true:
        return True
    if expr is sp.false or expr == sp.false:
        return False
    try:
        return is_tautology(expr, variables, strategy=strategy)
    except _RECOVERABLE_ERRORS:
        return False


def _semantic_implies(
    lhs: sp.Expr,
    rhs: sp.Expr,
    variables: Sequence[sp.Symbol],
    strategy: str | None,
) -> bool:
    if lhs == rhs:
        return True
    try:
        return implies(lhs, rhs, variables, strategy=strategy)
    except _RECOVERABLE_ERRORS:
        return False


def _complement_relational(expr: sp.Expr) -> sp.Expr | None:
    if isinstance(expr, sp.StrictGreaterThan):
        return expr.lhs <= expr.rhs
    if isinstance(expr, sp.StrictLessThan):
        return expr.lhs >= expr.rhs
    if isinstance(expr, sp.GreaterThan):
        return expr.lhs < expr.rhs
    if isinstance(expr, sp.LessThan):
        return expr.lhs > expr.rhs
    if isinstance(expr, sp.Equality):
        return sp.Ne(expr.lhs, expr.rhs)
    if isinstance(expr, sp.Unequality):
        return sp.Eq(expr.lhs, expr.rhs)
    return None


def _deduplicate(args: Iterable[sp.Expr]) -> list[sp.Expr]:
    out: list[sp.Expr] = []
    seen: set[str] = set()
    for arg in args:
        key = sp.sstr(arg)
        if key not in seen:
            seen.add(key)
            out.append(arg)
    return out


def _normalize_polynomial_side(
    expr: sp.Expr, variables: Sequence[sp.Symbol], *, normalize_sign: bool = True
) -> sp.Expr:
    """Return a primitive, sign-normalized polynomial expression when possible."""

    try:
        poly = sp.Poly(sp.expand(expr), *variables, domain=sp.QQ)
    except _RECOVERABLE_ERRORS:
        try:
            return sp.factor(sp.cancel(expr))
        except _RECOVERABLE_ERRORS:
            return expr
    if poly.is_zero:
        return sp.Integer(0)
    primitive = poly.primitive()[1]
    if normalize_sign and primitive.LC() < 0:
        primitive = -primitive
    return primitive.as_expr()


def _canonical_relational(atom: sp.Expr, variables: Sequence[sp.Symbol]) -> sp.Expr:
    """Canonicalize simple polynomial relations to ``p rel 0``.

    Inequalities are normalized to use ``>=``/``>`` where possible, with
    rational content removed and leading polynomial sign made positive. This
    makes duplicate detection and implication-based absorption much more stable
    without using scalar simplification on Boolean formulas.
    """

    if not isinstance(atom, sp.core.relational.Relational):
        return atom
    lhs, rhs = atom.lhs, atom.rhs
    if isinstance(atom, sp.GreaterThan):
        rel = ">="
        diff = lhs - rhs
    elif isinstance(atom, sp.StrictGreaterThan):
        rel = ">"
        diff = lhs - rhs
    elif isinstance(atom, sp.LessThan):
        rel = ">="
        diff = rhs - lhs
    elif isinstance(atom, sp.StrictLessThan):
        rel = ">"
        diff = rhs - lhs
    elif isinstance(atom, (sp.Equality, sp.Unequality)):
        rel = "==" if isinstance(atom, sp.Equality) else "!="
        diff = lhs - rhs
    else:
        return atom

    expr = _normalize_polynomial_side(diff, variables, normalize_sign=rel in {"==", "!="})
    if expr == 0:
        if rel in {">=", "=="}:
            return sp.true
        if rel in {">", "!="}:
            return sp.false
    if rel == ">=":
        return expr >= 0
    if rel == ">":
        return expr > 0
    if rel == "==":
        return sp.Eq(expr, 0)
    return sp.Ne(expr, 0)


def _set_to_formula(set_expr: sp.Set, var: sp.Symbol) -> sp.Expr | None:
    """Convert a one-dimensional real set into an equivalent Boolean formula."""

    if set_expr is sp.S.EmptySet or set_expr == sp.S.EmptySet:
        return sp.false
    if set_expr is sp.S.Reals or set_expr == sp.S.Reals:
        return sp.true
    if isinstance(set_expr, sp.ConditionSet):
        return None
    pieces = set_expr.args if isinstance(set_expr, sp.Union) else (set_expr,)
    formulas: list[sp.Expr] = []
    for piece in pieces:
        if isinstance(piece, sp.Interval):
            clauses: list[sp.Expr] = []
            if piece.start != -sp.oo:
                clauses.append(var > piece.start if piece.left_open else var >= piece.start)
            if piece.end != sp.oo:
                clauses.append(var < piece.end if piece.right_open else var <= piece.end)
            formulas.append(sp.And(*clauses) if clauses else sp.true)
        elif isinstance(piece, sp.FiniteSet):
            formulas.extend(sp.Eq(var, point) for point in sorted(piece, key=sp.default_sort_key))
        else:
            return None
    return sp.Or(*formulas) if len(formulas) > 1 else (formulas[0] if formulas else sp.false)


def _try_univariate_interval_simplify(
    expr: sp.Expr, variables: Sequence[sp.Symbol]
) -> sp.Expr | None:
    """Use SymPy's exact one-dimensional set logic for univariate formulas."""

    if len(variables) != 1:
        return None
    var = variables[0]
    try:
        reduced = sp.reduce_inequalities(
            list(expr.args) if isinstance(expr, sp.And) else [expr], var
        )
    except _RECOVERABLE_ERRORS:
        reduced = expr
    try:
        set_expr = reduced.as_set()
    except _RECOVERABLE_ERRORS:
        return None
    formula = _set_to_formula(set_expr, var)
    if formula is None:
        return None
    return _safe_simplify_logic(formula)


def _simplify_and(
    args: Sequence[sp.Expr], variables: Sequence[sp.Symbol], strategy: str | None
) -> sp.Expr:
    pieces: list[sp.Expr] = []
    for arg in args:
        simplified = _simplify_boolean_rec(arg, variables, strategy)
        if simplified is sp.false or simplified == sp.false:
            return sp.false
        if simplified is sp.true or simplified == sp.true:
            continue
        if isinstance(simplified, sp.And):
            pieces.extend(simplified.args)
        else:
            pieces.append(simplified)
    pieces = _deduplicate(pieces)
    if not pieces:
        return sp.true
    candidate = sp.And(*pieces)
    if _semantic_unsat(candidate, variables, strategy):
        return sp.false

    kept = pieces[:]
    changed = True
    while changed and len(kept) > 1:
        changed = False
        for arg in kept[:]:
            others = [piece for piece in kept if piece is not arg]
            premise = sp.And(*others) if others else sp.true
            if _semantic_implies(premise, arg, variables, strategy):
                kept.remove(arg)
                changed = True
                break
    return _safe_simplify_logic(sp.And(*kept) if kept else sp.true)


def _simplify_or(
    args: Sequence[sp.Expr], variables: Sequence[sp.Symbol], strategy: str | None
) -> sp.Expr:
    pieces: list[sp.Expr] = []
    for arg in args:
        simplified = _simplify_boolean_rec(arg, variables, strategy)
        if simplified is sp.true or simplified == sp.true:
            return sp.true
        if simplified is sp.false or simplified == sp.false:
            continue
        if isinstance(simplified, sp.Or):
            pieces.extend(simplified.args)
        else:
            pieces.append(simplified)
    pieces = _deduplicate(pieces)
    if not pieces:
        return sp.false
    candidate = sp.Or(*pieces)
    if _semantic_tautology(candidate, variables, strategy):
        return sp.true

    kept = pieces[:]
    changed = True
    while changed and len(kept) > 1:
        changed = False
        for arg in kept[:]:
            others = [piece for piece in kept if piece is not arg]
            cover = sp.Or(*others) if others else sp.false
            if _semantic_implies(arg, cover, variables, strategy):
                kept.remove(arg)
                changed = True
                break
    return _safe_simplify_logic(sp.Or(*kept) if kept else sp.false)


def _simplify_boolean_rec(
    expr: sp.Expr, variables: Sequence[sp.Symbol], strategy: str | None
) -> sp.Expr:
    expr = _normalize_formula(expr)
    if expr is sp.true or expr == sp.true or expr is sp.false or expr == sp.false:
        return expr
    if isinstance(expr, sp.And):
        return _simplify_and(expr.args, variables, strategy)
    if isinstance(expr, sp.Or):
        return _simplify_or(expr.args, variables, strategy)
    if isinstance(expr, sp.Not):
        inner = _simplify_boolean_rec(expr.args[0], variables, strategy)
        if inner is sp.true or inner == sp.true:
            return sp.false
        if inner is sp.false or inner == sp.false:
            return sp.true
        complemented = _complement_relational(inner)
        if complemented is not None:
            return complemented
        return _safe_simplify_logic(sp.Not(inner))
    simplified = _canonical_relational(expr, variables)
    simplified = _safe_simplify_logic(simplified)
    if _semantic_unsat(simplified, variables, strategy):
        return sp.false
    if _semantic_tautology(simplified, variables, strategy):
        return sp.true
    return simplified


def simplify_boole(
    expr: FormulaLike | Iterable[FormulaLike],
    variables: Sequence[sp.Symbol | str] | None = None,
    *,
    assumptions: FormulaLike | Iterable[FormulaLike] | None = None,
    strategy: str | None = None,
    form: str = "auto",
    semantic: bool = True,
    return_result: bool = False,
) -> sp.Expr | BooleanSimplificationResult:
    """Simplify a semialgebraic Boolean formula over the real numbers.

    This is a conservative public simplifier. It first performs ordinary SymPy
    Boolean simplification, then uses the CAD/QE-backed decision wrappers to
    remove impossible conjuncts/disjuncts and implication-redundant clauses.
    When assumptions are supplied, simplification is performed on the formula
    restricted to those assumptions and the returned expression is a formula in
    the original variable space.
    """

    formula = _normalize_formula(expr)
    original = formula
    if assumptions is not None:
        assumption_expr = _normalize_formula(assumptions)
        universe_expr = sp.And(assumption_expr, formula)
    else:
        assumption_expr = sp.true
        universe_expr = formula
    vars_ = _normalize_variables(variables, universe_expr)
    diagnostics: dict[str, object] = {"form": form, "semantic": semantic}

    def finish(
        value: sp.Expr, method: str = "semantic_boolean_simplification"
    ) -> sp.Expr | BooleanSimplificationResult:
        value = _safe_simplify_logic(value)
        if return_result:
            return BooleanSimplificationResult(
                value, original, vars_, assumption_expr, method=method, diagnostics=diagnostics
            )
        return value

    if form not in {"auto", "dnf", "cnf", "intervals", "cad"}:
        raise ValueError("form must be one of 'auto', 'dnf', 'cnf', 'intervals', or 'cad'")

    if assumptions is not None and semantic:
        if _semantic_unsat(assumption_expr, vars_, strategy):
            diagnostics["reason"] = "assumptions_inconsistent"
            return finish(sp.false, "assumption_contradiction")
        if _semantic_implies(assumption_expr, formula, vars_, strategy):
            diagnostics["reason"] = "formula_implied_by_assumptions"
            return finish(sp.true, "assumption_implication")
        if _semantic_unsat(sp.And(assumption_expr, formula), vars_, strategy):
            diagnostics["reason"] = "formula_inconsistent_with_assumptions"
            return finish(sp.false, "assumption_contradiction")

    simplified = _simplify_boolean_rec(formula, vars_, strategy if semantic else None)

    if form in {"auto", "intervals"} and len(vars_) == 1:
        interval_formula = _try_univariate_interval_simplify(simplified, vars_)
        if interval_formula is not None:
            simplified = interval_formula
            diagnostics["used_univariate_interval_simplification"] = True

    if assumptions is not None and simplified not in (sp.true, sp.false) and semantic:
        # Remove clauses that become redundant in the assumed region.
        restricted = _simplify_boolean_rec(sp.And(assumption_expr, simplified), vars_, strategy)
        if isinstance(restricted, sp.And):
            pieces = [
                arg
                for arg in restricted.args
                if not _semantic_implies(assumption_expr, arg, vars_, strategy)
            ]
            simplified = _simplify_boolean_rec(
                sp.And(*pieces) if pieces else sp.true, vars_, strategy
            )

    if form == "cnf":
        try:
            simplified = sp.simplify_logic(simplified, form="cnf")
        except _RECOVERABLE_ERRORS:
            pass
    elif form == "dnf":
        try:
            simplified = sp.simplify_logic(simplified, form="dnf")
        except _RECOVERABLE_ERRORS:
            pass

    return finish(simplified)


def _simplify_value(
    value: sp.Expr,
    variables: Sequence[sp.Symbol],
    strategy: str | None,
    assumptions: sp.Expr | None = None,
) -> sp.Expr:
    if isinstance(value, sp.Piecewise):
        return simplify_piecewise(value, variables, assumptions=assumptions, strategy=strategy)
    try:
        from .reasoning import simplify_under_assumptions

        if assumptions is not None:
            return simplify_under_assumptions(value, assumptions, variables, strategy=strategy)
    except _RECOVERABLE_ERRORS:
        pass
    try:
        return sp.simplify(value)
    except _RECOVERABLE_ERRORS:
        return value


def _merge_piecewise_branches(
    branches: list[tuple[sp.Expr, sp.Expr]], variables: Sequence[sp.Symbol], strategy: str | None
) -> list[tuple[sp.Expr, sp.Expr]]:
    merged: list[tuple[sp.Expr, sp.Expr]] = []
    for value, condition in branches:
        if merged and sp.simplify(merged[-1][0] - value) == 0:
            prev_value, prev_condition = merged[-1]
            merged[-1] = (
                prev_value,
                simplify_boole(sp.Or(prev_condition, condition), variables, strategy=strategy),
            )
        else:
            merged.append((value, condition))
    return merged


def simplify_piecewise(
    expr: sp.Expr,
    variables: Sequence[sp.Symbol | str] | None = None,
    *,
    assumptions: FormulaLike | Iterable[FormulaLike] | None = None,
    strategy: str | None = None,
    return_result: bool = False,
) -> sp.Expr | PiecewiseSimplificationResult:
    """Simplify a Piecewise expression using semialgebraic branch conditions.

    Branches whose effective condition is unsatisfiable are removed. A branch
    becomes unconditional when its effective condition is a tautology. Adjacent
    branches with equal values are merged by disjoining their conditions.
    """

    original_expr = sp.sympify(expr)
    if not isinstance(original_expr, sp.Piecewise):
        simplified = _simplify_value(
            original_expr,
            (),
            strategy,
            _normalize_formula(assumptions) if assumptions is not None else None,
        )
        result = PiecewiseSimplificationResult(
            simplified,
            original_expr,
            (),
            _normalize_formula(assumptions) if assumptions is not None else sp.true,
        )
        return result if return_result else simplified
    expr = original_expr
    removed_unreachable: list[tuple[sp.Expr, sp.Expr]] = []
    value_simplifications = 0

    all_conditions = [cond for _, cond in expr.args]
    variable_expr = (
        sp.And(*[cond for cond in all_conditions if cond is not True and cond != sp.true])
        if all_conditions
        else sp.true
    )
    if assumptions is not None:
        variable_expr = sp.And(variable_expr, _normalize_formula(assumptions))
    vars_ = _normalize_variables(variables, variable_expr)
    assumption_expr = _normalize_formula(assumptions) if assumptions is not None else sp.true

    covered = sp.false
    branches: list[tuple[sp.Expr, sp.Expr]] = []
    for value, condition in expr.args:
        condition_expr = sp.true if condition is True else _normalize_formula(condition)
        if _semantic_tautology(covered, vars_, strategy):
            removed_unreachable.append((value, condition_expr))
            continue
        effective = simplify_boole(
            sp.And(assumption_expr, condition_expr, sp.Not(covered)),
            vars_,
            strategy=strategy,
        )
        if _semantic_unsat(effective, vars_, strategy):
            removed_unreachable.append((value, condition_expr))
            covered = simplify_boole(sp.Or(covered, condition_expr), vars_, strategy=strategy)
            continue

        new_covered = simplify_boole(sp.Or(covered, condition_expr), vars_, strategy=strategy)
        if (
            condition_expr is sp.true
            or condition_expr == sp.true
            or _semantic_tautology(new_covered, vars_, strategy)
        ):
            output_condition = sp.true
        else:
            output_condition = simplify_boole(
                sp.And(condition_expr, sp.Not(covered)), vars_, strategy=strategy
            )
        if assumptions is not None:
            if _semantic_implies(assumption_expr, output_condition, vars_, strategy):
                output_condition = sp.true
            elif _semantic_unsat(sp.And(assumption_expr, output_condition), vars_, strategy):
                removed_unreachable.append((value, condition_expr))
                covered = simplify_boole(sp.Or(covered, condition_expr), vars_, strategy=strategy)
                continue

        branch_assumptions = simplify_boole(
            sp.And(assumption_expr, effective), vars_, strategy=strategy
        )
        simplified_value = _simplify_value(value, vars_, strategy, branch_assumptions)
        if simplified_value != value:
            value_simplifications += 1
        branches.append((simplified_value, output_condition))
        covered = new_covered

    if not branches:
        result = PiecewiseSimplificationResult(
            sp.nan,
            original_expr,
            vars_,
            assumption_expr,
            tuple(removed_unreachable),
            diagnostics={"reason": "no_reachable_branches"},
        )
        return result if return_result else sp.nan
    before_merge = len(branches)
    branches = _merge_piecewise_branches(branches, vars_, strategy)
    merged_count = before_merge - len(branches)
    if assumptions is not None:
        for candidate, _ in reversed(branches):
            if all(
                _semantic_implies(
                    sp.And(assumption_expr, condition),
                    sp.Eq(value, candidate),
                    vars_,
                    strategy,
                )
                for value, condition in branches
            ):
                result = PiecewiseSimplificationResult(
                    candidate,
                    original_expr,
                    vars_,
                    assumption_expr,
                    tuple(removed_unreachable),
                    merged_count,
                    value_simplifications,
                )
                return result if return_result else candidate
    if len(branches) == 1 and (branches[0][1] is sp.true or branches[0][1] == sp.true):
        final = branches[0][0]
    else:
        final = sp.Piecewise(*branches, evaluate=False)
    result = PiecewiseSimplificationResult(
        final,
        original_expr,
        vars_,
        assumption_expr,
        tuple(removed_unreachable),
        merged_count,
        value_simplifications,
        {"input_branch_count": len(original_expr.args), "output_branch_count": len(branches)},
    )
    return result if return_result else final


__all__ = [
    "BooleanSimplificationResult",
    "PiecewiseSimplificationResult",
    "simplify_boole",
    "simplify_piecewise",
]
