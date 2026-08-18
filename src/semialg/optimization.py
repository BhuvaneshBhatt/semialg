from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from itertools import combinations

import sympy as sp
from sympy.core.relational import (
    Equality,
    GreaterThan,
    LessThan,
    StrictGreaterThan,
    StrictLessThan,
    Unequality,
)
from sympy.core.sympify import SympifyError
from sympy.logic.boolalg import And as SymAnd
from sympy.logic.boolalg import BooleanFalse, BooleanTrue
from sympy.logic.boolalg import Not as SymNot
from sympy.logic.boolalg import Or as SymOr
from sympy.polys.polyerrors import PolynomialError

from .formula import parse_formula, to_sympy

_EXPECTED_ERRORS = (
    TypeError,
    ValueError,
    ArithmeticError,
    NotImplementedError,
    SympifyError,
    PolynomialError,
)


FormulaLike = sp.Expr | sp.logic.boolalg.Boolean | bool


def qe_by_complete_cad(*args, **kwargs):
    """Import the CAD backend only when a range computation needs it."""

    from .qe import qe_by_complete_cad as impl

    return impl(*args, **kwargs)


def satisfies_formula(*args, **kwargs):
    """Import instance checking only for optimization candidate validation."""

    from .instances.real_fallbacks import satisfies_formula as impl

    return impl(*args, **kwargs)


def _finite_real_roots(*args, **kwargs):
    from .measure import _finite_real_roots as impl

    return impl(*args, **kwargs)


def _one_dimensional_intervals(*args, **kwargs):
    from .measure import _one_dimensional_intervals as impl

    return impl(*args, **kwargs)


def _relational_polynomials(*args, **kwargs):
    from .measure import _relational_polynomials as impl

    return impl(*args, **kwargs)


@dataclass(frozen=True)
class FunctionRangeResult:
    """Exact range summary for a supported semialgebraic image problem.

    ``formula`` is the primary answer: a quantifier-free condition on
    ``value_symbol`` describing the range. The bound fields are populated
    when the range has an interval summary and may be ``None`` for disconnected
    ranges or formula-only image results.
    """

    expression: sp.Expr
    formula: sp.Expr
    value_symbol: sp.Symbol
    variables: tuple[sp.Symbol, ...]
    infimum: sp.Expr | None
    supremum: sp.Expr | None
    minimum_attained: bool | None
    maximum_attained: bool | None
    minimizers: tuple[Mapping[sp.Symbol, sp.Expr], ...] = ()
    maximizers: tuple[Mapping[sp.Symbol, sp.Expr], ...] = ()
    method: str = "qe_image"
    diagnostics: Mapping[str, object] = field(default_factory=dict)
    is_interval: bool | None = None
    interval_count: int | None = None

    @property
    def range_condition(self) -> sp.Expr:
        """Alias for ``formula`` emphasizing image semantics."""

        return self.formula

    @property
    def lower_bound(self) -> sp.Expr | None:
        """Alias for the infimum of the range, when known."""

        return self.infimum

    @property
    def upper_bound(self) -> sp.Expr | None:
        """Alias for the supremum of the range, when known."""

        return self.supremum

    @property
    def lower_bound_attained(self) -> bool | None:
        """Whether the lower bound is attained as an actual value."""

        return self.minimum_attained

    @property
    def upper_bound_attained(self) -> bool | None:
        """Whether the upper bound is attained as an actual value."""

        return self.maximum_attained


@dataclass(frozen=True)
class OptimizationResult:
    """Exact optimum summary for supported semialgebraic problems."""

    objective: sp.Expr
    variables: tuple[sp.Symbol, ...]
    value: sp.Expr
    points: tuple[Mapping[sp.Symbol, sp.Expr], ...]
    attained: bool
    kind: str
    method: str = "critical_point_enumeration"
    diagnostics: Mapping[str, object] = field(default_factory=dict)

    @property
    def point(self) -> Mapping[sp.Symbol, sp.Expr] | None:
        return self.points[0] if self.points else None


@dataclass(frozen=True)
class _Candidate:
    value: sp.Expr
    point: Mapping[sp.Symbol, sp.Expr] | None
    attained: bool


def _as_real_symbol(var: sp.Symbol | str) -> sp.Symbol:
    return sp.Symbol(var, real=True) if isinstance(var, str) else var


def _normalize_variables(
    variables: Sequence[sp.Symbol | str] | None, expr: sp.Expr
) -> tuple[sp.Symbol, ...]:
    out: list[sp.Symbol] = []
    seen: set[sp.Symbol] = set()
    if variables is not None:
        for var in variables:
            sym = _as_real_symbol(var)
            if sym not in seen:
                out.append(sym)
                seen.add(sym)
    for sym in sorted(expr.free_symbols, key=lambda item: item.name):
        if sym not in seen:
            out.append(sym)
            seen.add(sym)
    return tuple(out)


def _normalize_formula(formula: FormulaLike | Iterable[FormulaLike] | None) -> sp.Expr:
    if formula is None:
        return sp.true
    if isinstance(formula, (list, tuple, set, frozenset)):
        pieces = [sp.sympify(piece) for piece in formula]
        return sp.And(*pieces) if pieces else sp.true
    if formula is True:
        return sp.true
    if formula is False:
        return sp.false
    if isinstance(formula, (sp.Basic, sp.logic.boolalg.Boolean)):
        return formula  # type: ignore[return-value]
    return to_sympy(formula)  # type: ignore[arg-type]


def _relation_parts(atom: sp.Expr) -> tuple[sp.Expr, str]:
    if isinstance(atom, Equality):
        return sp.expand(atom.lhs - atom.rhs), "=="
    if isinstance(atom, Unequality):
        return sp.expand(atom.lhs - atom.rhs), "!="
    if isinstance(atom, StrictLessThan):
        return sp.expand(atom.lhs - atom.rhs), "<"
    if isinstance(atom, LessThan):
        return sp.expand(atom.lhs - atom.rhs), "<="
    if isinstance(atom, StrictGreaterThan):
        return sp.expand(atom.lhs - atom.rhs), ">"
    if isinstance(atom, GreaterThan):
        return sp.expand(atom.lhs - atom.rhs), ">="
    raise TypeError(f"expected a relational atom, got {atom!r}")


def _atoms(condition: sp.Expr) -> tuple[sp.Expr, ...]:
    if condition is sp.true or isinstance(condition, BooleanTrue):
        return ()
    if condition is sp.false or isinstance(condition, BooleanFalse):
        return (sp.false,)
    if isinstance(condition, SymAnd):
        out: list[sp.Expr] = []
        for arg in condition.args:
            out.extend(_atoms(arg))
        return tuple(out)
    if isinstance(condition, (SymOr, SymNot)):
        raise NotImplementedError(
            "optimization candidate enumeration currently supports conjunctions"
        )
    if getattr(condition, "is_Relational", False):
        return (condition,)
    raise TypeError(f"unsupported formula expression: {condition!r}")


def _is_feasible(condition: sp.Expr, point: Mapping[sp.Symbol, sp.Expr]) -> bool:
    try:
        return bool(satisfies_formula(condition, point, strict=False))
    except _EXPECTED_ERRORS:
        value = sp.simplify(condition.subs(point))
        if value is sp.true or isinstance(value, BooleanTrue):
            return True
        if value is sp.false or isinstance(value, BooleanFalse):
            return False
        return bool(value)


def _finite_compare(a: sp.Expr, b: sp.Expr) -> int:
    a = sp.sympify(a)
    b = sp.sympify(b)
    if a == b or sp.simplify(a - b) == 0:
        return 0
    if a in (sp.oo, -sp.oo) or b in (sp.oo, -sp.oo):
        return -1 if a < b else 1
    diff = sp.simplify(a - b)
    if diff.is_negative:
        return -1
    if diff.is_positive:
        return 1
    numeric = sp.N(diff, 80)
    if numeric < 0:
        return -1
    if numeric > 0:
        return 1
    return 0


def _best_candidates(
    candidates: Iterable[_Candidate], *, kind: str
) -> tuple[sp.Expr, tuple[Mapping[sp.Symbol, sp.Expr], ...], bool]:
    cand = list(candidates)
    if not cand:
        raise ValueError("no feasible candidate points were found")
    best = cand[0]
    for item in cand[1:]:
        cmp = _finite_compare(item.value, best.value)
        if (kind == "min" and cmp < 0) or (kind == "max" and cmp > 0):
            best = item
    tied = [item for item in cand if _finite_compare(item.value, best.value) == 0]
    attained = any(item.attained for item in tied)
    points: list[Mapping[sp.Symbol, sp.Expr]] = []
    seen: set[tuple[tuple[str, str], ...]] = set()
    for item in tied:
        if item.attained and item.point is not None:
            key = tuple(
                sorted((sp.sstr(k), sp.sstr(sp.simplify(v))) for k, v in item.point.items())
            )
            if key not in seen:
                points.append(dict(item.point))
                seen.add(key)
    return sp.simplify(best.value), tuple(points), attained


def _limits_for_interval(
    objective: sp.Expr, variable: sp.Symbol, lo: sp.Expr, hi: sp.Expr
) -> tuple[_Candidate, ...]:
    out: list[_Candidate] = []
    if lo != -sp.oo:
        try:
            out.append(
                _Candidate(sp.simplify(sp.limit(objective, variable, lo, dir="+")), None, False)
            )
        except _EXPECTED_ERRORS:
            pass
    else:
        try:
            out.append(_Candidate(sp.simplify(sp.limit(objective, variable, -sp.oo)), None, False))
        except _EXPECTED_ERRORS:
            pass
    if hi != sp.oo:
        try:
            out.append(
                _Candidate(sp.simplify(sp.limit(objective, variable, hi, dir="-")), None, False)
            )
        except _EXPECTED_ERRORS:
            pass
    else:
        try:
            out.append(_Candidate(sp.simplify(sp.limit(objective, variable, sp.oo)), None, False))
        except _EXPECTED_ERRORS:
            pass
    return tuple(out)


def _univariate_candidates(
    objective: sp.Expr, condition: sp.Expr, variable: sp.Symbol
) -> tuple[_Candidate, ...]:
    candidates: list[_Candidate] = []
    intervals = _one_dimensional_intervals(condition, variable, None)
    cuts: set[str] = set()
    cut_values: list[sp.Expr] = []
    for poly in _relational_polynomials(condition):
        if poly.free_symbols <= {variable}:
            for root in _finite_real_roots(poly, variable):
                key = sp.sstr(root)
                if key not in cuts:
                    cuts.add(key)
                    cut_values.append(root)
    # Zero-dimensional feasible pieces and closed interval endpoints.
    for point_value in cut_values:
        point = {variable: sp.simplify(point_value)}
        if _is_feasible(condition, point):
            candidates.append(_Candidate(sp.simplify(objective.subs(point)), point, True))
    derivative = sp.diff(objective, variable)
    derivative_roots: tuple[sp.Expr, ...]
    if derivative == 0:
        derivative_roots = ()
    else:
        try:
            derivative_roots = _finite_real_roots(derivative, variable)
        except _EXPECTED_ERRORS:
            derivative_roots = ()
    for lo, hi in intervals:
        candidates.extend(_limits_for_interval(objective, variable, lo, hi))
        for root in derivative_roots:
            if (lo == -sp.oo or _finite_compare(lo, root) < 0) and (
                hi == sp.oo or _finite_compare(root, hi) < 0
            ):
                point = {variable: sp.simplify(root)}
                if _is_feasible(condition, point):
                    candidates.append(_Candidate(sp.simplify(objective.subs(point)), point, True))
        # Include a simple interior point for constant objectives.
        if derivative == 0:
            if lo == -sp.oo and hi == sp.oo:
                sample = sp.Integer(0)
            elif lo == -sp.oo:
                sample = hi - 1
            elif hi == sp.oo:
                sample = lo + 1
            else:
                sample = sp.simplify((lo + hi) / 2)
            point = {variable: sample}
            if _is_feasible(condition, point):
                candidates.append(_Candidate(sp.simplify(objective.subs(point)), point, True))
    if condition is sp.true or isinstance(condition, BooleanTrue):
        # No relational polynomials means the whole real line.
        if not intervals:
            candidates.extend(_limits_for_interval(objective, variable, -sp.oo, sp.oo))
            for root in derivative_roots:
                point = {variable: sp.simplify(root)}
                candidates.append(_Candidate(sp.simplify(objective.subs(point)), point, True))
    return tuple(candidates)


def _active_boundary_polys(
    condition: sp.Expr, variables: Sequence[sp.Symbol]
) -> tuple[sp.Expr, ...]:
    polys: list[sp.Expr] = []
    seen: set[str] = set()
    for atom in _atoms(condition):
        if atom is sp.false:
            return ()
        expr, op = _relation_parts(atom)
        if op == "!=":
            continue
        if expr.free_symbols <= set(variables):
            key = sp.sstr(sp.factor(expr))
            if key not in seen:
                polys.append(expr)
                seen.add(key)
    return tuple(polys)


def _solutions_to_points(
    solutions: object, variables: Sequence[sp.Symbol]
) -> tuple[dict[sp.Symbol, sp.Expr], ...]:
    points: list[dict[sp.Symbol, sp.Expr]] = []
    if isinstance(solutions, dict):
        solutions = [solutions]
    if isinstance(solutions, (list, tuple, set)):
        for sol in solutions:
            if isinstance(sol, dict):
                if all(var in sol for var in variables):
                    points.append({var: sp.simplify(sol[var]) for var in variables})
            elif isinstance(sol, (list, tuple)) and len(sol) >= len(variables):
                points.append({var: sp.simplify(sol[i]) for i, var in enumerate(variables)})
    return tuple(points)


def _real_point(point: Mapping[sp.Symbol, sp.Expr]) -> bool:
    for value in point.values():
        value = sp.simplify(value)
        # The first public optimizer only accepts isolated candidate points.
        # Parametric solution families require recursive boundary optimization.
        if value.free_symbols:
            return False
        if value.has(sp.I):
            if sp.simplify(sp.im(value)) != 0:
                return False
        if value.is_real is False:
            return False
    return True


def _multivariate_candidates(
    objective: sp.Expr, condition: sp.Expr, variables: tuple[sp.Symbol, ...]
) -> tuple[_Candidate, ...]:
    if len(variables) != 2:
        raise NotImplementedError("multivariate optimization currently supports two variables")
    x, y = variables
    candidates: list[_Candidate] = []
    grad = [sp.diff(objective, var) for var in variables]
    try:
        candidates.extend(
            _Candidate(sp.simplify(objective.subs(point)), point, True)
            for point in _solutions_to_points(sp.solve(grad, variables, dict=True), variables)
            if _real_point(point) and _is_feasible(condition, point)
        )
    except _EXPECTED_ERRORS:
        pass
    boundaries = _active_boundary_polys(condition, variables)
    # Critical points on one active boundary via Lagrange multipliers.
    for boundary in boundaries:
        lam = sp.Symbol("lambda_semialg", real=True)
        equations = [
            boundary,
            *(sp.diff(objective, var) - lam * sp.diff(boundary, var) for var in variables),
        ]
        try:
            for point in _solutions_to_points(
                sp.solve(equations, (*variables, lam), dict=True), variables
            ):
                if _real_point(point) and _is_feasible(condition, point):
                    candidates.append(_Candidate(sp.simplify(objective.subs(point)), point, True))
        except _EXPECTED_ERRORS:
            pass
    # Boundary intersections/vertices.
    for b1, b2 in combinations(boundaries, 2):
        try:
            for point in _solutions_to_points(sp.solve([b1, b2], variables, dict=True), variables):
                if _real_point(point) and _is_feasible(condition, point):
                    candidates.append(_Candidate(sp.simplify(objective.subs(point)), point, True))
        except _EXPECTED_ERRORS:
            pass
    return tuple(candidates)


def _optimization_candidates(
    objective: sp.Expr, constraints: sp.Expr, variables: tuple[sp.Symbol, ...]
) -> tuple[_Candidate, ...]:
    if constraints is sp.false or isinstance(constraints, BooleanFalse):
        return ()
    if len(variables) == 1:
        return _univariate_candidates(objective, constraints, variables[0])
    if len(variables) == 2:
        return _multivariate_candidates(objective, constraints, variables)
    raise NotImplementedError("optimization currently supports one or two variables")


def _optimize(
    objective: sp.Expr,
    constraints: FormulaLike | Iterable[FormulaLike] | None,
    variables: Sequence[sp.Symbol | str] | None,
    *,
    kind: str,
    domain: str = "reals",
    return_result: bool = True,
) -> OptimizationResult | sp.Expr:
    if domain.lower() not in {"real", "reals", "r", "rr"}:
        raise NotImplementedError(
            "semialgebraic optimization currently supports only the real domain"
        )
    obj = sp.sympify(objective)
    condition = _normalize_formula(constraints)
    vars_ = _normalize_variables(variables, sp.And(condition, sp.Eq(sp.Symbol("_dummy"), obj)))
    vars_ = tuple(var for var in vars_ if var.name != "_dummy")
    if not vars_:
        value = sp.simplify(obj)
        result = OptimizationResult(obj, vars_, value, (), True, kind, "constant")
        return result if return_result else result.value
    candidates = _optimization_candidates(obj, condition, vars_)
    value, points, attained = _best_candidates(candidates, kind=kind)
    result = OptimizationResult(
        obj,
        vars_,
        value,
        points,
        attained,
        kind,
        "critical_point_enumeration",
        {"candidate_count": len(candidates), "constraints": sp.sstr(condition)},
    )
    return result if return_result else result.value


def semialgebraic_minimize(
    objective: sp.Expr,
    constraints: FormulaLike | Iterable[FormulaLike] | None = None,
    variables: Sequence[sp.Symbol | str] | None = None,
    *,
    domain: str = "reals",
    return_result: bool = True,
) -> OptimizationResult | sp.Expr:
    """Return the exact minimum/infimum for supported semialgebraic problems."""

    return _optimize(
        objective, constraints, variables, kind="min", domain=domain, return_result=return_result
    )


def semialgebraic_maximize(
    objective: sp.Expr,
    constraints: FormulaLike | Iterable[FormulaLike] | None = None,
    variables: Sequence[sp.Symbol | str] | None = None,
    *,
    domain: str = "reals",
    return_result: bool = True,
) -> OptimizationResult | sp.Expr:
    """Return the exact maximum/supremum for supported semialgebraic problems."""

    return _optimize(
        objective, constraints, variables, kind="max", domain=domain, return_result=return_result
    )


def _range_formula(
    symbol: sp.Symbol, lo: sp.Expr, hi: sp.Expr, lo_attained: bool, hi_attained: bool
) -> sp.Expr:
    pieces: list[sp.Expr] = []
    if lo != -sp.oo:
        pieces.append(symbol >= lo if lo_attained else symbol > lo)
    if hi != sp.oo:
        pieces.append(symbol <= hi if hi_attained else symbol < hi)
    if not pieces:
        return sp.true
    return sp.And(*pieces)


def _relation_for_function_graph(expr: sp.Expr, value_symbol: sp.Symbol) -> tuple[sp.Expr, sp.Expr]:
    """Return ``(relation, domain)`` for ``value_symbol == expr``.

    Polynomial expressions produce ``Eq(value_symbol - expr, 0)``. Rational
    expressions are handled by clearing denominators and adding the denominator
    nonzero condition, keeping the image problem semialgebraic.
    """

    numerator, denominator = sp.fraction(sp.together(expr))
    numerator = sp.expand(numerator)
    denominator = sp.expand(denominator)
    relation = sp.Eq(sp.expand(value_symbol * denominator - numerator), 0)
    domain = sp.true if denominator == 1 else sp.Ne(denominator, 0)
    return relation, domain


def _is_supported_sqrt(expr: sp.Expr) -> bool:
    """Return True for the square-root form supported by graph conversion."""

    return isinstance(expr, sp.Pow) and expr.exp == sp.Rational(1, 2)


def _is_semialgebraic_special(expr: sp.Expr) -> bool:
    """Return True for expression heads handled by semialgebraic graph conversion."""

    return (
        expr.func is sp.Abs
        or expr.func is sp.Max
        or expr.func is sp.Min
        or isinstance(expr, sp.Piecewise)
        or _is_supported_sqrt(expr)
    )


def _has_semialgebraic_special(expr: sp.Expr) -> bool:
    """Whether ``expr`` contains supported non-rational semialgebraic heads."""

    return any(_is_semialgebraic_special(item) for item in sp.preorder_traversal(expr))


def _fresh_graph_symbol(counter: list[int]) -> sp.Symbol:
    symbol = sp.Symbol(f"_semialg_graph_aux_{counter[0]}", real=True)
    counter[0] += 1
    return symbol


def _graph_formula_for_expression(
    expr: sp.Expr,
    target: sp.Symbol,
    counter: list[int],
) -> tuple[sp.Expr, tuple[sp.Symbol, ...]]:
    """Return a semialgebraic graph formula for ``target == expr``.

    The graph conversion supports common expression heads whose graphs are semialgebraic:
    ``Abs``, square roots, ``Min``, ``Max``, and simple ``Piecewise``. General
    arithmetic combinations are handled by replacing special subexpressions by
    auxiliary graph variables and then clearing denominators for the remaining
    rational expression.
    """

    expr = sp.sympify(expr)
    aux_symbols: list[sp.Symbol] = []

    if expr.func is sp.Abs:
        arg_value = _fresh_graph_symbol(counter)
        arg_formula, arg_aux = _graph_formula_for_expression(expr.args[0], arg_value, counter)
        aux_symbols.extend((arg_value, *arg_aux))
        formula = sp.And(
            arg_formula,
            target >= 0,
            sp.Or(sp.Eq(target, arg_value), sp.Eq(target, -arg_value)),
        )
        return formula, tuple(aux_symbols)

    if _is_supported_sqrt(expr):
        base_value = _fresh_graph_symbol(counter)
        base_formula, base_aux = _graph_formula_for_expression(expr.base, base_value, counter)
        aux_symbols.extend((base_value, *base_aux))
        formula = sp.And(
            base_formula,
            base_value >= 0,
            target >= 0,
            sp.Eq(target**2, base_value),
        )
        return formula, tuple(aux_symbols)

    if expr.func in {sp.Max, sp.Min}:
        arg_values: list[sp.Symbol] = []
        arg_formulas: list[sp.Expr] = []
        for arg in expr.args:
            arg_value = _fresh_graph_symbol(counter)
            arg_formula, arg_aux = _graph_formula_for_expression(arg, arg_value, counter)
            arg_values.append(arg_value)
            arg_formulas.append(arg_formula)
            aux_symbols.extend((arg_value, *arg_aux))
        branch_formulas: list[sp.Expr] = []
        for _index, arg_value in enumerate(arg_values):
            if expr.func is sp.Max:
                order = sp.And(*(arg_value >= other for other in arg_values))
            else:
                order = sp.And(*(arg_value <= other for other in arg_values))
            branch_formulas.append(sp.And(sp.Eq(target, arg_value), order))
        return sp.And(*arg_formulas, sp.Or(*branch_formulas)), tuple(aux_symbols)

    if isinstance(expr, sp.Piecewise):
        branches: list[sp.Expr] = []
        previous_conditions: list[sp.Expr] = []
        for branch_expr, branch_condition in expr.args:
            effective_condition = sp.And(
                branch_condition, *(sp.Not(cond) for cond in previous_conditions)
            )
            branch_formula, branch_aux = _graph_formula_for_expression(branch_expr, target, counter)
            aux_symbols.extend(branch_aux)
            branches.append(sp.And(effective_condition, branch_formula))
            previous_conditions.append(branch_condition)
        return sp.Or(*branches), tuple(aux_symbols)

    if not expr.args:
        relation, domain = _relation_for_function_graph(expr, target)
        return sp.And(domain, relation), ()

    constraints: list[sp.Expr] = []
    replacements: dict[sp.Expr, sp.Symbol] = {}
    all_aux: list[sp.Symbol] = []

    def replace_specials(item: sp.Expr) -> sp.Expr:
        item = sp.sympify(item)
        if _is_semialgebraic_special(item):
            if item in replacements:
                return replacements[item]
            aux = _fresh_graph_symbol(counter)
            formula, nested_aux = _graph_formula_for_expression(item, aux, counter)
            replacements[item] = aux
            constraints.append(formula)
            all_aux.extend((aux, *nested_aux))
            return aux
        if not item.args:
            return item
        new_args = tuple(replace_specials(arg) for arg in item.args)
        if new_args == item.args:
            return item
        return item.func(*new_args)

    transformed = replace_specials(expr)
    relation, domain = _relation_for_function_graph(transformed, target)
    aux_symbols.extend(all_aux)
    return sp.And(*(constraints + [domain, relation])), tuple(dict.fromkeys(aux_symbols))


def _range_formula_from_set(range_set: sp.Set, value_symbol: sp.Symbol) -> sp.Expr | None:
    """Convert a one-dimensional SymPy set into a range condition."""

    try:
        intervals = _intervals_from_range_set(range_set)
    except _EXPECTED_ERRORS:
        intervals = None
    if intervals is not None:
        pieces: list[sp.Expr] = []
        for left, right, left_closed, right_closed in intervals:
            if left == right:
                pieces.append(sp.Eq(value_symbol, left))
                continue
            lower = (
                sp.true
                if left == -sp.oo
                else (value_symbol >= left if left_closed else value_symbol > left)
            )
            upper = (
                sp.true
                if right == sp.oo
                else (value_symbol <= right if right_closed else value_symbol < right)
            )
            pieces.append(sp.And(lower, upper))
        return sp.Or(*pieces) if pieces else sp.false
    try:
        return range_set.as_relational(value_symbol)
    except _EXPECTED_ERRORS:
        return None


def _domain_set_from_condition(condition: sp.Expr, variable: sp.Symbol) -> sp.Set | None:
    """Return the one-variable domain represented by simple interval constraints."""

    if condition is sp.true or isinstance(condition, BooleanTrue):
        return sp.S.Reals
    if condition is sp.false or isinstance(condition, BooleanFalse):
        return sp.S.EmptySet
    atoms = condition.args if isinstance(condition, SymAnd) else (condition,)
    lower = -sp.oo
    upper = sp.oo
    left_open = False
    right_open = False
    for atom in atoms:
        if not getattr(atom, "is_Relational", False):
            return None
        expr, op = _relation_parts(atom)
        try:
            poly = sp.Poly(expr, variable)
        except _EXPECTED_ERRORS:
            return None
        if poly.degree() != 1:
            return None
        coeff = sp.simplify(poly.coeff_monomial(variable))
        const = sp.simplify(poly.eval(0))
        if coeff == 0:
            truth = bool(_relation_from_scalar(const, op))
            if not truth:
                return sp.S.EmptySet
            continue
        bound = sp.simplify(-const / coeff)
        if coeff.is_negative:
            op = {"<": ">", "<=": ">=", ">": "<", ">=": "<="}.get(op, op)
        if op in {">", ">="}:
            if lower == -sp.oo or _finite_compare(lower, bound) < 0:
                lower = bound
                left_open = op == ">"
            elif sp.simplify(lower - bound) == 0:
                left_open = left_open or op == ">"
        elif op in {"<", "<="}:
            if upper == sp.oo or _finite_compare(bound, upper) < 0:
                upper = bound
                right_open = op == "<"
            elif sp.simplify(upper - bound) == 0:
                right_open = right_open or op == "<"
        elif op == "==":
            lower = upper = bound
            left_open = right_open = False
        else:
            return None
    if lower != -sp.oo and upper != sp.oo and _finite_compare(lower, upper) > 0:
        return sp.S.EmptySet
    return sp.Interval(lower, upper, left_open=left_open, right_open=right_open)


def _relation_from_scalar(value: sp.Expr, op: str) -> bool:
    if op == "<":
        return bool(value < 0)
    if op == "<=":
        return bool(value <= 0)
    if op == ">":
        return bool(value > 0)
    if op == ">=":
        return bool(value >= 0)
    if op == "==":
        return bool(sp.simplify(value) == 0)
    if op == "!=":
        return bool(sp.simplify(value) != 0)
    raise ValueError(op)


def _interval_from_bounds(lo: sp.Expr, hi: sp.Expr) -> sp.Interval:
    if _finite_compare(lo, hi) > 0:
        lo, hi = hi, lo
    return sp.Interval(sp.simplify(lo), sp.simplify(hi))


def _poly_range_on_interval(
    expr: sp.Expr, variable: sp.Symbol, domain: sp.Interval
) -> sp.Set | None:
    """Return the range of a univariate polynomial on a closed interval."""

    if not isinstance(domain, sp.Interval):
        return None
    try:
        poly = sp.Poly(sp.expand(expr), variable)
    except _EXPECTED_ERRORS:
        return None
    if domain.start in (-sp.oo, sp.oo) or domain.end in (-sp.oo, sp.oo):
        return None
    candidates = [domain.start, domain.end]
    derivative = sp.diff(poly.as_expr(), variable)
    try:
        for root in sp.solve(sp.Eq(derivative, 0), variable):
            if root.is_real is False:
                continue
            if bool(root >= domain.start) and bool(root <= domain.end):
                candidates.append(root)
    except _EXPECTED_ERRORS:
        pass
    values = [sp.simplify(expr.subs(variable, item)) for item in candidates]
    lo = min(values, key=_bound_sort_key)
    hi = max(values, key=_bound_sort_key)
    return _interval_from_bounds(lo, hi)


def _range_via_abs_affine(
    expr: sp.Expr,
    condition: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    value_symbol: sp.Symbol,
) -> sp.Expr | None:
    """Fast range for affine expressions in ``Abs(variable)``."""

    if len(variables) != 1:
        return None
    variable = variables[0]
    abs_expr = sp.Abs(variable)
    if not expr.has(abs_expr):
        return None
    domain = _domain_set_from_condition(condition, variable)
    if domain is None:
        return None
    u = sp.Symbol("_abs_value", real=True)
    transformed = expr.xreplace({abs_expr: u})
    try:
        poly = sp.Poly(sp.expand(transformed), u)
    except _EXPECTED_ERRORS:
        return None
    if poly.degree() != 1 or any(sym != u for sym in transformed.free_symbols):
        return None
    coeff = sp.simplify(poly.coeff_monomial(u))
    const = sp.simplify(poly.eval(0))
    if coeff == 0:
        return _range_formula_from_set(sp.FiniteSet(const), value_symbol)
    if domain is sp.S.Reals:
        abs_range = sp.Interval(0, sp.oo)
    elif isinstance(domain, sp.Interval):
        values = []
        if domain.start != -sp.oo:
            values.append(sp.Abs(domain.start))
        if domain.end != sp.oo:
            values.append(sp.Abs(domain.end))
        crosses_zero = domain.start <= 0 and domain.end >= 0
        lower = sp.Integer(0) if bool(crosses_zero) else min(values, key=_bound_sort_key)
        upper = (
            sp.oo
            if domain.start == -sp.oo or domain.end == sp.oo
            else max(values, key=_bound_sort_key)
        )
        abs_range = sp.Interval(lower, upper)
    else:
        return None
    left = sp.simplify(coeff * abs_range.start + const) if abs_range.start != -sp.oo else -sp.oo
    right = (
        sp.simplify(coeff * abs_range.end + const)
        if abs_range.end != sp.oo
        else (sp.oo if coeff.is_positive is not False else -sp.oo)
    )
    return _range_formula_from_set(_interval_from_bounds(left, right), value_symbol)


def _range_via_sqrt_quadratic(
    expr: sp.Expr,
    condition: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    value_symbol: sp.Symbol,
) -> sp.Expr | None:
    """Fast range for square roots of nonnegative quadratic radicands."""

    if len(variables) != 1 or not _is_supported_sqrt(expr):
        return None
    variable = variables[0]
    domain = _domain_set_from_condition(condition, variable)
    if domain is None:
        return None
    try:
        base_poly = sp.Poly(sp.expand(expr.base), variable)
    except _EXPECTED_ERRORS:
        return None
    if base_poly.degree() > 2:
        return None
    coeff2 = sp.simplify(base_poly.coeff_monomial(variable**2))
    coeff1 = sp.simplify(base_poly.coeff_monomial(variable))
    coeff0 = sp.simplify(base_poly.coeff_monomial(1))
    if coeff2 == -1 and coeff1 == 0 and coeff0.is_real:
        natural = sp.Interval(-sp.sqrt(coeff0), sp.sqrt(coeff0))
    else:
        return None
    active_domain = natural if domain is sp.S.Reals else domain.intersect(natural)
    if not isinstance(active_domain, sp.Interval):
        return None
    base_range = _poly_range_on_interval(expr.base, variable, active_domain)
    if not isinstance(base_range, sp.Interval):
        return None
    upper = sp.sqrt(sp.simplify(base_range.end))
    return _range_formula_from_set(sp.Interval(0, upper), value_symbol)


def _range_via_sympy_calculus(
    expr: sp.Expr,
    condition: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    value_symbol: sp.Symbol,
) -> sp.Expr | None:
    """Fast image computation for common univariate expression heads."""

    if len(variables) != 1:
        return None
    variable = variables[0]
    domain = _domain_set_from_condition(condition, variable)
    if domain is None:
        return None
    try:
        from sympy.calculus.util import function_range as sympy_range

        range_set = sympy_range(expr, variable, domain)
    except _EXPECTED_ERRORS:
        return None
    return _range_formula_from_set(range_set, value_symbol)


def _bound_sort_key(value: sp.Expr) -> float:
    if value == -sp.oo:
        return -float("inf")
    if value == sp.oo:
        return float("inf")
    return float(sp.N(value, 50))


def _range_via_max_min(
    expr: sp.Expr,
    condition: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    value_symbol: sp.Symbol,
) -> sp.Expr | None:
    """Fast range for two-argument Max/Min with affine arguments."""

    if len(variables) != 1 or expr.func not in {sp.Max, sp.Min} or len(expr.args) != 2:
        return None
    variable = variables[0]
    domain = _domain_set_from_condition(condition, variable)
    if not isinstance(domain, sp.Interval):
        return None
    arg_ranges: list[sp.Set] = []
    for arg in expr.args:
        try:
            poly = sp.Poly(sp.expand(arg), variable)
        except _EXPECTED_ERRORS:
            return None
        if poly.degree() > 1:
            return None
        coeff = sp.simplify(poly.coeff_monomial(variable))
        const = sp.simplify(poly.eval(0))
        if coeff == 0:
            arg_ranges.append(sp.FiniteSet(const))
            continue
        left = -sp.oo if domain.start == -sp.oo else sp.simplify(coeff * domain.start + const)
        right = sp.oo if domain.end == sp.oo else sp.simplify(coeff * domain.end + const)
        lo, hi = (left, right) if _finite_compare(left, right) <= 0 else (right, left)
        arg_ranges.append(sp.Interval(lo, hi))
    if expr.func is sp.Max:
        lower = max((item.inf for item in arg_ranges), key=_bound_sort_key)
        upper = max((item.sup for item in arg_ranges), key=_bound_sort_key)
    else:
        lower = min((item.inf for item in arg_ranges), key=_bound_sort_key)
        upper = min((item.sup for item in arg_ranges), key=_bound_sort_key)
    return _range_formula_from_set(sp.Interval(lower, upper), value_symbol)


def _range_via_piecewise(
    expr: sp.Expr,
    condition: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    value_symbol: sp.Symbol,
) -> sp.Expr | None:
    """Fast range for univariate Piecewise expressions by branch images."""

    if len(variables) != 1 or not isinstance(expr, sp.Piecewise):
        return None
    variable = variables[0]
    pieces: list[sp.Set] = []
    previous: list[sp.Expr] = []
    for branch_expr, branch_cond in expr.args:
        effective = sp.And(condition, branch_cond, *(sp.Not(item) for item in previous))
        previous.append(branch_cond)
        domain = _domain_set_from_condition(effective, variable)
        if domain is None or domain is sp.S.EmptySet:
            continue
        if isinstance(domain, sp.Union):
            domains = domain.args
        else:
            domains = (domain,)
        for branch_domain in domains:
            image = _poly_range_on_interval(branch_expr, variable, branch_domain)
            if image is None:
                return None
            pieces.append(image)
    if not pieces:
        return sp.false
    return _range_formula_from_set(sp.Union(*pieces), value_symbol)


def _try_direct_special_image(
    expr: sp.Expr,
    condition: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    value_symbol: sp.Symbol,
) -> sp.Expr | None:
    """Use fast univariate range formulas before constructing a CAD image problem."""

    for helper in (
        _range_via_abs_affine,
        _range_via_sqrt_quadratic,
        _range_via_max_min,
        _range_via_piecewise,
        _range_via_sympy_calculus,
    ):
        result = helper(expr, condition, variables, value_symbol)
        if result is not None:
            return result
    return None


def _try_semialgebraic_graph_image(
    expr: sp.Expr,
    condition: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    value_symbol: sp.Symbol,
) -> sp.Expr | None:
    """Use CAD/QE on a semialgebraic graph relation for semialgebraic graph ranges."""

    if not _has_semialgebraic_special(expr):
        return None
    try:
        graph_formula, aux_symbols = _graph_formula_for_expression(expr, value_symbol, [0])
        image_formula = sp.And(condition, graph_formula)
        elimination_variables = tuple(dict.fromkeys((*variables, *aux_symbols)))
        parsed = parse_formula(image_formula)
        result = qe_by_complete_cad(
            (value_symbol, *elimination_variables),
            tuple(("exists", variable) for variable in elimination_variables),
            parsed,
            free_variables=(value_symbol,),
        )
    except _EXPECTED_ERRORS:
        return None
    return sp.simplify(result.formula)


def _substitute_formula(expr: sp.Expr, substitution: Mapping[sp.Symbol, sp.Expr]) -> sp.Expr:
    if expr is sp.true or isinstance(expr, BooleanTrue):
        return sp.true
    if expr is sp.false or isinstance(expr, BooleanFalse):
        return sp.false
    return sp.simplify(expr.subs(substitution))


def _try_affine_univariate_image(
    expr: sp.Expr,
    condition: sp.Expr,
    variable: sp.Symbol,
    value_symbol: sp.Symbol,
) -> sp.Expr | None:
    """Fast exact image for an affine one-variable map.

    This covers important disconnected ranges without needing a full CAD run,
    e.g. the image of ``x`` over ``x <= -1 or x >= 1``.
    """

    try:
        poly = sp.Poly(sp.expand(expr), variable)
    except _EXPECTED_ERRORS:
        return None
    if poly.degree() != 1:
        return None
    coefficient = poly.coeff_monomial(variable)
    constant = poly.eval(0)
    if sp.simplify(coefficient) == 0:
        sample = {variable: sp.Integer(0)}
        return (
            sp.Eq(value_symbol, sp.simplify(constant))
            if _is_feasible(condition, sample)
            else sp.false
        )
    inverse = sp.simplify((value_symbol - constant) / coefficient)
    return _substitute_formula(condition, {variable: inverse})


def _try_solved_graph_image(
    expr: sp.Expr,
    condition: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    value_symbol: sp.Symbol,
) -> sp.Expr | None:
    """Try exact image computation by solving the graph equation.

    This conservative helper supports univariate polynomial/rational maps when
    SymPy can solve the graph equation explicitly enough to substitute branches
    back into the domain condition.
    """

    if len(variables) != 1:
        return None
    variable = variables[0]
    affine = _try_affine_univariate_image(expr, condition, variable, value_symbol)
    if affine is not None:
        return sp.simplify(affine)
    numerator, denominator = sp.fraction(sp.together(expr))
    equation = sp.expand(value_symbol * denominator - numerator)
    try:
        solutions = sp.solve(equation, variable)
    except _EXPECTED_ERRORS:
        return None
    if not solutions:
        return None
    branches: list[sp.Expr] = []
    for solution in solutions:
        if variable in solution.free_symbols or solution.has(sp.I):
            return None
        if any(
            isinstance(pow_expr, sp.Pow) and pow_expr.exp.is_integer is False
            for pow_expr in solution.atoms(sp.Pow)
        ):
            return None
        branch_condition = _substitute_formula(condition, {variable: solution})
        branch_domain = _substitute_formula(sp.Ne(denominator, 0), {variable: solution})
        branches.append(sp.And(branch_condition, branch_domain))
    return sp.simplify(sp.Or(*branches))


def _try_complete_cad_image(
    expr: sp.Expr,
    condition: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    value_symbol: sp.Symbol,
) -> sp.Expr | None:
    """Run the CAD/QE image formulation for polynomial/rational expressions."""

    relation, domain_constraint = _relation_for_function_graph(expr, value_symbol)
    image_formula = sp.And(condition, domain_constraint, relation)
    try:
        parsed = parse_formula(image_formula)
        result = qe_by_complete_cad(
            (value_symbol, *variables),
            tuple(("exists", variable) for variable in variables),
            parsed,
            free_variables=(value_symbol,),
        )
    except _EXPECTED_ERRORS:
        return None
    return sp.simplify(result.formula)


def _simplify_range_condition(formula: sp.Expr, value_symbol: sp.Symbol) -> sp.Expr:
    if formula is sp.true or formula is sp.false:
        return formula
    if isinstance(formula, SymAnd):
        return sp.And(*(_simplify_range_condition(arg, value_symbol) for arg in formula.args))
    if isinstance(formula, SymOr):
        return sp.Or(*(_simplify_range_condition(arg, value_symbol) for arg in formula.args))
    if getattr(formula, "is_Relational", False) and formula.free_symbols <= {value_symbol}:
        return formula
    return formula


def _intervals_from_range_set(
    range_set: sp.Set,
) -> tuple[tuple[sp.Expr, sp.Expr, bool, bool], ...] | None:
    """Return intervals as ``(left, right, left_closed, right_closed)``.

    The helper intentionally accepts only one-dimensional set forms that have
    direct range-bound semantics: intervals, finite point sets, and finite
    unions of those. More exotic sets leave metadata unknown without changing
    the primary range formula.
    """

    if range_set is sp.S.EmptySet:
        return ()
    if range_set is sp.S.Reals or range_set is sp.S.UniversalSet:
        return ((-sp.oo, sp.oo, False, False),)
    if isinstance(range_set, sp.Interval):
        return (
            (range_set.start, range_set.end, not range_set.left_open, not range_set.right_open),
        )
    if isinstance(range_set, sp.FiniteSet):
        return tuple(
            (item, item, True, True) for item in sorted(range_set, key=sp.default_sort_key)
        )
    if isinstance(range_set, sp.Union):
        pieces: list[tuple[sp.Expr, sp.Expr, bool, bool]] = []
        for arg in range_set.args:
            intervals = _intervals_from_range_set(arg)
            if intervals is None:
                return None
            pieces.extend(intervals)
        return tuple(sorted(pieces, key=lambda item: sp.default_sort_key(item[0])))
    return None


def _compare_bounds_for_sort(a: sp.Expr, b: sp.Expr) -> int:
    if a == b:
        return 0
    if a == -sp.oo or b == sp.oo:
        return -1
    if a == sp.oo or b == -sp.oo:
        return 1
    return _finite_compare(a, b)


def _range_metadata_from_formula(
    formula: sp.Expr,
    value_symbol: sp.Symbol,
) -> tuple[sp.Expr | None, sp.Expr | None, bool | None, bool | None, bool | None, int | None]:
    """Extract interval metadata from a one-variable range condition.

    The range formula remains authoritative. These fields are summaries for
    interval-like answers and are intentionally ``None`` when extraction is not
    reliable.
    """

    range_set = _domain_set_from_condition(formula, value_symbol)
    if range_set is None:
        try:
            range_set = formula.as_set()
        except _EXPECTED_ERRORS:
            return None, None, None, None, None, None
    intervals = _intervals_from_range_set(range_set)
    if intervals is None:
        return None, None, None, None, None, None
    interval_count = len(intervals)
    if interval_count == 0:
        return None, None, False, False, True, 0
    left = intervals[0][0]
    right = intervals[0][1]
    left_closed = intervals[0][2]
    right_closed = intervals[0][3]
    for item in intervals[1:]:
        if _compare_bounds_for_sort(item[0], left) < 0:
            left, left_closed = item[0], item[2]
        if _compare_bounds_for_sort(right, item[1]) < 0:
            right, right_closed = item[1], item[3]
    lower_attained = bool(left_closed) if left not in (-sp.oo, sp.oo) else False
    upper_attained = bool(right_closed) if right not in (-sp.oo, sp.oo) else False
    is_interval = interval_count == 1
    return left, right, lower_attained, upper_attained, is_interval, interval_count


def _range_result_from_formula(
    expression: sp.Expr,
    formula: sp.Expr,
    value_symbol: sp.Symbol,
    variables: tuple[sp.Symbol, ...],
    method: str,
    diagnostics: Mapping[str, object],
) -> FunctionRangeResult:
    simplified_formula = _simplify_range_condition(formula, value_symbol)
    lo, hi, lo_attained, hi_attained, is_interval, interval_count = _range_metadata_from_formula(
        simplified_formula,
        value_symbol,
    )
    return FunctionRangeResult(
        expression,
        simplified_formula,
        value_symbol,
        variables,
        lo,
        hi,
        lo_attained,
        hi_attained,
        (),
        (),
        method,
        {
            **dict(diagnostics),
            "metadata_source": "range_formula_as_set" if interval_count is not None else "unknown",
        },
        is_interval,
        interval_count,
    )


def _range_via_optimization_bounds(
    expr: sp.Expr,
    condition: sp.Expr,
    vars_: tuple[sp.Symbol, ...],
    val_sym: sp.Symbol,
    domain: str,
) -> FunctionRangeResult:
    minimum = semialgebraic_minimize(expr, condition, vars_, domain=domain, return_result=True)
    maximum = semialgebraic_maximize(expr, condition, vars_, domain=domain, return_result=True)
    assert isinstance(minimum, OptimizationResult)
    assert isinstance(maximum, OptimizationResult)
    formula = sp.simplify(
        _range_formula(val_sym, minimum.value, maximum.value, minimum.attained, maximum.attained)
    )
    _, _, _, _, is_interval, interval_count = _range_metadata_from_formula(formula, val_sym)
    return FunctionRangeResult(
        expr,
        formula,
        val_sym,
        vars_,
        minimum.value,
        maximum.value,
        minimum.attained,
        maximum.attained,
        minimum.points,
        maximum.points,
        "optimization_bounds",
        {"constraints": sp.sstr(condition), "metadata_source": "optimization_bounds"},
        is_interval,
        interval_count,
    )


def function_range(
    expression: sp.Expr,
    constraints: FormulaLike | Iterable[FormulaLike] | None = None,
    variables: Sequence[sp.Symbol | str] | None = None,
    *,
    value_symbol: sp.Symbol | str | None = None,
    domain: str = "reals",
    method: str = "qe",
    return_result: bool = False,
) -> sp.Expr | FunctionRangeResult:
    """Return a quantifier-free formula describing a real function range.

    The preferred direct backend uses the semialgebraic image formulation
    ``exists variables. constraints and value_symbol == expression``. It first
    applies guarded exact graph-elimination shortcuts for common univariate
    images, then tries complete CAD/QE, and finally falls back to the earlier
    optimization-bound summary for ``method='auto'`` or ``method='bounds'``.
    """

    if domain.lower() not in {"real", "reals", "r", "rr"}:
        raise NotImplementedError("function_range currently supports only the real domain")
    requested_method = method.lower()
    if requested_method not in {"qe", "cad", "auto", "bounds", "optimization"}:
        raise ValueError("method must be 'qe', 'cad', 'auto', 'bounds', or 'optimization'")
    expr = sp.sympify(expression)
    condition = _normalize_formula(constraints)
    vars_ = _normalize_variables(variables, sp.And(condition, sp.Eq(sp.Symbol("_dummy"), expr)))
    vars_ = tuple(var for var in vars_ if var.name != "_dummy")
    val_sym = _as_real_symbol(value_symbol or "t")

    if not vars_:
        formula = sp.Eq(val_sym, expr)
        result = FunctionRangeResult(
            expr,
            formula,
            val_sym,
            vars_,
            expr,
            expr,
            True,
            True,
            (),
            (),
            "constant",
            {"metadata_source": "constant"},
            True,
            1,
        )
        return result if return_result else result.formula

    if requested_method in {"bounds", "optimization"}:
        result = _range_via_optimization_bounds(expr, condition, vars_, val_sym, domain)
        return result if return_result else result.formula

    direct_special = _try_direct_special_image(expr, condition, vars_, val_sym)
    if direct_special is not None:
        direct_method = "direct_univariate_image"
        try:
            if len(vars_) == 1 and sp.Poly(sp.expand(expr), vars_[0]).degree() >= 2:
                direct_method = "optimization_bounds_direct_univariate_image"
        except _EXPECTED_ERRORS:
            pass
        result = _range_result_from_formula(
            expr,
            direct_special,
            val_sym,
            vars_,
            direct_method,
            {"constraints": sp.sstr(condition), "requested_method": requested_method},
        )
        return result if return_result else result.formula

    semialgebraic_graph = _try_semialgebraic_graph_image(expr, condition, vars_, val_sym)
    if semialgebraic_graph is not None:
        result = _range_result_from_formula(
            expr,
            semialgebraic_graph,
            val_sym,
            vars_,
            "qe_image_semialgebraic_graph",
            {"constraints": sp.sstr(condition), "requested_method": requested_method},
        )
        return result if return_result else result.formula

    direct = _try_solved_graph_image(expr, condition, vars_, val_sym)
    if direct is not None:
        result = _range_result_from_formula(
            expr,
            direct,
            val_sym,
            vars_,
            "qe_image_solved_graph",
            {"constraints": sp.sstr(condition), "requested_method": requested_method},
        )
        return result if return_result else result.formula

    cad_formula = (
        _try_complete_cad_image(expr, condition, vars_, val_sym)
        if requested_method == "cad"
        else None
    )
    if cad_formula is not None:
        result = _range_result_from_formula(
            expr,
            cad_formula,
            val_sym,
            vars_,
            "qe_image_complete_cad",
            {"constraints": sp.sstr(condition), "requested_method": requested_method},
        )
        return result if return_result else result.formula

    if requested_method == "cad":
        raise NotImplementedError(
            "CAD/QE-image range computation failed for this expression/domain"
        )

    result = _range_via_optimization_bounds(expr, condition, vars_, val_sym, domain)
    result = FunctionRangeResult(
        result.expression,
        result.formula,
        result.value_symbol,
        result.variables,
        result.infimum,
        result.supremum,
        result.minimum_attained,
        result.maximum_attained,
        result.minimizers,
        result.maximizers,
        "optimization_bounds_after_qe_image_fallback",
        {**dict(result.diagnostics), "requested_method": requested_method},
        result.is_interval,
        result.interval_count,
    )
    return result if return_result else result.formula


__all__ = [
    "FunctionRangeResult",
    "OptimizationResult",
    "function_range",
    "semialgebraic_maximize",
    "semialgebraic_minimize",
]
