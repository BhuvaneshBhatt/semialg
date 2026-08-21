from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from functools import cmp_to_key

import sympy as sp
from sympy.logic.boolalg import And as SymAnd
from sympy.logic.boolalg import BooleanFalse, BooleanTrue
from sympy.logic.boolalg import Not as SymNot
from sympy.logic.boolalg import Or as SymOr

from .exact_arithmetic import compare_exact_reals, exact_truth
from .interval_decomposition import finite_real_roots as _finite_real_roots
from .normalization import normalize_formula
from .normalization import normalize_problem_variables as _shared_variables
from .relations import make_zero_relation
from .relations import split_relation as _relational_difference


@dataclass(frozen=True)
class ImplicitFormulaPiece:
    """One disjunctive piece of a semialgebraic formula.

    ``inequalities`` are normalized as expressions required to be ``<= 0``.
    Strict inequalities are relaxed to weak inequalities, matching the common
    region-closure/implicit-description convention used by CAD region code.
    ``equalities`` are normalized as expressions required to be ``== 0``.
    Unequalities are intentionally ignored because they do not change ordinary
    closure or ambient Lebesgue-measure computations.
    """

    inequalities: tuple[sp.Expr, ...]
    equalities: tuple[sp.Expr, ...]

    def as_formula(self) -> sp.Expr:
        parts: list[sp.Expr] = [ineq <= 0 for ineq in self.inequalities]
        parts.extend(sp.Eq(eq, 0) for eq in self.equalities)
        return sp.And(*parts) if parts else sp.true


@dataclass(frozen=True)
class SymbolicBoxBounds:
    """Explicit axis-aligned bounds for a conjunction of simple inequalities."""

    limits: tuple[tuple[sp.Symbol, sp.Expr, sp.Expr], ...]

    def as_dict(self) -> dict[sp.Symbol, tuple[sp.Expr, sp.Expr]]:
        return {var: (lo, hi) for var, lo, hi in self.limits}


@dataclass(frozen=True)
class VerticalBoundCell2D:
    """One explicit vertical cell description in two variables.

    ``x_interval`` is ``(lo, hi)``. If ``lo == hi`` the cell is a vertical
    fiber over a point and contributes zero to ambient two-dimensional measure.
    ``y_bounds`` is a tuple of ``(lower(x), upper(x))`` intervals. The current
    cell object records the exact source formula for semantic checks and uses
    weak inequalities in generated structural formulas; callers that need exact
    open/closed boundary information should use ``source_formula``.
    """

    x_variable: sp.Symbol
    y_variable: sp.Symbol
    x_interval: tuple[sp.Expr, sp.Expr]
    y_bounds: tuple[tuple[sp.Expr, sp.Expr], ...]
    x_condition: sp.Expr = sp.true
    source_formula: sp.Expr | None = None

    @property
    def is_full_dimensional(self) -> bool:
        lo, hi = self.x_interval
        return sp.simplify(hi - lo) != 0 and bool(lo != hi)

    @property
    def dimension(self) -> int:
        if not self.is_full_dimensional:
            return 1 if self.y_bounds else 0
        if all(sp.simplify(upper - lower) == 0 for lower, upper in self.y_bounds):
            return 1
        return 2

    @property
    def bounded(self) -> bool:
        lo, hi = self.x_interval
        if lo == -sp.oo or hi == sp.oo:
            return False
        return all(lower != -sp.oo and upper != sp.oo for lower, upper in self.y_bounds)

    def as_formula(self, *, use_source: bool = True) -> sp.Expr:
        """Return a formula describing the cell.

        When ``use_source`` is true and the original disjunct is available, this
        returns that source formula so strict boundary semantics are preserved.
        Otherwise it returns a weak structural formula from the stored bounds.
        """

        if use_source and self.source_formula is not None:
            return self.source_formula
        x, y = self.x_variable, self.y_variable
        xlo, xhi = self.x_interval
        parts: list[sp.Expr] = []
        if xlo == xhi:
            parts.append(sp.Eq(x, xlo))
        else:
            if xlo != -sp.oo:
                parts.append(x >= xlo)
            if xhi != sp.oo:
                parts.append(x <= xhi)
        y_parts: list[sp.Expr] = []
        for lower, upper in self.y_bounds:
            if sp.simplify(upper - lower) == 0:
                y_parts.append(sp.Eq(y, lower))
            else:
                yp: list[sp.Expr] = []
                if lower != -sp.oo:
                    yp.append(y >= lower)
                if upper != sp.oo:
                    yp.append(y <= upper)
                y_parts.append(sp.And(*yp) if yp else sp.true)
        if y_parts:
            parts.append(sp.Or(*y_parts) if len(y_parts) > 1 else y_parts[0])
        return sp.And(*parts) if parts else sp.true

    def sample_point(self) -> dict[sp.Symbol, sp.Expr]:
        """Return a deterministic interior-style sample for the stored bounds."""

        xlo, xhi = self.x_interval
        if xlo == xhi:
            x_value = xlo
        elif xlo == -sp.oo and xhi == sp.oo:
            x_value = sp.Integer(0)
        elif xlo == -sp.oo:
            x_value = sp.simplify(xhi - 1)
        elif xhi == sp.oo:
            x_value = sp.simplify(xlo + 1)
        else:
            x_value = sp.simplify((xlo + xhi) / 2)
        if not self.y_bounds:
            y_value = sp.Integer(0)
        else:
            lower, upper = self.y_bounds[0]
            lower_value = sp.simplify(lower.subs(self.x_variable, x_value))
            upper_value = sp.simplify(upper.subs(self.x_variable, x_value))
            if lower_value == upper_value:
                y_value = lower_value
            elif lower_value == -sp.oo and upper_value == sp.oo:
                y_value = sp.Integer(0)
            elif lower_value == -sp.oo:
                y_value = sp.simplify(upper_value - 1)
            elif upper_value == sp.oo:
                y_value = sp.simplify(lower_value + 1)
            else:
                y_value = sp.simplify((lower_value + upper_value) / 2)
        return {self.x_variable: x_value, self.y_variable: y_value}


def _normalize_formula(condition: object) -> sp.Expr:
    """Normalize a geometry condition using the shared public-API rules."""

    return normalize_formula(condition)


def _normalize_variables(
    variables: Sequence[sp.Symbol | str] | None,
    formula: sp.Expr | None = None,
) -> tuple[sp.Symbol, ...]:
    """Resolve geometry variables and append symbols from the formula."""

    context = () if formula is None else (formula,)
    return _shared_variables(variables, *context)


def _as_leq_expression(expr: sp.Expr, op: str) -> sp.Expr | None:
    if op in ("<", "<="):
        return sp.simplify(expr)
    if op in (">", ">="):
        return sp.simplify(-expr)
    if op == "!=":
        return None
    raise ValueError(op)


def _formula_to_nnf(expr: sp.Expr) -> sp.Expr:
    """Push negations to atoms for relational Boolean formulas."""

    if expr is sp.true or isinstance(expr, BooleanTrue):
        return sp.true
    if expr is sp.false or isinstance(expr, BooleanFalse):
        return sp.false
    if isinstance(expr, SymNot):
        arg = expr.args[0]
        if isinstance(arg, SymAnd):
            return sp.Or(*[_formula_to_nnf(sp.Not(part)) for part in arg.args])
        if isinstance(arg, SymOr):
            return sp.And(*[_formula_to_nnf(sp.Not(part)) for part in arg.args])
        if isinstance(arg, SymNot):
            return _formula_to_nnf(arg.args[0])
        if getattr(arg, "is_Relational", False):
            return arg.negated
        return sp.Not(_formula_to_nnf(arg))
    if isinstance(expr, SymAnd):
        return sp.And(*[_formula_to_nnf(arg) for arg in expr.args])
    if isinstance(expr, SymOr):
        return sp.Or(*[_formula_to_nnf(arg) for arg in expr.args])
    return expr


def _distribute_or_over_and(expr: sp.Expr) -> sp.Expr:
    expr = _formula_to_nnf(expr)
    if isinstance(expr, SymAnd):
        args = [_distribute_or_over_and(arg) for arg in expr.args]
        ors = [arg for arg in args if isinstance(arg, SymOr)]
        if not ors:
            return sp.And(*args)
        first_or = ors[0]
        rest: list[sp.Expr] = []
        used = False
        for arg in args:
            if arg is first_or and not used:
                used = True
                continue
            rest.append(arg)
        return sp.Or(*[_distribute_or_over_and(sp.And(option, *rest)) for option in first_or.args])
    if isinstance(expr, SymOr):
        return sp.Or(*[_distribute_or_over_and(arg) for arg in expr.args])
    return expr


def _as_disjuncts(expr: sp.Expr) -> tuple[sp.Expr, ...]:
    if expr is sp.false or isinstance(expr, BooleanFalse):
        return ()
    dnf = _distribute_or_over_and(expr)
    if dnf is sp.true or isinstance(dnf, BooleanTrue):
        return (sp.true,)
    if isinstance(dnf, SymOr):
        return tuple(dnf.args)
    return (dnf,)


def _as_conjuncts(expr: sp.Expr) -> tuple[sp.Expr, ...]:
    if expr is sp.true or isinstance(expr, BooleanTrue):
        return ()
    if isinstance(expr, SymAnd):
        return tuple(expr.args)
    return (expr,)


def semialgebraic_level_function(
    condition: object, variables: Sequence[sp.Symbol | str] | None = None
) -> sp.Expr:
    """Return a scalar level function whose sublevel set represents ``condition``.

    The returned expression follows the standard convention: conjunctions map
    to ``Max`` and disjunctions map to ``Min``. A relation ``a <= b`` maps to
    ``a - b``; strict inequalities are relaxed to their weak counterparts.
    Equalities are intentionally unsupported because a single ``f <= 0`` level
    set cannot represent an arbitrary equality without changing dimension.
    Unequalities are ignored by mapping them to ``0``, matching closure-style
    semialgebraic region utilities.
    """

    expr = _normalize_formula(condition)
    _normalize_variables(variables, expr)  # validates symbols and preserves API symmetry

    def rec(node: sp.Expr) -> sp.Expr:
        if node is sp.true or isinstance(node, BooleanTrue):
            return sp.Integer(0)
        if node is sp.false or isinstance(node, BooleanFalse):
            return sp.oo
        if isinstance(node, SymAnd):
            return sp.Max(*[rec(arg) for arg in node.args])
        if isinstance(node, SymOr):
            return sp.Min(*[rec(arg) for arg in node.args])
        if isinstance(node, SymNot):
            raise NotImplementedError(
                "level functions require formulas with negations pushed to relational atoms"
            )
        if getattr(node, "is_Relational", False):
            diff, op = _relational_difference(node)
            if op == "==":
                raise NotImplementedError(
                    "level functions for equality constraints are not represented as one inequality"
                )
            if op == "!=":
                return sp.Integer(0)
            leq = _as_leq_expression(diff, op)
            assert leq is not None
            return sp.simplify(leq)
        raise TypeError(f"unsupported Boolean expression in level function: {node!r}")

    return sp.simplify(rec(_formula_to_nnf(expr)))


def decompose_implicit_formula(
    condition: object, variables: Sequence[sp.Symbol | str] | None = None
) -> tuple[ImplicitFormulaPiece, ...]:
    """Decompose a Boolean semialgebraic condition into DNF implicit pieces.

    Each disjunctive piece contains relaxed inequalities ``f <= 0`` and
    equalities ``g == 0``. Unequalities are ignored. This is primarily a
    structural utility for closure/boundary/integration code, not a prettifier.
    """

    expr = _normalize_formula(condition)
    _normalize_variables(variables, expr)
    pieces: list[ImplicitFormulaPiece] = []
    for disjunct in _as_disjuncts(expr):
        inequalities: list[sp.Expr] = []
        equalities: list[sp.Expr] = []
        impossible = False
        for atom in _as_conjuncts(disjunct):
            if atom is sp.true or isinstance(atom, BooleanTrue):
                continue
            if atom is sp.false or isinstance(atom, BooleanFalse):
                impossible = True
                break
            if not getattr(atom, "is_Relational", False):
                raise TypeError(f"expected relational atom in implicit decomposition, got {atom!r}")
            diff, op = _relational_difference(atom)
            if op == "==":
                equalities.append(sp.simplify(diff))
            elif op == "!=":
                continue
            else:
                leq = _as_leq_expression(diff, op)
                if leq is not None:
                    inequalities.append(sp.simplify(leq))
        if not impossible:
            pieces.append(ImplicitFormulaPiece(tuple(inequalities), tuple(equalities)))
    return tuple(pieces)


def _merge_min(current: sp.Expr, candidate: sp.Expr) -> sp.Expr:
    if current == sp.oo:
        return sp.simplify(candidate)
    return sp.Min(current, candidate)


def _merge_max(current: sp.Expr, candidate: sp.Expr) -> sp.Expr:
    if current == -sp.oo:
        return sp.simplify(candidate)
    return sp.Max(current, candidate)


def _linear_bound_from_relation(
    atom: sp.Expr, variable: sp.Symbol, later_variables: Sequence[sp.Symbol]
) -> tuple[str, sp.Expr] | None:
    diff, op = _relational_difference(atom)
    if op not in ("<", "<=", ">", ">="):
        return None
    if diff.free_symbols & set(later_variables):
        return None
    if variable not in diff.free_symbols:
        return None
    try:
        poly = sp.Poly(diff, variable)
    except Exception:
        return None
    if poly.degree() != 1:
        return None
    coeff = sp.simplify(poly.coeff_monomial(variable))
    rest = sp.simplify(poly.as_expr() - coeff * variable)
    if coeff == 0:
        return None
    boundary = sp.simplify(-rest / coeff)
    normalized_op = op
    if coeff.is_negative:
        normalized_op = {"<": ">", "<=": ">=", ">": "<", ">=": "<="}[op]
    elif not coeff.is_positive:
        # Keep symbolic-sign handling conservative.
        return None
    if normalized_op in ("<", "<="):
        return "upper", boundary
    if normalized_op in (">", ">="):
        return "lower", boundary
    return None


def extract_symbolic_box_bounds(
    condition: object, variables: Sequence[sp.Symbol | str]
) -> SymbolicBoxBounds:
    """Extract explicit axis-aligned bounds from a box-like conjunction.

    Bounds may depend on earlier variables but not on the variable being bounded
    or later variables. For a pure box, they are constants or parameters.
    Raises ``NotImplementedError`` when the formula is not an explicit box.
    """

    expr = _normalize_formula(condition)
    vars_ = _normalize_variables(variables, expr)
    if isinstance(expr, SymOr):
        raise NotImplementedError(
            "explicit box extraction expects one conjunctive cell, not a disjunction"
        )
    lower = {var: -sp.oo for var in vars_}
    upper = {var: sp.oo for var in vars_}
    used_atoms: set[sp.Expr] = set()
    atoms = _as_conjuncts(_formula_to_nnf(expr))
    for i, var in enumerate(vars_):
        later = vars_[i:]
        for atom in atoms:
            if atom in used_atoms or not getattr(atom, "is_Relational", False):
                continue
            bound = _linear_bound_from_relation(atom, var, later_variables=later[1:])
            if bound is None:
                continue
            kind, value = bound
            if kind == "lower":
                lower[var] = _merge_max(lower[var], value)
            else:
                upper[var] = _merge_min(upper[var], value)
            used_atoms.add(atom)
    # Reject unconsumed nontrivial relations involving variables.
    for atom in atoms:
        if atom in used_atoms:
            continue
        if getattr(atom, "is_Relational", False):
            diff, op = _relational_difference(atom)
            if op == "!=" or not (diff.free_symbols & set(vars_)):
                continue
            raise NotImplementedError("condition is not an explicit axis-aligned/symbolic box")
        if atom not in (sp.true, True):
            raise NotImplementedError("condition is not an explicit axis-aligned/symbolic box")
    return SymbolicBoxBounds(
        tuple((var, sp.simplify(lower[var]), sp.simplify(upper[var])) for var in vars_)
    )


def _truth_at(condition: sp.Expr, subs: Mapping[sp.Symbol, sp.Expr]) -> bool:
    """Evaluate a specialized cell condition with exact real arithmetic."""

    return exact_truth(sp.simplify(condition.subs(subs)))


def _x_intervals_from_condition(
    x_condition: sp.Expr, x: sp.Symbol
) -> tuple[tuple[sp.Expr, sp.Expr], ...]:
    if x_condition is sp.true or isinstance(x_condition, BooleanTrue):
        return ((-sp.oo, sp.oo),)
    if x_condition is sp.false or isinstance(x_condition, BooleanFalse):
        return ()
    atoms = _as_conjuncts(_formula_to_nnf(x_condition))
    cuts: list[sp.Expr] = [-sp.oo, sp.oo]
    for atom in atoms:
        if not getattr(atom, "is_Relational", False):
            raise NotImplementedError("x-cell conditions must be relational")
        diff, _ = _relational_difference(atom)
        if diff.free_symbols - {x}:
            raise NotImplementedError("x-cell conditions must be univariate in the base variable")
        if x in diff.free_symbols:
            cuts.extend(_finite_real_roots(diff, x))
    ordered: list[sp.Expr] = []
    seen: set[str] = set()

    def compare_cuts(left: sp.Expr, right: sp.Expr) -> int:
        if left == right:
            return 0
        if left == -sp.oo or right == sp.oo:
            return -1
        if left == sp.oo or right == -sp.oo:
            return 1
        return compare_exact_reals(left, right)

    for cut in sorted(cuts, key=cmp_to_key(compare_cuts)):
        key = sp.sstr(cut)
        if key not in seen:
            ordered.append(cut)
            seen.add(key)
    intervals: list[tuple[sp.Expr, sp.Expr]] = []
    for lo, hi in zip(ordered, ordered[1:], strict=False):
        if lo == hi:
            continue
        sample = (
            sp.Integer(0)
            if lo == -sp.oo and hi == sp.oo
            else (
                hi - 1 if lo == -sp.oo else (lo + 1 if hi == sp.oo else sp.simplify((lo + hi) / 2))
            )
        )
        if _truth_at(x_condition, {x: sample}):
            intervals.append((lo, hi))
    # Point cells are useful for boundary but not ambient integration.
    for cut in ordered:
        if cut in (-sp.oo, sp.oo):
            continue
        if _truth_at(x_condition, {x: cut}):
            intervals.append((cut, cut))
    return tuple(sorted(intervals, key=cmp_to_key(lambda a, b: compare_cuts(a[0], b[0]))))


def _vertical_bound_from_atom(
    atom: sp.Expr, x: sp.Symbol, y: sp.Symbol
) -> tuple[str, sp.Expr] | sp.Expr | None:
    if not getattr(atom, "is_Relational", False):
        return None
    diff, op = _relational_difference(atom)
    if op == "!=":
        return sp.true
    if y not in diff.free_symbols:
        return make_zero_relation(diff, op)
    try:
        poly_y = sp.Poly(diff, y)
    except Exception:
        return None
    if poly_y.degree() == 2:
        a2 = sp.simplify(poly_y.coeff_monomial(y**2))
        a1 = sp.simplify(poly_y.coeff_monomial(y))
        rest2 = sp.simplify(poly_y.as_expr() - a2 * y**2 - a1 * y)
        if a1 == 0 and a2.is_positive and op in ("<", "<="):
            radicand = sp.simplify(-rest2 / a2)
            if not (radicand.free_symbols - {x}):
                root = sp.sqrt(radicand)
                return "between", sp.simplify(-root), sp.simplify(root), sp.simplify(radicand >= 0)
        return None
    if poly_y.degree() != 1:
        return None
    coeff = sp.simplify(poly_y.coeff_monomial(y))
    rest = sp.simplify(poly_y.as_expr() - coeff * y)
    bound = sp.simplify(-rest / coeff)
    if bound.free_symbols - {x}:
        return None
    normalized_op = op
    if coeff.is_negative:
        normalized_op = {"<": ">", "<=": ">=", ">": "<", ">=": "<="}.get(op, op)
    elif not coeff.is_positive:
        return None
    if normalized_op == "==":
        return "equal", bound
    if normalized_op in ("<", "<="):
        return "upper", bound
    if normalized_op in (">", ">="):
        return "lower", bound
    return None


def decompose_cylindrical_formula_to_vertical_bounds_2d(
    condition: object,
    variables: Sequence[sp.Symbol | str],
) -> tuple[VerticalBoundCell2D, ...]:
    """Convert supported 2D cylindrical formulas into vertical bounds.

    The accepted formulas are disjunctions of cells whose constraints are
    cylindrical in ``(x, y)``: base conditions in ``x`` and linear-in-``y``
    stack constraints such as ``g(x) <= y <= h(x)`` or ``Eq(y, g(x))``.
    This is intentionally a descriptive, public helper for a 2D CAD stack
    extraction utility; unsupported formulas raise ``NotImplementedError``.
    """

    expr = _normalize_formula(condition)
    vars_ = _normalize_variables(variables, expr)
    if len(vars_) != 2:
        raise ValueError("vertical-bound decomposition requires exactly two variables")
    x, y = vars_
    cells: list[VerticalBoundCell2D] = []
    for disjunct in _as_disjuncts(expr):
        lower: list[sp.Expr] = []
        upper: list[sp.Expr] = []
        equal: list[sp.Expr] = []
        x_conditions: list[sp.Expr] = []
        for atom in _as_conjuncts(disjunct):
            if atom is sp.true or isinstance(atom, BooleanTrue):
                continue
            parsed = _vertical_bound_from_atom(atom, x, y)
            if parsed is None:
                raise NotImplementedError(
                    f"formula is not a supported cylindrical 2D cell: {atom!r}"
                )
            if isinstance(parsed, tuple):
                kind = parsed[0]
                if kind == "lower":
                    lower.append(parsed[1])
                elif kind == "upper":
                    upper.append(parsed[1])
                elif kind == "equal":
                    equal.append(parsed[1])
                elif kind == "between":
                    _, lower_value, upper_value, condition_value = parsed
                    lower.append(lower_value)
                    upper.append(upper_value)
                    x_conditions.append(condition_value)
            else:
                x_conditions.append(parsed)
        if equal:
            if len(equal) > 1 and any(sp.simplify(item - equal[0]) != 0 for item in equal[1:]):
                continue
            y_bounds = ((sp.simplify(equal[0]), sp.simplify(equal[0])),)
        else:
            lo = sp.Max(*lower) if len(lower) > 1 else (lower[0] if lower else -sp.oo)
            hi = sp.Min(*upper) if len(upper) > 1 else (upper[0] if upper else sp.oo)
            if lo == -sp.oo or hi == sp.oo:
                raise NotImplementedError(
                    "unbounded vertical cells are not supported by this reduction"
                )
            width_condition = sp.simplify(lo <= hi)
            try:
                _ = (
                    sp.Poly(sp.expand(width_condition.lhs - width_condition.rhs), x)
                    if getattr(width_condition, "is_Relational", False)
                    else None
                )
            except Exception:
                width_condition = sp.true
            if width_condition is not sp.true and not isinstance(width_condition, BooleanTrue):
                x_conditions.append(width_condition)
            y_bounds = ((sp.simplify(lo), sp.simplify(hi)),)
        x_condition = sp.And(*x_conditions) if x_conditions else sp.true
        for interval in _x_intervals_from_condition(x_condition, x):
            cells.append(
                VerticalBoundCell2D(
                    x_variable=x,
                    y_variable=y,
                    x_interval=interval,
                    y_bounds=y_bounds,
                    x_condition=x_condition,
                    source_formula=disjunct,
                )
            )
    return tuple(cells)


__all__ = [
    "ImplicitFormulaPiece",
    "SymbolicBoxBounds",
    "VerticalBoundCell2D",
    "semialgebraic_level_function",
    "decompose_implicit_formula",
    "extract_symbolic_box_bounds",
    "decompose_cylindrical_formula_to_vertical_bounds_2d",
]
