from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field

import sympy as sp
from sympy.core.relational import Relational
from sympy.logic.boolalg import BooleanFalse, BooleanTrue

from .._zero_testing import certified_zero
from ..exact_arithmetic import compare_exact_reals
from ..structural_keys import symbol_identity_key

_RECOVERABLE_ERRORS = (
    ArithmeticError,
    TypeError,
    ValueError,
    NotImplementedError,
    sp.PolynomialError,
)


@dataclass(frozen=True)
class FallbackAttempt:
    """One attempted non-CAD real-instance method."""

    name: str
    status: str
    details: Mapping[str, object] = field(default_factory=dict)


@dataclass(frozen=True)
class FallbackInstanceResult:
    """Result of the lightweight real algebraic fallback pipeline."""

    instances: tuple[dict[sp.Symbol, sp.Expr], ...]
    status: str
    method: str
    attempts: tuple[FallbackAttempt, ...]
    exact: bool = True

    @property
    def found(self) -> bool:
        return bool(self.instances)

    def first(self) -> dict[sp.Symbol, sp.Expr] | None:
        return self.instances[0] if self.instances else None


@dataclass(frozen=True)
class CoordinateBounds:
    """Conservative coordinate bounds implied by a formula.

    ``bounds`` stores ``(variable, lower, upper)`` triples. ``strictness``
    records whether each finite endpoint is open
    as ``(variable, lower_strict, upper_strict)``.
    """

    bounds: tuple[tuple[sp.Symbol, sp.Expr, sp.Expr], ...]
    inconsistent: bool = False
    complete: bool = False
    strictness: tuple[tuple[sp.Symbol, bool, bool], ...] = ()


@dataclass(frozen=True)
class LinearEliminationResult:
    """Result of safe linear equation elimination."""

    equations: tuple[sp.Expr, ...]
    variables: tuple[sp.Symbol, ...]
    replacements: tuple[tuple[sp.Symbol, sp.Expr], ...]


@dataclass(frozen=True)
class MethodSearchResult:
    """Outcome of a staged method search."""

    value: object
    status: str
    partial_results: tuple[object, ...] = ()
    successful_method: str | None = None


# ---------------------------------------------------------------------------
# Numeric and symbolic utility predicates


def is_valid_numeric_value(value: object) -> bool:
    """Return True for finite usable numeric values.

    This rejects NaN, complex infinities, and symbolic infinities. It is used
    by the permissive validator and bounded sampling layer before a candidate is
    trusted as a witness.
    """

    try:
        value = sp.sympify(value)
    except _RECOVERABLE_ERRORS:
        return False
    if value in {sp.oo, -sp.oo, sp.zoo, sp.nan}:
        return False
    if value.has(sp.oo, -sp.oo, sp.zoo, sp.nan):
        return False
    if value.is_number:
        # This helper is used by the real-instance pipeline: finite complex
        # numbers are not valid real witnesses. Unknown realness is rejected
        # conservatively rather than accepted from a floating conversion.
        if value.is_real is not True:
            return False
        try:
            numeric = value.evalf(30)
            return bool(numeric.is_finite is True)
        except _RECOVERABLE_ERRORS:
            return False
    return False


def is_reliably_zero(expr: object) -> bool:
    """Conservative exact/algebraic zero test."""

    try:
        value = sp.sympify(expr)
    except _RECOVERABLE_ERRORS:
        return False
    if value == 0:
        return True
    try:
        simplified = sp.simplify(value)
        if simplified == 0:
            return True
    except _RECOVERABLE_ERRORS:
        pass
    try:
        return certified_zero(value) is True
    except _RECOVERABLE_ERRORS:
        return False


def could_be_zero(expr: object) -> bool:
    """Return True when a numeric expression is zero or too close to trust."""

    try:
        value = sp.sympify(expr)
    except _RECOVERABLE_ERRORS:
        return True
    if is_reliably_zero(value):
        return True
    if not value.is_number:
        return True
    try:
        approx = complex(value.evalf(30))
    except _RECOVERABLE_ERRORS:
        return True
    return abs(approx) < 1.0e-12


def _truth_value(value: object) -> bool | None:
    if value is True or value is sp.true or isinstance(value, BooleanTrue):
        return True
    if value is False or value is sp.false or isinstance(value, BooleanFalse):
        return False
    try:
        if value == True:  # noqa: E712 - intentional SymPy coercion point
            return True
        if value == False:  # noqa: E712 - intentional SymPy coercion point
            return False
    except _RECOVERABLE_ERRORS:
        pass
    return None


def _safe_simplify(expr: sp.Expr) -> sp.Expr:
    try:
        return sp.simplify(expr)
    except _RECOVERABLE_ERRORS:
        return expr


def _relation_delta(rel: Relational) -> sp.Expr:
    return sp.expand(rel.lhs - rel.rhs)


def _atoms(expr: sp.Expr) -> tuple[sp.Expr, ...]:
    if expr is sp.true or expr is True:
        return (sp.true,)
    if expr is sp.false or expr is False:
        return (sp.false,)
    if isinstance(expr, sp.And):
        return tuple(expr.args)
    return (expr,)


def _relations(expr: sp.Expr) -> tuple[Relational, ...]:
    return tuple(atom for atom in _atoms(expr) if isinstance(atom, Relational))


# ---------------------------------------------------------------------------
# Formula normalization and candidate validation


def _eval_atom(
    atom: sp.Expr, assignment: Mapping[sp.Symbol, object], *, strict: bool
) -> bool | None:
    try:
        value = atom.subs(assignment)
    except _RECOVERABLE_ERRORS:
        return None
    value = _safe_simplify(value)
    truth = _truth_value(value)
    if truth is not None:
        return truth
    if not strict:
        if not hasattr(value, "evalf"):
            return None
        try:
            numeric = value.evalf(50)
            truth = _truth_value(numeric)
            if truth is not None:
                return truth
        except _RECOVERABLE_ERRORS:
            pass
    return None


def satisfies_formula(
    formula: sp.Expr,
    assignment: Mapping[sp.Symbol, object],
    *,
    strict: bool = True,
) -> bool:
    """Return whether ``assignment`` satisfies ``formula`` over the reals.

    Strict mode rejects unresolved atoms. Permissive mode still prefers exact
    truth values but accepts numerical truth values when SymPy can evaluate
    them cleanly.
    """

    if isinstance(formula, sp.Or):
        return any(satisfies_formula(arg, assignment, strict=strict) for arg in formula.args)
    for atom in _atoms(formula):
        truth = _eval_atom(atom, assignment, strict=strict)
        if truth is not True:
            return False
    return True


def relation_to_zero_rhs(atom: Relational) -> Relational:
    """Normalize a binary relation to a zero right-hand side."""

    delta = _relation_delta(atom)
    if isinstance(atom, sp.Equality):
        return sp.Eq(delta, 0)
    if isinstance(atom, sp.Unequality):
        return sp.Ne(delta, 0)
    if isinstance(atom, sp.StrictLessThan):
        return sp.Lt(delta, 0)
    if isinstance(atom, sp.LessThan):
        return sp.Le(delta, 0)
    if isinstance(atom, sp.StrictGreaterThan):
        return sp.Lt(-delta, 0)
    if isinstance(atom, sp.GreaterThan):
        return sp.Le(-delta, 0)
    return atom


def relations_to_zero_rhs(formula: sp.Expr) -> sp.Expr:
    """Normalize relation atoms in a Boolean formula."""

    if isinstance(formula, sp.And):
        return sp.And(*(relations_to_zero_rhs(arg) for arg in formula.args))
    if isinstance(formula, sp.Or):
        return sp.Or(*(relations_to_zero_rhs(arg) for arg in formula.args))
    if isinstance(formula, Relational):
        return relation_to_zero_rhs(formula)
    return formula


def to_dnf_formula(formula: sp.Expr) -> sp.Expr:
    """Distribute conjunction over disjunction while preserving SymPy atoms."""

    try:
        return sp.to_dnf(formula, simplify=False)
    except _RECOVERABLE_ERRORS:
        return formula


def _apply_equalities_to_later_atoms(args: Sequence[sp.Expr]) -> tuple[sp.Expr, ...]:
    rules: list[tuple[sp.Expr, sp.Expr]] = []
    for atom in args:
        if isinstance(atom, sp.Equality) and not atom.lhs.is_number:
            rules.append((atom.lhs, atom.rhs))
    repl = dict(rules)
    out: list[sp.Expr] = []
    for atom in args:
        current = atom
        if repl and not isinstance(atom, sp.Equality):
            try:
                current = current.xreplace(repl)
            except _RECOVERABLE_ERRORS:
                try:
                    current = current.subs(repl)
                except _RECOVERABLE_ERRORS:
                    pass
        out.append(current)
    return tuple(out)


def expand_with_subs(formula: sp.Expr) -> sp.Expr:
    """Return a DNF-like formula after propagating earlier equations.

    This adapts a useful formula-normalization idea:
    equations that bind a nonnumeric left-hand side are applied to later atoms
    in the same conjunction before distributing disjunctions.
    """

    formula = to_dnf_formula(formula)
    if isinstance(formula, sp.Or):
        return sp.Or(*(expand_with_subs(arg) for arg in formula.args))
    if isinstance(formula, sp.And):
        return sp.And(*_apply_equalities_to_later_atoms(formula.args))
    return formula


def formula_to_rule_sets(
    formula: sp.Expr,
) -> tuple[tuple[tuple[sp.Expr, sp.Expr], ...], ...] | None:
    """Convert a pure equality formula to replacement rule sets when possible."""

    formula = to_dnf_formula(formula)
    branches = formula.args if isinstance(formula, sp.Or) else (formula,)
    all_rules: list[tuple[tuple[sp.Expr, sp.Expr], ...]] = []
    for branch in branches:
        atoms = branch.args if isinstance(branch, sp.And) else (branch,)
        rules: list[tuple[sp.Expr, sp.Expr]] = []
        for atom in atoms:
            if atom in {sp.true, True}:
                continue
            if not isinstance(atom, sp.Equality) or atom.lhs.is_number:
                return None
            rules.append((atom.lhs, atom.rhs))
        all_rules.append(tuple(rules))
    return tuple(all_rules)


# ---------------------------------------------------------------------------
# Rationalization, factors, algebraic variables, and ordering helpers


def is_rational_number(value: object) -> bool:
    try:
        value = sp.sympify(value)
    except _RECOVERABLE_ERRORS:
        return False
    return value.is_Integer is True or value.is_Rational is True


def as_rational_if_exact(value: object, *, max_denominator: int = 10_000) -> sp.Rational | None:
    """Heuristically recognize a value as a small rational number."""

    try:
        value = sp.sympify(value)
    except _RECOVERABLE_ERRORS:
        return None
    if value.is_Rational:
        return sp.Rational(value)
    try:
        rational = sp.Rational(str(value.evalf(30))).limit_denominator(max_denominator)
        if is_reliably_zero(value - rational):
            return rational
    except _RECOVERABLE_ERRORS:
        return None
    return None


def rational_bound(value: object, direction: int) -> sp.Expr:
    """Return a rational approximation below or above a real value.

    ``direction=-1`` requests a lower rational bound and ``direction=1`` an
    upper rational bound. Infinities are returned when no trustworthy rational
    approximation is available.
    """

    if direction not in {-1, 1}:
        raise ValueError("direction must be -1 or 1")
    try:
        value = sp.sympify(value)
    except _RECOVERABLE_ERRORS:
        return -sp.oo if direction < 0 else sp.oo
    if value in {sp.oo, -sp.oo} or value.is_Rational:
        return value
    if value.is_real is False:
        return -sp.oo if direction < 0 else sp.oo
    # Approximation proposes a compact rational candidate; exact comparison
    # certifies and, if necessary, adjusts it.  Fixed-precision truth values
    # never decide whether the result is actually a lower/upper bound.
    try:
        approx = value.evalf(40)
        if approx.is_finite is False:
            return -sp.oo if direction < 0 else sp.oo
        base = sp.Rational(str(approx)).limit_denominator(10_000)
        cmp = compare_exact_reals(base, value)
    except _RECOVERABLE_ERRORS:
        return -sp.oo if direction < 0 else sp.oo
    if direction < 0 and cmp <= 0:
        return base
    if direction > 0 and cmp >= 0:
        return base

    # A single denominator-sized step is normally sufficient; retain exact
    # certification in the loop for values exceptionally close to a candidate.
    step = sp.Rational(1, 10_000)
    for _ in range(20_000):
        base = base - step if direction < 0 else base + step
        cmp = compare_exact_reals(base, value)
        if direction < 0 and cmp <= 0:
            return base
        if direction > 0 and cmp >= 0:
            return base
    return -sp.oo if direction < 0 else sp.oo


def algebraic_variables(expr: sp.Expr) -> tuple[sp.Expr, ...]:
    """Return variables appearing at the algebraic level of an expression."""

    expr = sp.sympify(expr)
    ignored_heads = (sp.Add, sp.Mul, sp.Pow, Relational, sp.And, sp.Or)
    if expr.is_number:
        return ()
    if isinstance(expr, Relational):
        found: set[sp.Expr] = set()
        for side in (expr.lhs, expr.rhs):
            found.update(algebraic_variables(side))
        return tuple(
            sorted(
                found,
                key=lambda item: (
                    symbol_identity_key(item)
                    if isinstance(item, sp.Symbol)
                    else (sp.sstr(item), sp.srepr(item))
                ),
            )
        )
    if isinstance(expr, (sp.And, sp.Or)):
        found: set[sp.Expr] = set()
        for arg in expr.args:
            found.update(algebraic_variables(arg))
        return tuple(
            sorted(
                found,
                key=lambda item: (
                    symbol_identity_key(item)
                    if isinstance(item, sp.Symbol)
                    else (sp.sstr(item), sp.srepr(item))
                ),
            )
        )
    if isinstance(expr, sp.Pow) and expr.exp.is_Rational:
        return algebraic_variables(expr.base)
    if isinstance(expr, ignored_heads):
        found: set[sp.Expr] = set()
        for arg in expr.args:
            found.update(algebraic_variables(arg))
        return tuple(
            sorted(
                found,
                key=lambda item: (
                    symbol_identity_key(item)
                    if isinstance(item, sp.Symbol)
                    else (sp.sstr(item), sp.srepr(item))
                ),
            )
        )
    if expr.is_Symbol:
        return (expr,)
    return tuple(sorted(expr.free_symbols, key=symbol_identity_key))


def is_algebraic_condition(condition: sp.Expr, variables: Sequence[sp.Symbol]) -> bool:
    """Return whether a relation is algebraic in the requested variables."""

    allowed = set(variables)
    return set(algebraic_variables(condition)).issubset(allowed)


def fast_factor_list(
    expr: sp.Expr, *, max_power_cost: int = 100
) -> tuple[tuple[sp.Expr, int], ...]:
    """Return factor pairs without forcing expensive expansion.

    Expressions with very high powers are decomposed structurally, avoiding a
    potentially explosive call to ``factor_list``.
    """

    expr = sp.sympify(expr)

    def power_cost(term: sp.Expr) -> int:
        if isinstance(term, sp.Pow) and term.exp.is_Integer:
            return abs(int(term.exp)) * power_cost(term.base)
        if isinstance(term, sp.Mul):
            return sum(power_cost(arg) for arg in term.args)
        if isinstance(term, sp.Add):
            return max((power_cost(arg) for arg in term.args), default=1)
        return 1

    if power_cost(expr) <= max_power_cost:
        try:
            coeff, factors = sp.factor_list(expr)
            out = [] if coeff == 1 else [(sp.sympify(coeff), 1)]
            out.extend((sp.sympify(factor), int(exp)) for factor, exp in factors)
            return tuple(out)
        except _RECOVERABLE_ERRORS:
            pass
    if isinstance(expr, sp.Mul):
        out: list[tuple[sp.Expr, int]] = []
        for arg in expr.args:
            out.extend(fast_factor_list(arg, max_power_cost=max_power_cost))
        return tuple(out)
    if isinstance(expr, sp.Pow) and expr.exp.is_Integer:
        return ((expr.base, int(expr.exp)),)
    return ((expr, 1),)


def sort_conditions(formula: sp.Expr, variables: Sequence[sp.Symbol]) -> sp.Expr:
    """Sort conjunction atoms into a stable, variable-aware order."""

    if isinstance(formula, sp.Or):
        return sp.Or(*(sort_conditions(arg, variables) for arg in formula.args))
    if not isinstance(formula, sp.And):
        return formula
    positions = {var: idx for idx, var in enumerate(variables)}

    def key(atom: sp.Expr) -> tuple[int, int, str]:
        if isinstance(atom, Relational):
            syms = sorted(
                atom.free_symbols, key=lambda s: (positions.get(s, 10_000), symbol_identity_key(s))
            )
            pos = positions.get(syms[0], 10_000) if syms else 10_000
            return (2, pos, sp.sstr(atom))
        if isinstance(atom, sp.Or):
            syms = sorted(
                atom.free_symbols, key=lambda s: (positions.get(s, 10_000), symbol_identity_key(s))
            )
            pos = positions.get(syms[0], 10_000) if syms else 10_000
            return (1, pos, sp.sstr(atom))
        return (0, 10_000, sp.sstr(atom))

    return sp.And(*sorted(formula.args, key=key))


__all__ = [
    "FallbackAttempt",
    "FallbackInstanceResult",
    "CoordinateBounds",
    "LinearEliminationResult",
    "MethodSearchResult",
    "is_valid_numeric_value",
    "is_reliably_zero",
    "could_be_zero",
    "satisfies_formula",
    "relation_to_zero_rhs",
    "relations_to_zero_rhs",
    "to_dnf_formula",
    "expand_with_subs",
    "formula_to_rule_sets",
    "is_rational_number",
    "as_rational_if_exact",
    "rational_bound",
    "algebraic_variables",
    "is_algebraic_condition",
    "fast_factor_list",
    "sort_conditions",
]
