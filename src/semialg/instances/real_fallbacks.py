from __future__ import annotations

import math
import signal
from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from itertools import product

import sympy as sp
from sympy.core.relational import Relational
from sympy.logic.boolalg import BooleanFalse, BooleanTrue

from .random_sections import find_random_section_wit
from .witness_generation import sample_free_assignments


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
    """Conservative coordinate bounds implied by a formula."""

    bounds: tuple[tuple[sp.Symbol, sp.Expr, sp.Expr], ...]
    inconsistent: bool = False
    complete: bool = False


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
    except Exception:
        return False
    if value in {sp.oo, -sp.oo, sp.zoo, sp.nan}:
        return False
    if value.has(sp.oo, -sp.oo, sp.zoo, sp.nan):
        return False
    if value.is_number:
        try:
            complex(value.evalf(30))
            return True
        except Exception:
            return False
    return False


def is_reliably_zero(expr: object) -> bool:
    """Conservative exact/algebraic zero test."""

    try:
        value = sp.sympify(expr)
    except Exception:
        return False
    if value == 0:
        return True
    try:
        simplified = sp.simplify(value)
        if simplified == 0:
            return True
    except Exception:
        pass
    try:
        return bool(value.equals(0))
    except Exception:
        return False


def could_be_zero(expr: object) -> bool:
    """Return True when a numeric expression is zero or too close to trust."""

    try:
        value = sp.sympify(expr)
    except Exception:
        return True
    if is_reliably_zero(value):
        return True
    if not value.is_number:
        return True
    try:
        approx = complex(value.evalf(30))
    except Exception:
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
    except Exception:
        pass
    return None


def _safe_simplify(expr: sp.Expr) -> sp.Expr:
    try:
        return sp.simplify(expr)
    except Exception:
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
    except Exception:
        return None
    value = _safe_simplify(value)
    truth = _truth_value(value)
    if truth is not None:
        return truth
    if not strict:
        try:
            numeric = value.evalf(50)
            truth = _truth_value(numeric)
            if truth is not None:
                return truth
        except Exception:
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


def normalize_relation(atom: Relational) -> Relational:
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


def normalize_relations(formula: sp.Expr) -> sp.Expr:
    """Normalize relation atoms in a Boolean formula."""

    if isinstance(formula, sp.And):
        return sp.And(*(normalize_relations(arg) for arg in formula.args))
    if isinstance(formula, sp.Or):
        return sp.Or(*(normalize_relations(arg) for arg in formula.args))
    if isinstance(formula, Relational):
        return normalize_relation(formula)
    return formula


def to_dnf_formula(formula: sp.Expr) -> sp.Expr:
    """Distribute conjunction over disjunction while preserving SymPy atoms."""

    try:
        return sp.to_dnf(formula, simplify=False)
    except Exception:
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
            except Exception:
                try:
                    current = current.subs(repl)
                except Exception:
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
    except Exception:
        return False
    return value.is_Integer is True or value.is_Rational is True


def as_rational_if_exact(value: object, *, max_denominator: int = 10_000) -> sp.Rational | None:
    """Heuristically recognize a value as a small rational number."""

    try:
        value = sp.sympify(value)
    except Exception:
        return None
    if value.is_Rational:
        return sp.Rational(value)
    try:
        rational = sp.Rational(str(value.evalf(30))).limit_denominator(max_denominator)
        if is_reliably_zero(value - rational):
            return rational
    except Exception:
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
    except Exception:
        return -sp.oo if direction < 0 else sp.oo
    if value in {sp.oo, -sp.oo} or value.is_Rational:
        return value
    if value.is_real is False:
        return -sp.oo if direction < 0 else sp.oo
    try:
        approx = float(value.evalf(30))
    except Exception:
        return -sp.oo if direction < 0 else sp.oo
    if not math.isfinite(approx):
        return -sp.oo if direction < 0 else sp.oo
    base = sp.Rational(approx).limit_denominator(10_000)
    if direction < 0:
        while not bool(sp.N(base <= value, 30)):
            base -= sp.Rational(1, 10_000)
    else:
        while not bool(sp.N(base >= value, 30)):
            base += sp.Rational(1, 10_000)
    return base


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
        return tuple(sorted(found, key=sp.sstr))
    if isinstance(expr, (sp.And, sp.Or)):
        found: set[sp.Expr] = set()
        for arg in expr.args:
            found.update(algebraic_variables(arg))
        return tuple(sorted(found, key=sp.sstr))
    if isinstance(expr, sp.Pow) and expr.exp.is_Rational:
        return algebraic_variables(expr.base)
    if isinstance(expr, ignored_heads):
        found: set[sp.Expr] = set()
        for arg in expr.args:
            found.update(algebraic_variables(arg))
        return tuple(sorted(found, key=sp.sstr))
    if expr.is_Symbol:
        return (expr,)
    return tuple(sorted(expr.free_symbols, key=sp.sstr))


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
        except Exception:
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
            syms = sorted(atom.free_symbols, key=lambda s: positions.get(s, 10_000))
            pos = positions.get(syms[0], 10_000) if syms else 10_000
            return (2, pos, sp.sstr(atom))
        if isinstance(atom, sp.Or):
            syms = sorted(atom.free_symbols, key=lambda s: positions.get(s, 10_000))
            pos = positions.get(syms[0], 10_000) if syms else 10_000
            return (1, pos, sp.sstr(atom))
        return (0, 10_000, sp.sstr(atom))

    return sp.And(*sorted(formula.args, key=key))


# ---------------------------------------------------------------------------
# Bounds, boundedness, vector inequalities, and linear elimination


def _linear_coeffs(
    expr: sp.Expr, variables: Sequence[sp.Symbol]
) -> tuple[list[sp.Expr], sp.Expr] | None:
    expr = sp.expand(expr)
    coeffs: list[sp.Expr] = []
    remainder = expr
    for var in variables:
        coeff = sp.diff(expr, var)
        if coeff.free_symbols & set(variables):
            return None
        coeffs.append(coeff)
        remainder -= coeff * var
    if remainder.free_symbols & set(variables):
        return None
    return coeffs, sp.simplify(remainder)


def _relation_to_linear_bound(
    rel: Relational, variables: Sequence[sp.Symbol]
) -> tuple[sp.Symbol, sp.Expr, str] | None:
    norm = normalize_relation(rel)
    if not isinstance(norm, (sp.StrictLessThan, sp.LessThan, sp.StrictGreaterThan, sp.GreaterThan)):
        return None
    delta = _relation_delta(norm)
    data = _linear_coeffs(delta, variables)
    if data is None:
        return None
    coeffs, const = data
    nonzero = [(variables[i], coeff) for i, coeff in enumerate(coeffs) if coeff != 0]
    if len(nonzero) != 1:
        return None
    var, coeff = nonzero[0]
    bound = sp.simplify(-const / coeff)
    upper = isinstance(norm, (sp.StrictLessThan, sp.LessThan))
    if coeff.could_extract_minus_sign():
        upper = not upper
    return var, bound, "upper" if upper else "lower"


def coordinate_bounds(
    formula: sp.Expr,
    variables: Sequence[sp.Symbol],
) -> CoordinateBounds:
    """Extract inexpensive coordinate bounds from linear inequalities."""

    variables = tuple(variables)
    lower: dict[sp.Symbol, sp.Expr] = {var: -sp.oo for var in variables}
    upper: dict[sp.Symbol, sp.Expr] = {var: sp.oo for var in variables}
    for atom in _relations(normalize_relations(formula)):
        bound = _relation_to_linear_bound(atom, variables)
        if bound is None:
            continue
        var, value, side = bound
        if side == "upper":
            try:
                if upper[var] is sp.oo or bool(value < upper[var]):
                    upper[var] = value
            except Exception:
                upper[var] = value
        else:
            try:
                if lower[var] is -sp.oo or bool(value > lower[var]):
                    lower[var] = value
            except Exception:
                lower[var] = value
    inconsistent = False
    for var in variables:
        try:
            if bool(lower[var] > upper[var]):
                inconsistent = True
        except Exception:
            pass
    return CoordinateBounds(
        tuple((var, lower[var], upper[var]) for var in variables),
        inconsistent=inconsistent,
        complete=True,
    )


def is_bounded_solution_set(formula: sp.Expr, variables: Sequence[sp.Symbol]) -> bool | None:
    """Try to prove coordinate-boundedness or obvious unboundedness."""

    bounds = coordinate_bounds(formula, variables)
    if bounds.inconsistent:
        return True
    finite_all = all(lo is not -sp.oo and hi is not sp.oo for _, lo, hi in bounds.bounds)
    if finite_all:
        return True
    if formula in {sp.true, True}:
        return False
    return None


def vector_relations(lhs: Sequence[sp.Expr], rhs: Sequence[sp.Expr], relation: str) -> sp.Expr:
    """Convert componentwise vector inequalities to scalar inequalities."""

    if len(lhs) != len(rhs):
        raise ValueError("vector inequality sides must have the same length")
    if relation == "lt":
        return sp.And(*(sp.Lt(a, b) for a, b in zip(lhs, rhs, strict=True)))
    if relation == "le":
        return sp.And(*(sp.Le(a, b) for a, b in zip(lhs, rhs, strict=True)))
    if relation == "gt":
        return sp.And(*(sp.Gt(a, b) for a, b in zip(lhs, rhs, strict=True)))
    if relation == "ge":
        return sp.And(*(sp.Ge(a, b) for a, b in zip(lhs, rhs, strict=True)))
    raise ValueError("relation must be one of 'lt', 'le', 'gt', or 'ge'")


def eliminate_linear_equations(
    equations: Sequence[sp.Expr],
    variables: Sequence[sp.Symbol],
    reference_assignment: Mapping[sp.Symbol, object] | None = None,
) -> LinearEliminationResult:
    """Eliminate variables using equations that are linear in one variable.

    A variable is eliminated only when a relation can be solved with a
    coefficient that is provably nonzero at the optional reference assignment.
    """

    remaining = [
        sp.expand(eq.lhs - eq.rhs) if isinstance(eq, sp.Equality) else sp.expand(eq)
        for eq in equations
    ]
    vars_left = list(variables)
    replacements: list[tuple[sp.Symbol, sp.Expr]] = []
    reference_assignment = dict(reference_assignment or {})
    changed = True
    while changed:
        changed = False
        for var in list(vars_left):
            for idx, eq in enumerate(list(remaining)):
                try:
                    poly = sp.Poly(eq, var)
                except Exception:
                    continue
                if poly.degree() != 1:
                    continue
                coeff = poly.coeff_monomial(var)
                const = poly.coeff_monomial(1)
                test_coeff = coeff.subs(reference_assignment)
                if could_be_zero(test_coeff):
                    continue
                replacement = sp.simplify(-const / coeff)
                replacements.append((var, replacement))
                vars_left.remove(var)
                remaining = [
                    sp.simplify(other.subs(var, replacement))
                    for j, other in enumerate(remaining)
                    if j != idx
                ]
                changed = True
                break
            if changed:
                break
    remaining = [eq for eq in remaining if not is_reliably_zero(eq)]
    return LinearEliminationResult(tuple(remaining), tuple(vars_left), tuple(replacements))


# ---------------------------------------------------------------------------
# Witness search pipeline


def _dedupe_instances(
    instances: Iterable[Mapping[sp.Symbol, object]], variables: Sequence[sp.Symbol]
) -> list[dict[sp.Symbol, sp.Expr]]:
    seen: set[tuple[str, ...]] = set()
    out: list[dict[sp.Symbol, sp.Expr]] = []
    for inst in instances:
        full = {var: sp.sympify(inst.get(var, 0)) for var in variables}
        key = tuple(sp.sstr(sp.simplify(full[var])) for var in variables)
        if key not in seen:
            seen.add(key)
            out.append(full)
    return out


def find_nonzero_polynomial_witness(
    poly: sp.Expr, variables: Sequence[sp.Symbol]
) -> dict[sp.Symbol, sp.Expr] | None:
    """Find a small integer point where a nonzero polynomial is nonzero."""

    poly = sp.expand(poly)
    used = tuple(var for var in variables if var in poly.free_symbols)
    if not used:
        return {var: sp.Integer(0) for var in variables} if not is_reliably_zero(poly) else None
    for radius in range(0, 6):
        grid = range(-radius, radius + 1)
        for values in product(grid, repeat=len(used)):
            assn = dict(zip(used, map(sp.Integer, values), strict=True))
            try:
                if not is_reliably_zero(poly.subs(assn)):
                    return {var: sp.sympify(assn.get(var, 0)) for var in variables}
            except Exception:
                continue
    return None


def _solve_univariate_rel(rel: Relational, var: sp.Symbol) -> list[sp.Expr]:
    """Return candidate sample values for a univariate real relation."""

    rel = normalize_relation(rel)
    poly = _relation_delta(rel)
    candidates: list[sp.Expr] = []
    try:
        solveset = sp.solveset(rel, var, domain=sp.S.Reals)
    except Exception:
        solveset = None
    if isinstance(solveset, sp.FiniteSet):
        candidates.extend(list(solveset))
    elif solveset is not None:
        candidates.extend([sp.Integer(0), sp.Integer(1), -sp.Integer(1)])
    try:
        roots = sp.solve(sp.Eq(poly, 0), var)
    except Exception:
        roots = []
    for root in roots:
        if root.is_real is not False:
            candidates.extend([root, root + sp.Rational(1, 3), root - sp.Rational(1, 3)])
    candidates.extend(
        [sp.Integer(0), sp.Integer(1), -sp.Integer(1), sp.Rational(1, 2), -sp.Rational(1, 2)]
    )
    out: list[sp.Expr] = []
    seen: set[str] = set()
    for candidate in candidates:
        try:
            candidate = sp.nsimplify(candidate)
        except Exception:
            candidate = sp.sympify(candidate)
        key = sp.sstr(candidate)
        if key not in seen:
            seen.add(key)
            out.append(candidate)
    return out


def _single_atom_instance(
    atom: sp.Expr,
    variables: Sequence[sp.Symbol],
    *,
    strict: bool,
) -> tuple[dict[sp.Symbol, sp.Expr] | None, str]:
    if atom is sp.true or atom is True:
        return {var: sp.Integer(0) for var in variables}, "tautology"
    if atom is sp.false or atom is False:
        return None, "contradiction"
    if not isinstance(atom, Relational) or len(variables) == 0:
        return None, "unsupported_atom"
    rel = normalize_relation(atom)
    poly = _relation_delta(rel)
    if isinstance(rel, sp.Unequality):
        witness = find_nonzero_polynomial_witness(poly, variables)
        if witness is not None and satisfies_formula(atom, witness, strict=strict):
            return witness, "nonzero_polynomial"
        return None, "unequal_no_small_witness"
    if len(variables) == 1:
        var = variables[0]
        for value in _solve_univariate_rel(rel, var):
            witness = {var: value}
            if satisfies_formula(atom, witness, strict=strict):
                return witness, "univariate_relation"
    if isinstance(rel, (sp.StrictLessThan, sp.LessThan, sp.StrictGreaterThan, sp.GreaterThan)):
        poly = sp.expand(poly)
        for var in variables:
            degree = sp.degree(poly, gen=var)
            if degree is not None and degree >= 0 and degree % 2 == 1:
                for value in (-4, -2, -1, 0, 1, 2, 4):
                    partial = {var: sp.Integer(value)}
                    rest = [v for v in variables if v != var]
                    base = find_nonzero_polynomial_witness(poly.subs(partial), rest) if rest else {}
                    candidate = {**partial, **(base or {})}
                    candidate = {v: sp.sympify(candidate.get(v, 0)) for v in variables}
                    if satisfies_formula(atom, candidate, strict=strict):
                        return candidate, "odd_degree_inequality"
        witness = find_nonzero_polynomial_witness(poly, variables)
        if witness is not None and satisfies_formula(atom, witness, strict=strict):
            return witness, "nonzero_probe"
    if isinstance(rel, sp.Equality) and len(variables) > 1:
        for value in (0, 1, -1, 2, -2):
            candidate = {var: sp.Integer(0) for var in variables}
            candidate[variables[0]] = sp.Integer(value)
            if satisfies_formula(atom, candidate, strict=strict):
                return candidate, "axis_probe"
    return None, "no_single_atom_witness"


def try_fast_witness(
    formula: sp.Expr,
    variables: Sequence[sp.Symbol],
    *,
    count: int = 1,
    strict: bool = True,
) -> FallbackInstanceResult:
    """Try cheap symbolic and algebraic real-instance heuristics."""

    variables = tuple(variables)
    attempts: list[FallbackAttempt] = []
    normalized = normalize_relations(formula)
    atoms = _atoms(normalized)
    if len(atoms) == 1:
        inst, reason = _single_atom_instance(atoms[0], variables, strict=strict)
        status = (
            "satisfied"
            if inst is not None
            else ("unsat" if reason == "contradiction" else "unknown")
        )
        attempts.append(FallbackAttempt("single_atom", status, {"reason": reason}))
        return FallbackInstanceResult(
            tuple([inst] if inst else []), status, "fast_real_heuristics", tuple(attempts)
        )
    equations = [atom for atom in atoms if isinstance(atom, sp.Equality)]
    instances: list[dict[sp.Symbol, sp.Expr]] = []
    if equations:
        try:
            sol = sp.solve(equations, variables, dict=True)
        except Exception as exc:
            sol = []
            attempts.append(FallbackAttempt("solve_equations", "unknown", {"error": str(exc)}))
        for raw in sol[: max(count * 4, 8)]:
            candidate = {var: sp.sympify(raw.get(var, 0)) for var in variables}
            if satisfies_formula(normalized, candidate, strict=strict):
                instances.append(candidate)
                if len(instances) >= count:
                    break
        if instances:
            attempts.append(
                FallbackAttempt("solve_equations", "satisfied", {"candidate_count": len(instances)})
            )
            return FallbackInstanceResult(
                tuple(_dedupe_instances(instances, variables)),
                "satisfied",
                "fast_real_heuristics",
                tuple(attempts),
            )
    origin = {var: sp.Integer(0) for var in variables}
    if satisfies_formula(normalized, origin, strict=strict):
        attempts.append(FallbackAttempt("origin_probe", "satisfied"))
        return FallbackInstanceResult(
            (origin,), "satisfied", "fast_real_heuristics", tuple(attempts)
        )
    attempts.append(FallbackAttempt("origin_probe", "unknown"))
    return FallbackInstanceResult((), "unknown", "fast_real_heuristics", tuple(attempts))


def sample_bounded_witnesses(
    formula: sp.Expr,
    variables: Sequence[sp.Symbol],
    *,
    count: int = 1,
    seed: int | None = None,
    strict: bool = True,
    sample_count: int = 96,
) -> FallbackInstanceResult:
    """Try deterministic grid probes followed by bounded random rational samples."""

    variables = tuple(variables)
    attempts: list[FallbackAttempt] = []
    candidates: list[dict[sp.Symbol, sp.Expr]] = []
    grid_values = [
        sp.Integer(0),
        sp.Integer(1),
        -sp.Integer(1),
        sp.Integer(2),
        -sp.Integer(2),
        sp.Rational(1, 2),
        -sp.Rational(1, 2),
    ]
    max_grid = min(len(grid_values) ** max(len(variables), 1), 512)
    if variables:
        for idx, values in enumerate(product(grid_values, repeat=len(variables))):
            if idx >= max_grid:
                break
            assn = dict(zip(variables, values, strict=True))
            if satisfies_formula(formula, assn, strict=strict):
                candidates.append(assn)
                if len(candidates) >= count:
                    attempts.append(FallbackAttempt("small_grid", "satisfied", {"tested": idx + 1}))
                    return FallbackInstanceResult(
                        tuple(_dedupe_instances(candidates, variables)),
                        "satisfied",
                        "bounded_sampling",
                        tuple(attempts),
                        exact=True,
                    )
    else:
        if satisfies_formula(formula, {}, strict=strict):
            return FallbackInstanceResult(
                ({},), "satisfied", "bounded_sampling", tuple(attempts), exact=True
            )
    attempts.append(FallbackAttempt("small_grid", "unknown", {"tested": max_grid}))
    rng_seed = seed if seed is not None else 1234
    for batch in range(4):
        assignments = sample_free_assignments(
            variables, sample_count=sample_count // 4, seed=rng_seed + batch
        )
        for assn in assignments:
            if satisfies_formula(formula, assn, strict=strict):
                candidates.append({var: sp.sympify(assn[var]) for var in variables})
                if len(candidates) >= count:
                    attempts.append(
                        FallbackAttempt(
                            "bounded_random", "satisfied", {"batches": batch + 1, "seed": rng_seed}
                        )
                    )
                    return FallbackInstanceResult(
                        tuple(_dedupe_instances(candidates, variables)),
                        "satisfied",
                        "bounded_sampling",
                        tuple(attempts),
                        exact=False,
                    )
    attempts.append(
        FallbackAttempt(
            "bounded_random", "unknown", {"seed": rng_seed, "sample_count": sample_count}
        )
    )
    return FallbackInstanceResult((), "unknown", "bounded_sampling", tuple(attempts), exact=False)


def sample_section_witnesses(
    formula: sp.Expr,
    variables: Sequence[sp.Symbol],
    *,
    count: int = 1,
    seed: int = 7,
    attempts: int = 7,
    strict: bool = True,
) -> FallbackInstanceResult:
    """Search for witnesses on random one-dimensional affine sections."""

    variables = tuple(variables)
    wit = find_random_section_wit(formula, variables, seed=seed, attempts=attempts)
    if wit.assignment is None:
        return FallbackInstanceResult(
            (),
            "unknown",
            "random_section_sampling",
            (
                FallbackAttempt(
                    "random_section", "unknown", {"attempts": wit.attempts, "seed": seed}
                ),
            ),
            exact=False,
        )
    candidate = {var: sp.sympify(wit.assignment[var]) for var in variables}
    if satisfies_formula(formula, candidate, strict=strict):
        return FallbackInstanceResult(
            (candidate,),
            "satisfied",
            "random_section_sampling",
            (
                FallbackAttempt(
                    "random_section", "satisfied", {"attempts": wit.attempts, "seed": seed}
                ),
            ),
            exact=False,
        )
    return FallbackInstanceResult(
        (),
        "unknown",
        "random_section_sampling",
        (FallbackAttempt("random_section", "rejected", {"attempts": wit.attempts, "seed": seed}),),
        exact=False,
    )


def solve_univar_witness(
    formula: sp.Expr,
    variables: Sequence[sp.Symbol],
    *,
    strict: bool = True,
) -> FallbackInstanceResult:
    """Try one-variable candidates from roots and decomposed factors."""

    variables = tuple(variables)
    if len(variables) != 1:
        return FallbackInstanceResult(
            (), "unknown", "univariate_decomposition", (FallbackAttempt("arity_check", "skipped"),)
        )
    var = variables[0]
    atoms = _atoms(normalize_relations(formula))
    candidates: list[sp.Expr] = []
    for atom in atoms:
        if isinstance(atom, Relational):
            candidates.extend(_solve_univariate_rel(atom, var))
            delta = _relation_delta(atom)
            try:
                for component in sp.decompose(delta, var):
                    if component != delta:
                        try:
                            candidates.extend(sp.solve(sp.Eq(component, 0), var))
                        except Exception:
                            pass
            except Exception:
                pass
    for value in candidates:
        candidate = {var: sp.nsimplify(value)}
        if satisfies_formula(formula, candidate, strict=strict):
            return FallbackInstanceResult(
                (candidate,),
                "satisfied",
                "univariate_decomposition",
                (FallbackAttempt("candidate_validation", "satisfied"),),
            )
    return FallbackInstanceResult(
        (),
        "unknown",
        "univariate_decomposition",
        (FallbackAttempt("candidate_validation", "unknown", {"candidate_count": len(candidates)}),),
    )


def find_real_witnesses(
    formula: sp.Expr,
    variables: Sequence[sp.Symbol],
    *,
    count: int = 1,
    seed: int | None = None,
    strict: bool = True,
) -> FallbackInstanceResult:
    """Layered non-CAD real instance pipeline."""

    variables = tuple(variables)
    collected: list[dict[sp.Symbol, sp.Expr]] = []
    attempts: list[FallbackAttempt] = []
    normalized = normalize_relations(formula)
    methods = (
        lambda: try_fast_witness(normalized, variables, count=count, strict=strict),
        lambda: sample_bounded_witnesses(
            normalized, variables, count=count, seed=seed, strict=strict
        ),
        lambda: sample_section_witnesses(
            normalized, variables, count=count, seed=seed if seed is not None else 7, strict=strict
        ),
        lambda: solve_univar_witness(normalized, variables, strict=strict),
    )
    for method in methods:
        result = method()
        attempts.extend(result.attempts)
        if result.instances:
            collected.extend(result.instances)
            collected = _dedupe_instances(collected, variables)
            return FallbackInstanceResult(
                tuple(collected[:count]),
                "satisfied",
                result.method,
                tuple(attempts),
                exact=result.exact,
            )
        if result.status == "unsat":
            return FallbackInstanceResult(
                (), "unsat", result.method, tuple(attempts), exact=result.exact
            )
    return FallbackInstanceResult(
        tuple(collected),
        "satisfied" if collected else "unknown",
        "real_instance_fallback_pipeline",
        tuple(attempts),
        exact=False,
    )


# ---------------------------------------------------------------------------
# Staged method orchestration


class _Timeout(Exception):
    pass


def _timeout_handler(signum, frame):  # type: ignore[no-untyped-def]
    raise _Timeout()


def try_methods(
    methods: Sequence[Callable[..., object]],
    args: Sequence[object] = (),
    *,
    initial_seconds: int = 1,
    growth_factor: int = 10,
    failure_value: object = None,
) -> MethodSearchResult:
    """Run methods with increasing time budgets.

    A method result that does not contain ``failure_value`` is considered a full
    success. Partial results containing other data plus ``failure_value`` are
    recorded and returned if no method fully succeeds.
    """

    remaining = list(methods)
    partial: list[object] = []
    budget = initial_seconds
    while remaining:
        timed_out: list[Callable[..., object]] = []
        for method in remaining:
            old_handler = signal.signal(signal.SIGALRM, _timeout_handler)
            signal.alarm(max(1, int(budget)))
            try:
                result = method(*args)
            except _Timeout:
                timed_out.append(method)
                continue
            except Exception as exc:
                result = failure_value
                partial.append(exc)
            finally:
                signal.alarm(0)
                signal.signal(signal.SIGALRM, old_handler)
            if result is failure_value:
                continue
            if isinstance(result, (list, tuple, set)) and failure_value in result:
                partial.append(result)
                continue
            return MethodSearchResult(
                result, "satisfied", tuple(partial), getattr(method, "__name__", repr(method))
            )
        remaining = timed_out
        budget = budget * growth_factor if len(remaining) > 1 else 10**9
    if partial:
        return MethodSearchResult(tuple(partial), "partial", tuple(partial))
    return MethodSearchResult(failure_value, "unknown")


__all__ = [
    "CoordinateBounds",
    "FallbackAttempt",
    "FallbackInstanceResult",
    "LinearEliminationResult",
    "MethodSearchResult",
    "algebraic_variables",
    "as_rational_if_exact",
    "satisfies_formula",
    "could_be_zero",
    "rational_bound",
    "eliminate_linear_equations",
    "coordinate_bounds",
    "expand_with_subs",
    "fast_factor_list",
    "find_nonzero_polynomial_witness",
    "find_real_witnesses",
    "solve_univar_witness",
    "formula_to_rule_sets",
    "is_algebraic_condition",
    "is_rational_number",
    "is_reliably_zero",
    "is_valid_numeric_value",
    "normalize_relation",
    "normalize_relations",
    "is_bounded_solution_set",
    "try_fast_witness",
    "try_methods",
    "sample_bounded_witnesses",
    "sample_section_witnesses",
    "sort_conditions",
    "to_dnf_formula",
    "vector_relations",
]
