"""Shared normalization and bounded implication reasoning for symbolic assumptions."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from functools import lru_cache

import sympy as sp
from sympy import S

from .assumption_normalization import normalize_assumptions
from .exactness import certify_equal


class PropertyKnowledge(str, Enum):
    EXACT = "exact"
    SUFFICIENT = "sufficient"
    NECESSARY = "necessary"
    PARTIAL = "partial"


@dataclass(frozen=True)
class PropertyProvenance:
    source: str
    reference: str | None = None
    note: str | None = None


@dataclass(frozen=True)
class PropertyDecision:
    predicate: sp.Expr
    verdict: bool | None
    assumptions: sp.Expr = S.true
    knowledge: PropertyKnowledge = PropertyKnowledge.EXACT
    provenance: tuple[PropertyProvenance, ...] = ()
    reasons: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        object.__setattr__(self, "predicate", sp.sympify(self.predicate))
        object.__setattr__(self, "assumptions", normalize_assumptions(self.assumptions))

    @property
    def certified(self) -> bool:
        return (
            self.verdict is not None and self.knowledge is not PropertyKnowledge.PARTIAL
        )


class PropertyEnforcementError(ValueError):
    def __init__(self, decision: PropertyDecision, *, operation: str | None = None):
        self.decision = decision
        self.operation = operation
        state = "false" if decision.verdict is False else "unresolved"
        prefix = f"{operation}: " if operation else ""
        super().__init__(f"{prefix}required property is {state}: {decision.predicate}")


def _affine_atom(atom: sp.Expr) -> sp.Expr:
    if not isinstance(atom, (sp.Equality, sp.Unequality)):
        return atom
    expr = sp.expand(atom.lhs - atom.rhs)
    symbols = tuple(sorted(expr.free_symbols, key=sp.default_sort_key))
    if len(symbols) != 1:
        return atom
    symbol = symbols[0]
    try:
        poly = sp.Poly(expr, symbol)
    except (sp.PolynomialError, TypeError, ValueError):
        return atom
    if poly.degree() != 1:
        return atom
    a, b = poly.all_coeffs()
    if a == 0:
        return atom
    value = sp.cancel(-b / a)
    head = sp.Eq if isinstance(atom, sp.Equality) else sp.Ne
    return head(symbol, value, evaluate=False)


def _normalize_logic(expr: sp.Expr | bool) -> sp.Expr:
    value = normalize_assumptions(expr)
    if isinstance(value, (sp.And, sp.Or)):
        return value.func(
            *(_normalize_logic(arg) for arg in value.args), evaluate=False
        )
    if isinstance(value, sp.Not):
        arg = _normalize_logic(value.args[0])
        if isinstance(arg, sp.Equality):
            return sp.Ne(arg.lhs, arg.rhs, evaluate=False)
        if isinstance(arg, sp.Unequality):
            return sp.Eq(arg.lhs, arg.rhs, evaluate=False)
        return sp.Not(arg, evaluate=False)
    return _affine_atom(value)


def _atom_key(expr: sp.Expr) -> str:
    return sp.srepr(_normalize_logic(expr))


def _structural_entailment(condition: sp.Expr, assumptions: sp.Expr) -> bool | None:
    if condition is S.true or assumptions is S.false:
        return True
    if condition is S.false:
        return False
    condition = _normalize_logic(condition)
    assumptions = _normalize_logic(assumptions)
    if _atom_key(condition) == _atom_key(assumptions):
        return True
    clauses = sp.And.make_args(assumptions)
    keys = {_atom_key(clause) for clause in clauses}
    if _atom_key(condition) in keys:
        return True
    # Cheap ordered-real closure: an explicit comparison with zero proves
    # realness; a strict comparison also proves nonzeroness.
    if isinstance(condition, sp.Contains) and condition.args[1] is S.Reals:
        target = condition.args[0]
        for clause in clauses:
            if isinstance(
                clause,
                (sp.StrictGreaterThan, sp.StrictLessThan, sp.GreaterThan, sp.LessThan),
            ):
                if (clause.lhs == target and clause.rhs == 0) or (
                    clause.rhs == target and clause.lhs == 0
                ):
                    return True
    if isinstance(condition, sp.Unequality):
        target = (
            condition.lhs
            if condition.rhs == 0
            else condition.rhs
            if condition.lhs == 0
            else None
        )
        if target is not None:
            for clause in clauses:
                if isinstance(clause, (sp.StrictGreaterThan, sp.StrictLessThan)) and (
                    (clause.lhs == target and clause.rhs == 0)
                    or (clause.rhs == target and clause.lhs == 0)
                ):
                    return True
    neg = _normalize_logic(sp.Not(condition, evaluate=False))
    if _atom_key(neg) in keys:
        return False
    if isinstance(condition, sp.And):
        vals = tuple(_structural_entailment(arg, assumptions) for arg in condition.args)
        if all(v is True for v in vals):
            return True
        if any(v is False for v in vals):
            return False
    if isinstance(condition, sp.Or):
        vals = tuple(_structural_entailment(arg, assumptions) for arg in condition.args)
        if any(v is True for v in vals):
            return True
        if all(v is False for v in vals):
            return False
    return None


def _affine_equalities_consistent(expr: sp.Expr) -> bool | None:
    clauses = sp.And.make_args(_normalize_logic(expr))
    eqs = [c for c in clauses if isinstance(c, sp.Equality)]
    neqs = [c for c in clauses if isinstance(c, sp.Unequality)]
    if not eqs and not neqs:
        return None

    # First collapse direct symbol equalities.  This retains exact affine
    # consequences such as x = a, x = b => a = b without asking the zero
    # oracle to decide an unconstrained symbolic difference.
    parent: dict[sp.Symbol, sp.Symbol] = {}

    def find(symbol: sp.Symbol) -> sp.Symbol:
        parent.setdefault(symbol, symbol)
        while parent[symbol] != symbol:
            parent[symbol] = parent[parent[symbol]]
            symbol = parent[symbol]
        return symbol

    def union(left: sp.Symbol, right: sp.Symbol) -> None:
        a, b = find(left), find(right)
        if a != b:
            # Use a stable representative so normalization is determinstic.
            if sp.default_sort_key(a) <= sp.default_sort_key(b):
                parent[b] = a
            else:
                parent[a] = b

    for eq in eqs:
        if isinstance(eq.lhs, sp.Symbol) and isinstance(eq.rhs, sp.Symbol):
            union(eq.lhs, eq.rhs)

    def collapse(value):
        return sp.sympify(value).xreplace(
            {symbol: find(symbol) for symbol in tuple(parent)}
        )

    # After collapsing equality classes, a residual equality/inequality may
    # become a constant contradiction (for example a = b + 1 when a = b).
    # Only the certified zero boundary may turn that residual into a proof.
    for atom in (*eqs, *neqs):
        verdict = certify_equal(collapse(atom.lhs), collapse(atom.rhs))
        if isinstance(atom, sp.Equality) and verdict is False:
            return False
        if isinstance(atom, sp.Unequality) and verdict is True:
            return False

    values: dict[sp.Symbol, sp.Expr] = {}
    for eq in eqs:
        lhs, rhs = collapse(eq.lhs), collapse(eq.rhs)
        q = _affine_atom(sp.Eq(lhs, rhs, evaluate=False))
        if not isinstance(q, sp.Equality) or not isinstance(q.lhs, sp.Symbol):
            continue
        prior = values.get(q.lhs)
        if prior is not None and certify_equal(prior, q.rhs) is False:
            return False
        values[q.lhs] = q.rhs
    for neq in neqs:
        lhs, rhs = collapse(neq.lhs), collapse(neq.rhs)
        q = _affine_atom(sp.Ne(lhs, rhs, evaluate=False))
        if (
            isinstance(q, sp.Unequality)
            and isinstance(q.lhs, sp.Symbol)
            and q.lhs in values
            and certify_equal(values[q.lhs], q.rhs) is True
        ):
            return False
    return True if values or parent else None


def _ask_refine(condition: sp.Expr, assumptions: sp.Expr) -> bool | None:
    try:
        refined = sp.refine(condition, assumptions)
    except (TypeError, ValueError, NotImplementedError):
        refined = condition
    if refined is S.true:
        return True
    if refined is S.false:
        return False
    try:
        asked = sp.ask(condition, assumptions)
    except (TypeError, ValueError, NotImplementedError):
        asked = None
    return bool(asked) if asked is not None else None


@lru_cache(maxsize=4096)
def _entails_cached(condition: sp.Expr, assumptions: sp.Expr) -> bool | None:
    structural = _structural_entailment(condition, assumptions)
    if structural is not None:
        return structural

    # Small Boolean partitions are common in parameter stratification.  This
    # is propositional simplification only; it does not invoke general
    # arithmetic SAT/QE.  Eq/Ne pairs are understood by SymPy as logical
    # complements, which proves exhaustive zero/nonzero partitions cheaply.
    implication = sp.Implies(assumptions, condition, evaluate=False)
    if sp.count_ops(implication) <= 40:
        try:
            logical = sp.simplify_logic(implication, force=True)
        except (NotImplementedError, TypeError, ValueError):
            logical = implication
        if logical is S.true:
            return True
        if logical is S.false:
            return False

    positive_consistency = _affine_equalities_consistent(
        sp.And(assumptions, sp.Not(condition, evaluate=False), evaluate=False)
    )
    if positive_consistency is False:
        return True
    negative_consistency = _affine_equalities_consistent(
        sp.And(assumptions, condition, evaluate=False)
    )
    if negative_consistency is False:
        return False
    return _ask_refine(condition, assumptions)


def entails(condition: sp.Expr, assumptions: sp.Expr | bool = True) -> bool | None:
    """Return whether normalized assumptions entail ``condition``.

    The implementation is deliberately bounded: structural and affine Eq/Ne
    reasoning precede SymPy's assumptions engine; unresolved predicates remain
    unresolved.
    """
    return _entails_cached(_normalize_logic(condition), _normalize_logic(assumptions))


def decide(
    predicate,
    assumptions: sp.Expr | bool = True,
    *,
    knowledge=PropertyKnowledge.EXACT,
    provenance: tuple[PropertyProvenance, ...] = (),
    reasons: tuple[str, ...] = (),
) -> PropertyDecision:
    assumptions = _normalize_logic(assumptions)
    return PropertyDecision(
        predicate,
        entails(predicate, assumptions),
        assumptions,
        knowledge,
        provenance,
        reasons,
    )


def require_decision(
    decision: PropertyDecision,
    *,
    operation: str | None = None,
    allow_unknown: bool = False,
) -> PropertyDecision:
    if decision.verdict is False or (decision.verdict is None and not allow_unknown):
        raise PropertyEnforcementError(decision, operation=operation)
    return decision


def clear_entailment_cache() -> None:
    _entails_cached.cache_clear()


def entailment_cache_info():
    return _entails_cached.cache_info()
