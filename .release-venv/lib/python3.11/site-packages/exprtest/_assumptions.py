"""Utilities for applying and validating caller-supplied assumptions."""

from __future__ import annotations

import sympy as sp
from sympy.assumptions import assuming

from ._errors import RECOVERABLE_SYMPY_ERRORS


def normalize_assumptions(assumptions):
    """Return a stable SymPy representation of an assumptions argument."""
    if assumptions is None or assumptions is True or assumptions is sp.true:
        return sp.true
    if assumptions is False or assumptions is sp.false:
        return sp.false
    if isinstance(assumptions, (tuple, list, set, frozenset)):
        items = tuple(sp.sympify(item) for item in assumptions)
        return sp.And(*items) if items else sp.true
    return sp.sympify(assumptions)


def assumption_facts(assumptions) -> tuple:
    """Flatten a conjunction into facts suitable for ``sympy.assuming``."""
    assumptions = normalize_assumptions(assumptions)
    if assumptions is sp.true:
        return ()
    if isinstance(assumptions, sp.And):
        return tuple(assumptions.args)
    return (assumptions,)


def equality_substitutions(assumptions) -> dict:
    """Extract safe direct substitutions from conjunctive equality facts."""
    mapping = {}
    for fact in assumption_facts(assumptions):
        if not isinstance(fact, sp.Equality):
            continue
        lhs, rhs = fact.lhs, fact.rhs
        if lhs.is_Symbol and lhs not in rhs.free_symbols:
            mapping[lhs] = rhs
        elif rhs.is_Symbol and rhs not in lhs.free_symbols:
            mapping[rhs] = lhs
    if not mapping:
        return mapping
    # Resolve simple chains such as x=y, y=2 without invoking a solver.
    for _ in range(len(mapping)):
        updated = {key: value.xreplace(mapping) for key, value in mapping.items()}
        if updated == mapping:
            break
        mapping = updated
    return mapping


def refine_with_assumptions(expr: sp.Expr, assumptions) -> sp.Expr:
    """Apply assumption-aware, value-preserving SymPy refinement."""
    assumptions = normalize_assumptions(assumptions)
    if assumptions is sp.true:
        return expr
    if assumptions is sp.false:
        return expr
    substitutions = equality_substitutions(assumptions)
    if substitutions:
        expr = expr.xreplace(substitutions).subs(substitutions)
    try:
        refined = sp.refine(expr, assumptions)
    except RECOVERABLE_SYMPY_ERRORS:
        refined = expr
    return refined


def assumptions_hold(assumptions, substitutions: dict) -> bool | None:
    """Evaluate assumptions at an exact sample point.

    Returns ``True`` or ``False`` when SymPy can decide the instantiated
    assumptions and ``None`` when it cannot.  Unknown points are rejected by
    witness sampling rather than treated as valid.
    """
    assumptions = normalize_assumptions(assumptions)
    if assumptions is sp.true:
        return True
    if assumptions is sp.false:
        return False
    try:
        instantiated = assumptions.xreplace(substitutions).subs(substitutions)
    except RECOVERABLE_SYMPY_ERRORS:
        return None
    if instantiated is sp.true or instantiated is True:
        return True
    if instantiated is sp.false or instantiated is False:
        return False
    try:
        asked = sp.ask(instantiated)
        if asked is not None:
            return bool(asked)
    except RECOVERABLE_SYMPY_ERRORS:
        pass
    try:
        return bool(instantiated)
    except (TypeError, ValueError):
        return None


def has_nontrivial_assumptions(assumptions) -> bool:
    return normalize_assumptions(assumptions) is not sp.true


def ask_property(predicate, term: sp.Expr, assumptions=True) -> bool | None:
    """Ask one cheap SymPy predicate under caller assumptions.

    This is the shared bounded assumptions gateway used by identity,
    definedness, and nonzero reasoning.  It deliberately does not invoke a
    solver or general simplification.
    """
    assumptions = normalize_assumptions(assumptions)
    try:
        facts = assumption_facts(assumptions)
        query = predicate(sp.sympify(term))
        if facts:
            with assuming(*facts):
                value = sp.ask(query)
        else:
            value = sp.ask(query)
    except RECOVERABLE_SYMPY_ERRORS:
        return None
    return None if value is None else bool(value)


def is_positive(term: sp.Expr, assumptions=True) -> bool | None:
    """Return a cheap proof/refutation of positivity under assumptions."""
    value = sp.sympify(term).is_positive
    if value is not None:
        return bool(value)
    return ask_property(sp.Q.positive, term, assumptions)


def is_nonzero(term: sp.Expr, assumptions=True) -> bool | None:
    """Return a cheap proof/refutation of nonzeroness under assumptions."""
    value = sp.sympify(term).is_zero
    if value is not None:
        return not bool(value)
    return ask_property(sp.Q.nonzero, term, assumptions)


def is_finite(term: sp.Expr, assumptions=True) -> bool | None:
    """Return a cheap proof/refutation of finiteness under assumptions."""
    value = sp.sympify(term).is_finite
    if value is not None:
        return bool(value)
    return ask_property(sp.Q.finite, term, assumptions)


def is_negative(term: sp.Expr, assumptions=True) -> bool | None:
    """Return a cheap proof/refutation of negativity under assumptions."""
    value = sp.sympify(term).is_negative
    if value is not None:
        return bool(value)
    return ask_property(sp.Q.negative, term, assumptions)


def is_nonnegative(term: sp.Expr, assumptions=True) -> bool | None:
    """Return a cheap proof/refutation of nonnegativity under assumptions."""
    value = sp.sympify(term).is_nonnegative
    if value is not None:
        return bool(value)
    return ask_property(sp.Q.nonnegative, term, assumptions)


def is_nonpositive(term: sp.Expr, assumptions=True) -> bool | None:
    """Return a cheap proof/refutation of nonpositivity under assumptions."""
    value = sp.sympify(term).is_nonpositive
    if value is not None:
        return bool(value)
    return ask_property(sp.Q.nonpositive, term, assumptions)


def is_real_value(term: sp.Expr, assumptions=True) -> bool | None:
    """Return a cheap proof/refutation that a value is real."""
    value = sp.sympify(term).is_real
    if value is not None:
        return bool(value)
    return ask_property(sp.Q.real, term, assumptions)


def is_integer_value(term: sp.Expr, assumptions=True) -> bool | None:
    """Return a cheap proof/refutation that a value is integral."""
    value = sp.sympify(term).is_integer
    if value is not None:
        return bool(value)
    return ask_property(sp.Q.integer, term, assumptions)
