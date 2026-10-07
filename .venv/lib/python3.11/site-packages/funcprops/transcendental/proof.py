from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

import sympy as sp
from sympy import S

from ..exactness import certify_zero
from ..exceptions import SYMBOLIC_QUERY_ERRORS


class PerformanceGoal(str, Enum):
    SPEED = "speed"
    QUALITY = "quality"


@dataclass(frozen=True)
class ProofPolicy:
    goal: PerformanceGoal = PerformanceGoal.QUALITY
    allow_numerical_evidence: bool = True
    require_certificate: bool = True
    max_symbolic_recursion: int = 3


def extended_real_sign(value):
    value = sp.sympify(value)
    if value is sp.oo:
        return 1
    if value is -sp.oo:
        return -1
    if certify_zero(value) is True:
        return 0
    if value.is_positive is True:
        return 1
    if value.is_negative is True:
        return -1
    # Do not certify a sign from a fixed-precision approximation.  This helper
    # is used by proof infrastructure, so expressions too close to zero (or
    # whose exact comparison is simply unknown) must remain unknown.
    try:
        q = sp.ask(sp.Q.positive(value))
        if q is True:
            return 1
        q = sp.ask(sp.Q.negative(value))
        if q is True:
            return -1
        q = sp.ask(sp.Q.zero(value))
        if q is True:
            return 0
    except SYMBOLIC_QUERY_ERRORS:
        pass
    return None


def find_parameters(expressions, variables):
    """Return free symbolic parameters independent of ``variables``."""
    variables = set(map(sp.sympify, variables))
    if not isinstance(expressions, (tuple, list, set, frozenset)):
        expressions = (expressions,)
    free = set()
    for e in expressions:
        if isinstance(e, sp.Basic):
            free.update(e.free_symbols)
    return tuple(sorted(free - variables, key=sp.default_sort_key))


def condition_to_intervals(condition, var, *, integers=False):
    """Parse a solved univariate condition into ordered intervals/points.

    This is the SymPy analogue of the useful CADToIntervals helper in the
    supplied FunctionProperties source.  ``None`` means the condition is not
    in a supported solved form.
    """
    condition = sp.sympify(condition)
    var = sp.sympify(var)
    try:
        s = (
            condition.as_set()
            if condition not in (S.true, S.false)
            else (S.Reals if condition is S.true else S.EmptySet)
        )
    except (AttributeError, *SYMBOLIC_QUERY_ERRORS):
        return None
    if integers:
        s = sp.Intersection(s, S.Integers)
    out = []

    def visit(t):
        if t is S.EmptySet:
            return True
        if isinstance(t, sp.Interval):
            out.append(t)
            return True
        if isinstance(t, sp.FiniteSet):
            out.extend(sp.FiniteSet(x) for x in sorted(t, key=sp.default_sort_key))
            return True
        if isinstance(t, sp.Union):
            return all(visit(a) for a in t.args)
        if isinstance(t, sp.Range):
            out.append(t)
            return True
        return False

    return tuple(out) if visit(s) else None
