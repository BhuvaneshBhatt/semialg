"""Semantic predicates and parsing for the bundled property catalogue.

Catalogue conditions are stored as SymPy constructor expressions.  The custom
Boolean predicates below preserve mathematical assumptions that do not map
cleanly to ordinary ``Symbol`` assumptions until the proof engine has enough
context to reduce them.
"""

from __future__ import annotations

from functools import lru_cache

from sympy import S, sympify
from sympy.core.relational import Relational
from sympy.logic.boolalg import Boolean, BooleanFunction


class AIntegerQ(BooleanFunction):
    pass


class ARationalQ(BooleanFunction):
    pass


class ARealQ(BooleanFunction):
    pass


class ARealIfDefinedQ(BooleanFunction):
    pass


class AAlgebraicQ(BooleanFunction):
    pass


class ANonIntegerQ(BooleanFunction):
    pass


class APositive(BooleanFunction):
    pass


class ANegative(BooleanFunction):
    pass


class ANonPositive(BooleanFunction):
    pass


class ANonNegative(BooleanFunction):
    pass


class AEvenQ(BooleanFunction):
    pass


class AOddQ(BooleanFunction):
    pass


class AUnequalQ(BooleanFunction):
    pass


class AEqualQ(BooleanFunction):
    pass


class AAssumedQ(BooleanFunction):
    pass


class AAssumedIneqQ(BooleanFunction):
    pass


class NonTerminatingQ(BooleanFunction):
    pass


_PREDICATE_TYPES = (
    AIntegerQ,
    ARationalQ,
    ARealQ,
    ARealIfDefinedQ,
    AAlgebraicQ,
    ANonIntegerQ,
    APositive,
    ANegative,
    ANonPositive,
    ANonNegative,
    AEvenQ,
    AOddQ,
    AUnequalQ,
    AEqualQ,
    AAssumedQ,
    AAssumedIneqQ,
    NonTerminatingQ,
)
_PREDICATES = {cls.__name__: cls for cls in _PREDICATE_TYPES}


@lru_cache(maxsize=8192)
def parse_catalogue_expression(text: str):
    """Parse one bundled catalogue value or condition into a SymPy object.

    The catalogue uses ``srepr``-compatible constructor syntax, so parsing is
    deterministic and does not require a second symbolic-language grammar.
    Custom semantic predicates remain explicit ``BooleanFunction`` nodes.
    """
    source = text.strip()
    if not source:
        return S.true
    return sympify(source, locals=_PREDICATES, evaluate=False)


def is_boolean_expr(expr) -> bool:
    """Return whether ``expr`` is a Boolean or relational SymPy expression."""
    return isinstance(expr, (Boolean, Relational)) or expr in (S.true, S.false)
