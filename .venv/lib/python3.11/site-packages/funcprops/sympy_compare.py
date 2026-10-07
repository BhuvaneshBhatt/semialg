from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

from sympy import S
from sympy.calculus.util import continuous_domain

from .conditions import parse_catalogue_expression
from .exceptions import SYMBOLIC_QUERY_ERRORS


class ComparisonStatus(str, Enum):
    AGREE = "agree"
    DISAGREE = "disagree"
    SYMPY_UNKNOWN = "sympy_unknown"
    NOT_COMPARABLE = "not_comparable"
    ERROR = "error"


@dataclass(frozen=True)
class Comparison:
    status: ComparisonStatus
    registry_value: object
    sympy_value: object = None
    detail: str = ""


def compare_real_domain(expr, var, registry_condition):
    """Compare a catalogue real-domain rule with SymPy when both are expressible.

    SymPy's continuous_domain returns a Set while catalogue rules are Boolean
    conditions, so exact equality is only asserted for readily convertible cases.
    """
    try:
        s = continuous_domain(expr, var, S.Reals)
    except NotImplementedError as e:
        return Comparison(
            ComparisonStatus.SYMPY_UNKNOWN, registry_condition, None, str(e)
        )
    except SYMBOLIC_QUERY_ERRORS as e:
        return Comparison(ComparisonStatus.ERROR, registry_condition, None, str(e))
    c = (
        parse_catalogue_expression(registry_condition)
        if isinstance(registry_condition, str)
        else registry_condition
    )
    if c is S.true:
        return Comparison(
            ComparisonStatus.AGREE if s == S.Reals else ComparisonStatus.DISAGREE, c, s
        )
    if c is S.false:
        return Comparison(
            ComparisonStatus.AGREE if s == S.EmptySet else ComparisonStatus.DISAGREE,
            c,
            s,
        )
    return Comparison(
        ComparisonStatus.NOT_COMPARABLE,
        c,
        s,
        "Boolean-condition to Set equivalence requires quantifier/set normalization",
    )
