"""Expected exceptions at bounded SymPy fallback boundaries."""

from sympy.core.sympify import SympifyError
from sympy.polys.polyerrors import (
    CoercionFailed,
    DomainError,
    ExactQuotientFailed,
    GeneratorsNeeded,
    NotAlgebraic,
    PolynomialError,
)

EXACT_METHOD_ERRORS = (
    ArithmeticError,
    TypeError,
    ValueError,
    NotImplementedError,
    CoercionFailed,
    DomainError,
    GeneratorsNeeded,
    NotAlgebraic,
    PolynomialError,
)

RECOVERABLE_SYMPY_ERRORS = (
    TypeError,
    ValueError,
    AttributeError,
    NotImplementedError,
    SympifyError,
)

RECOVERABLE_POLY_ERRORS = RECOVERABLE_SYMPY_ERRORS + (
    PolynomialError,
    CoercionFailed,
    DomainError,
    ExactQuotientFailed,
    GeneratorsNeeded,
)
