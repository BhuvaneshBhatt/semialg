"""Expected failures from optional or undecidable symbolic computations.

Fallback code should catch the smallest operation-specific family that can
legitimately be raised by the SymPy operation being attempted.  The aggregate
``SYMBOLIC_QUERY_ERRORS`` remains for compatibility at broad orchestration
boundaries; proof code should prefer the narrower groups below.
"""

from sympy.core.function import PoleError
from sympy.core.sympify import SympifyError
from sympy.polys.polyerrors import PolynomialError

SYMPIFY_ERRORS = (SympifyError, TypeError)
POLYNOMIAL_ERRORS = (PolynomialError, ValueError)
SIMPLIFICATION_ERRORS = (TypeError, ValueError, NotImplementedError, RecursionError)
SET_OPERATION_ERRORS = (TypeError, ValueError, NotImplementedError)
LIMIT_ERRORS = (PoleError, ValueError, NotImplementedError)
BACKEND_QUERY_ERRORS = (ValueError, NotImplementedError, ArithmeticError)

# Compatibility umbrella for high-level fallback boundaries.  Deliberately
# excludes SyntaxError: malformed source is a programming/input bug, not an
# expected symbolic-decision failure.
SYMBOLIC_QUERY_ERRORS = (
    *SIMPLIFICATION_ERRORS,
    ArithmeticError,
    SympifyError,
    PolynomialError,
    PoleError,
)
