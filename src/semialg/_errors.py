from __future__ import annotations

import sympy as sp
from sympy.polys.polyerrors import CoercionFailed

EXACT_OPERATION_ERRORS = (
    ArithmeticError,
    TypeError,
    ValueError,
    NotImplementedError,
    sp.PolynomialError,
    CoercionFailed,
)

POLY_OPERATION_ERRORS = (
    ArithmeticError,
    TypeError,
    ValueError,
    NotImplementedError,
    sp.PolynomialError,
    CoercionFailed,
)

__all__ = ["EXACT_OPERATION_ERRORS", "POLY_OPERATION_ERRORS"]
