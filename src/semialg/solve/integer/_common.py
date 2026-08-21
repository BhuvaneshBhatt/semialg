from __future__ import annotations

import sympy as sp

RECOVERABLE_ERRORS = (
    ArithmeticError,
    TypeError,
    ValueError,
    NotImplementedError,
    RuntimeError,
    sp.PolynomialError,
)


def expr_complexity(expr: sp.Expr) -> int:
    """Return a cheap deterministic expression-complexity score."""

    try:
        return int(sp.count_ops(expr, visual=False))
    except RECOVERABLE_ERRORS:
        return len(sp.srepr(expr))


__all__ = ["RECOVERABLE_ERRORS", "expr_complexity"]
