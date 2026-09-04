"""Conservative adapters around SymPy's one-dimensional inequality reducer."""

from __future__ import annotations

import sympy as sp
from sympy.core.relational import Relational

from ._errors import EXACT_OPERATION_ERRORS


def reduce_conjunctive_inequalities(expr: sp.Expr, variable: sp.Symbol) -> sp.Expr | None:
    """Reduce a conjunction of relational atoms, or return ``None``.

    ``sympy.reduce_inequalities`` is not a general Boolean formula reducer and
    can return ``False`` for unsupported disjunctive input.  Callers must not
    interpret that failure mode as a proof of inconsistency.
    """

    if expr in (True, sp.true, False, sp.false):
        return sp.sympify(expr)
    atoms = expr.args if isinstance(expr, sp.And) else (expr,)
    if not atoms or not all(isinstance(atom, Relational) for atom in atoms):
        return None
    try:
        return sp.reduce_inequalities(list(atoms), variable)
    except EXACT_OPERATION_ERRORS:
        return None


__all__ = ["reduce_conjunctive_inequalities"]
