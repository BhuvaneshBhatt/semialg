"""Cheap exact certificates used before algebraic decision procedures."""

from __future__ import annotations

import sympy as sp
from sympy.core.relational import Relational
from sympy.logic.boolalg import And, Boolean


def propositional_truth_value(formula: sp.Expr | Boolean) -> bool | None:
    """Return a value when Boolean structure alone decides ``formula``.

    Relational atoms are represented by independent propositional symbols.
    Exact complementary relations share one symbol with opposite polarity, so
    normalization from ``Not(x > 0)`` to ``x <= 0`` preserves the certificate.
    """
    # Conjunctive algebraic systems are handled more cheaply by the linear,
    # univariate, and finite-system specialists.  Propositional minimization is
    # useful when Boolean branching can itself certify the result.
    if isinstance(formula, And) and all(isinstance(arg, Relational) for arg in formula.args):
        return None

    atoms = tuple(sorted(formula.atoms(Relational), key=sp.sstr))
    if not atoms:
        if formula is sp.true:
            return True
        if formula is sp.false:
            return False
        return None

    replacements: dict[sp.Expr, sp.Expr] = {}
    for atom in atoms:
        if atom in replacements:
            continue
        symbol = sp.Dummy(boolean=True)
        replacements[atom] = symbol
        complement = sp.Not(atom)
        if isinstance(complement, Relational):
            replacements[complement] = sp.Not(symbol)

    skeleton = formula.xreplace(replacements)
    reduced = sp.simplify_logic(skeleton, deep=False)
    if reduced is sp.true:
        return True
    if reduced is sp.false:
        return False
    return None
