"""Exact normalization and comparison for closed numeric constants."""

from __future__ import annotations

from dataclasses import dataclass

import sympy as sp
from sympy.polys.numberfields import to_number_field

from . import _config as cfg
from ._algebraic import (
    classify_algebraic_expression,
    normalize_radicals,
    simplify_root_sets,
)
from ._cost import stage_allowed, within_budget
from ._errors import EXACT_METHOD_ERRORS
from ._memo import exact_root, minpoly_for


@dataclass(frozen=True)
class AlgebraicOrder:
    """Exact comparison result for two real algebraic constants."""

    order: int
    difference: sp.Expr
    method: str = "exact-algebraic-order"


def _field_zero(value) -> bool:
    if isinstance(value, sp.AlgebraicNumber):
        coeffs = tuple(value.coeffs())
        if coeffs and all(part == 0 for part in coeffs):
            return True
        return value.as_expr() == 0
    return value == 0


def algebraic_sign(term: sp.Expr) -> int | None:
    """Return -1, 0, or 1 for a real algebraic constant, else ``None``."""
    if not stage_allowed(sp.sympify(term), "minpoly"):
        return None

    term = sp.sympify(term)
    if term.has(sp.sin, sp.cos, sp.tan):
        from ._cyclotomic import cyclotomic_sign

        sign = cyclotomic_sign(term)
        if sign is not None:
            return sign
    info = classify_algebraic_expression(term)
    if not info.is_algebraic or term.is_real is not True:
        return None
    try:
        value = to_number_field(term)
        if _field_zero(value):
            return 0
        exact = exact_root(term)
        if exact.is_positive is True:
            return 1
        if exact.is_negative is True:
            return -1
        pos = sp.ask(sp.Q.positive(exact))
        if pos is True:
            return 1
        neg = sp.ask(sp.Q.negative(exact))
        if neg is True:
            return -1
    except EXACT_METHOD_ERRORS:
        return None
    return None


def compare_algebraic(left: sp.Expr, right: sp.Expr) -> AlgebraicOrder | None:
    """Compare two real algebraic constants exactly when possible."""
    if not within_budget(sp.Add(sp.sympify(left), -sp.sympify(right), evaluate=False)):
        return None

    diff = sp.Add(sp.sympify(left), -sp.sympify(right), evaluate=False)
    sign = algebraic_sign(diff)
    if sign is None:
        return None
    return AlgebraicOrder(sign, sp.cancel(diff))


def _rational_pi(arg: sp.Expr) -> sp.Rational | None:
    coeff = arg.coeff(sp.pi)
    if arg == coeff * sp.pi and coeff.is_Rational:
        return sp.Rational(coeff)
    return None


def reduce_exact_trig(term: sp.Expr, degree_limit: int | None = None) -> sp.Expr:
    """Convert rational-angle trig constants to exact algebraic values.

    A cyclotomic quotient reduction is attempted first so identities among
    several trig constants can collapse before generic number-field work.
    """
    if not stage_allowed(sp.sympify(term), "minpoly"):
        return sp.sympify(term)

    term = sp.sympify(term)
    budget = cfg.TRIG_FIELD_MAX_DEGREE if degree_limit is None else int(degree_limit)

    from ._cyclotomic import cyclotomic_form

    form = cyclotomic_form(term)
    if form is not None and form.is_zero:
        return sp.Integer(0)

    # For nonzero individual constants SymPy's direct algebraic conversion is
    # substantially cheaper than rebuilding the element from an exponential
    # root-of-unity spelling. Cyclotomic arithmetic is still used first for
    # exact identity reduction and by the dedicated comparison routines.
    atoms = sorted(term.atoms(sp.sin, sp.cos, sp.tan), key=sp.default_sort_key)
    repl = {}
    marker = sp.Dummy("z")
    for atom in atoms:
        ratio = _rational_pi(atom.args[0])
        if ratio is None:
            continue
        try:
            poly = minpoly_for(atom, marker)
            if poly.degree() > budget:
                continue
            repl[atom] = to_number_field(atom)
        except EXACT_METHOD_ERRORS:
            continue
    return term.xreplace(repl)


def canonical_algebraic(term: sp.Expr, degree_limit: int | None = None) -> sp.Expr:
    """Return a compact exact representation of a closed algebraic constant."""
    if not stage_allowed(sp.sympify(term), "minpoly"):
        return sp.sympify(term)

    term = sp.sympify(term)
    info = classify_algebraic_expression(term)
    if not info.is_algebraic:
        return term
    budget = cfg.ALG_CANON_MAX_DEGREE if degree_limit is None else int(degree_limit)
    prepared = reduce_exact_trig(term, budget)
    prepared = normalize_radicals(simplify_root_sets(prepared))
    try:
        prepared = sp.sqrtdenest(sp.radsimp(prepared))
    except EXACT_METHOD_ERRORS:
        pass
    marker = sp.Dummy("z")
    try:
        poly = minpoly_for(prepared, marker)
        if poly.degree() > budget:
            return prepared
        value = to_number_field(prepared)
        if _field_zero(value):
            return sp.Integer(0)
        if isinstance(value, sp.AlgebraicNumber):
            root = value.to_root(radicals=True, minpoly=value.minpoly_of_element())
            return root
        return value
    except EXACT_METHOD_ERRORS:
        return prepared

