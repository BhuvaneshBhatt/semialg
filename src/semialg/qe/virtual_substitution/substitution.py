from __future__ import annotations

from collections.abc import Sequence

import sympy as sp
from sympy.core.relational import (
    Relational,
)
from sympy.logic.boolalg import And as SymAnd
from sympy.logic.boolalg import Or as SymOr

from .atoms import (
    _and,
    _canonical_atom,
    _coefficient_list,
    _iter_atoms,
    _or,
    _polynomial_degree,
    _relation,
)
from .types import VirtualSubstitutionError, _QuadraticPoint


def _multiply_radical_pairs(
    left: tuple[sp.Expr, sp.Expr],
    right: tuple[sp.Expr, sp.Expr],
    radical: sp.Expr,
) -> tuple[sp.Expr, sp.Expr]:
    """Multiply ``a + b*r`` and ``c + d*r`` with ``r**2=radical``."""

    left_constant, left_radical = left
    right_constant, right_radical = right
    return (
        sp.expand(left_constant * right_constant + left_radical * right_radical * radical),
        sp.expand(left_constant * right_radical + left_radical * right_constant),
    )


def _evaluate_polynomial_at_quadratic_point(
    polynomial: sp.Expr,
    variable: sp.Symbol,
    point: _QuadraticPoint,
) -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    """Evaluate a polynomial at a quadratic virtual-substitution point.

    For degree ``k`` this computes ``denominator**k * f(point)`` as
    ``constant + radical_coeff*sqrt(radical)``. Keeping the computation in the
    two-dimensional basis ``1, sqrt(radical)`` is faster and less fragile than
    expanding with a placeholder and rewriting powers of that placeholder.
    """

    degree = max(_polynomial_degree(polynomial, variable), 0)
    coeffs = _coefficient_list(polynomial, variable, degree + 1)
    numerator_pair = (sp.expand(point.constant), sp.expand(point.radical_sign))
    powers: list[tuple[sp.Expr, sp.Expr]] = [(sp.Integer(1), sp.Integer(0))]
    for _ in range(degree):
        powers.append(_multiply_radical_pairs(powers[-1], numerator_pair, point.radical))

    constant = sp.Integer(0)
    radical_coeff = sp.Integer(0)
    for power, coeff in enumerate(coeffs):
        denominator_scale = point.denominator ** (degree - power)
        constant += coeff * denominator_scale * powers[power][0]
        radical_coeff += coeff * denominator_scale * powers[power][1]
    denominator_sign = point.denominator if degree % 2 else sp.Integer(1)
    return sp.expand(constant), sp.expand(radical_coeff), denominator_sign


def substitute_quadratic_root_in_atom(
    atom: Relational, variable: sp.Symbol, point: _QuadraticPoint
) -> sp.Expr:
    """Substitute a quadratic algebraic root into one relation without radicals."""

    polynomial, operator = _canonical_atom(atom)
    constant, radical_coeff, denominator_sign = _evaluate_polynomial_at_quadratic_point(
        polynomial, variable, point
    )
    norm = sp.expand(constant**2 - radical_coeff**2 * point.radical)

    if point.radical_sign == 0:
        if operator == "=":
            return sp.Eq(constant, 0)
        if operator == "!=":
            return sp.Ne(constant, 0)
        if operator == "<":
            return constant * denominator_sign < 0
        if operator == "<=":
            return constant * denominator_sign <= 0

    if operator == "=":
        return _and(constant * radical_coeff <= 0, sp.Eq(norm, 0))
    if operator == "!=":
        return _or(constant * radical_coeff > 0, sp.Ne(norm, 0))
    if operator == "<":
        return _or(
            _and(constant * denominator_sign < 0, norm > 0),
            _and(
                radical_coeff * denominator_sign <= 0,
                _or(constant * denominator_sign < 0, norm < 0),
            ),
        )
    if operator == "<=":
        return _or(
            _and(constant * denominator_sign <= 0, norm >= 0),
            _and(radical_coeff * denominator_sign <= 0, norm <= 0),
        )
    raise VirtualSubstitutionError(f"unsupported relation operator: {operator!r}")


def substitute_quadratic_root(
    formula: sp.Expr, variable: sp.Symbol, point: _QuadraticPoint
) -> sp.Expr:
    """Apply exact virtual substitution at a root candidate."""

    if formula is True or formula == sp.true:
        return sp.true
    if formula is False or formula == sp.false:
        return sp.false
    if isinstance(formula, Relational):
        return substitute_quadratic_root_in_atom(formula, variable, point)
    if isinstance(formula, SymAnd):
        return _and(*(substitute_quadratic_root(arg, variable, point) for arg in formula.args))
    if isinstance(formula, SymOr):
        return _or(*(substitute_quadratic_root(arg, variable, point) for arg in formula.args))
    raise VirtualSubstitutionError(f"unsupported Boolean formula node: {formula!r}")


def _negative_side_formula(polynomial: sp.Expr, variable: sp.Symbol, side: int) -> sp.Expr:
    """Return the virtual sign condition for ``polynomial < 0`` at root+side*eps."""

    polynomial = sp.expand(polynomial)
    if variable not in polynomial.free_symbols:
        return polynomial < 0
    derivative = sp.diff(polynomial, variable)
    return _or(
        polynomial < 0,
        _and(sp.Eq(polynomial, 0), _negative_side_formula(side * derivative, variable, side)),
    )


def substitute_perturbed_root_atom(
    atom: Relational,
    variable: sp.Symbol,
    point: _QuadraticPoint,
    side: int,
) -> sp.Expr:
    polynomial, operator = _canonical_atom(atom)
    coeffs = _coefficient_list(polynomial, variable, _polynomial_degree(polynomial, variable) + 1)
    if operator == "=":
        return _and(*(sp.Eq(coeff, 0) for coeff in coeffs))
    if operator == "!=":
        return _or(*(sp.Ne(coeff, 0) for coeff in coeffs))
    if operator == "<":
        return substitute_quadratic_root(
            _negative_side_formula(polynomial, variable, side), variable, point
        )
    if operator == "<=":
        identically_zero = _and(*(sp.Eq(coeff, 0) for coeff in coeffs))
        infinitesimally_negative = substitute_perturbed_root_atom(
            _relation(polynomial, "<"),
            variable,
            point,
            side,
        )
        return _or(identically_zero, infinitesimally_negative)
    raise VirtualSubstitutionError(f"unsupported relation operator: {operator!r}")


def substitute_perturbed_quadratic_root(
    formula: sp.Expr,
    variable: sp.Symbol,
    point: _QuadraticPoint,
    side: int = 1,
) -> sp.Expr:
    """Apply virtual substitution at ``point + side*epsilon``."""

    if formula is True or formula == sp.true:
        return sp.true
    if formula is False or formula == sp.false:
        return sp.false
    if isinstance(formula, Relational):
        return substitute_perturbed_root_atom(formula, variable, point, side)
    if isinstance(formula, SymAnd):
        return _and(
            *(
                substitute_perturbed_quadratic_root(arg, variable, point, side)
                for arg in formula.args
            )
        )
    if isinstance(formula, SymOr):
        return _or(
            *(
                substitute_perturbed_quadratic_root(arg, variable, point, side)
                for arg in formula.args
            )
        )
    raise VirtualSubstitutionError(f"unsupported Boolean formula node: {formula!r}")


def _infinity_negative_formula(coeffs: Sequence[sp.Expr], side: int) -> sp.Expr:
    pieces: list[sp.Expr] = []
    for degree in range(len(coeffs) - 1, -1, -1):
        coeff = coeffs[degree]
        higher_zero = [sp.Eq(coeffs[j], 0) for j in range(degree + 1, len(coeffs))]
        pieces.append(_and(*higher_zero, side**degree * coeff < 0))
    return _or(*pieces)


def substitute_infinity_in_atom(atom: Relational, variable: sp.Symbol, side: int) -> sp.Expr:
    polynomial, operator = _canonical_atom(atom)
    coeffs = _coefficient_list(polynomial, variable, _polynomial_degree(polynomial, variable) + 1)
    if operator == "=":
        return _and(*(sp.Eq(coeff, 0) for coeff in coeffs))
    if operator == "!=":
        return _or(*(sp.Ne(coeff, 0) for coeff in coeffs))
    if operator == "<":
        return _infinity_negative_formula(coeffs, side)
    if operator == "<=":
        return _or(
            substitute_infinity_in_atom(_relation(polynomial, "="), variable, side),
            substitute_infinity_in_atom(_relation(polynomial, "<"), variable, side),
        )
    raise VirtualSubstitutionError(f"unsupported relation operator: {operator!r}")


def substitute_infinity(formula: sp.Expr, variable: sp.Symbol, side: int = -1) -> sp.Expr:
    """Apply virtual substitution at ``side*infinity``."""

    if formula is True or formula == sp.true:
        return sp.true
    if formula is False or formula == sp.false:
        return sp.false
    if isinstance(formula, Relational):
        return substitute_infinity_in_atom(formula, variable, side)
    if isinstance(formula, SymAnd):
        return _and(*(substitute_infinity(arg, variable, side) for arg in formula.args))
    if isinstance(formula, SymOr):
        return _or(*(substitute_infinity(arg, variable, side) for arg in formula.args))
    raise VirtualSubstitutionError(f"unsupported Boolean formula node: {formula!r}")


def _root_candidates_from_polynomial(
    polynomial: sp.Expr, variable: sp.Symbol
) -> tuple[tuple[sp.Expr, _QuadraticPoint], ...]:
    degree = _polynomial_degree(polynomial, variable)
    if degree < 1:
        return ()
    if degree > 2:
        raise VirtualSubstitutionError(
            f"degree {degree} in {variable} exceeds the supported quadratic virtual-substitution fragment"
        )

    coeffs = _coefficient_list(polynomial, variable, 3)
    constant, linear, quadratic = coeffs
    linear_point = _QuadraticPoint(-constant, sp.Integer(0), sp.Integer(0), linear)
    linear_guard = sp.Ne(linear, 0)
    if degree == 1:
        return ((linear_guard, linear_point),)

    discriminant = sp.expand(linear**2 - 4 * quadratic * constant)
    quadratic_guard = _and(sp.Ne(quadratic, 0), discriminant >= 0)
    minus_root = _QuadraticPoint(-linear, sp.Integer(-1), discriminant, 2 * quadratic)
    plus_root = _QuadraticPoint(-linear, sp.Integer(1), discriminant, 2 * quadratic)
    degenerate_guard = _and(sp.Eq(quadratic, 0), sp.Ne(linear, 0))
    return (
        (degenerate_guard, linear_point),
        (quadratic_guard, minus_root),
        (quadratic_guard, plus_root),
    )


def _unique_polynomials(formula: sp.Expr, variable: sp.Symbol) -> tuple[sp.Expr, ...]:
    polynomials: list[sp.Expr] = []
    seen: set[str] = set()
    for atom in _iter_atoms(formula):
        polynomial, _ = _canonical_atom(atom)
        degree = _polynomial_degree(polynomial, variable)
        if degree > 2:
            raise VirtualSubstitutionError(
                f"degree {degree} in {variable} exceeds the supported quadratic virtual-substitution fragment"
            )
        if degree <= 0:
            continue
        key = sp.sstr(sp.expand(polynomial))
        if key not in seen:
            polynomials.append(sp.expand(polynomial))
            seen.add(key)
    return tuple(polynomials)
