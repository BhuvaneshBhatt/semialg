"""Exact normalization of radicals and related algebraic root sets."""

from __future__ import annotations

from collections.abc import Sequence
from itertools import combinations

import sympy as sp

from . import _config as cfg
from ._cost import within_budget
from ._errors import EXACT_METHOD_ERRORS


def _root_groups(term: sp.Expr) -> dict[tuple, list[sp.Expr]]:
    groups: dict[tuple, list[sp.Expr]] = {}
    for node in term.atoms(sp.CRootOf):
        root_poly = node.poly
        coeff_key = tuple(root_poly.all_coeffs())
        key = (coeff_key, int(root_poly.degree()))
        groups.setdefault(key, []).append(node)
    return groups


def _sym_poly_value(root_poly: sp.Poly, order: int) -> sp.Expr:
    coeffs = root_poly.all_coeffs()
    lead = coeffs[0]
    return sp.cancel(((-1) ** order) * coeffs[order] / lead)


def _subset_sym_values(
    root_poly: sp.Poly, omitted: Sequence[sp.Expr], order: int
) -> list[sp.Expr]:
    """Elementary symmetric values for selected roots via their complement."""
    omitted = tuple(omitted)
    comp = [sp.Integer(1)]
    for degree in range(1, order + 1):
        if degree > len(omitted):
            comp.append(sp.Integer(0))
        else:
            total = sum(
                (sp.Mul(*pick) for pick in combinations(omitted, degree)), sp.Integer(0)
            )
            comp.append(total)
    chosen = [sp.Integer(1)]
    for degree in range(1, order + 1):
        total = _sym_poly_value(root_poly, degree)
        correction = sum(
            (chosen[degree - idx] * comp[idx] for idx in range(1, degree + 1)),
            sp.Integer(0),
        )
        chosen.append(sp.cancel(total - correction))
    return chosen


def _rewrite_root_group(term: sp.Expr, roots: list[sp.Expr]) -> sp.Expr:
    """Rewrite a majority root subset through complementary symmetric data."""
    root_poly = roots[0].poly
    degree = int(root_poly.degree())
    selected = {int(root.index): root for root in roots}
    count = len(selected)
    if count <= degree // 2:
        return term
    all_roots = [sp.CRootOf(root_poly.as_expr(), pos) for pos in range(degree)]
    chosen = [selected[pos] for pos in sorted(selected)]
    omitted = [root for pos, root in enumerate(all_roots) if pos not in selected]
    markers = tuple(sp.Dummy(f"r{pos}") for pos in range(count))
    replaced = term.xreplace(dict(zip(chosen, markers)))
    try:
        reduced = sp.cancel(replaced)
        top, bottom = sp.fraction(reduced)
        top_sym, top_rest, top_map = sp.symmetrize(top, markers, formal=True)
        bot_sym, bot_rest, bot_map = sp.symmetrize(bottom, markers, formal=True)
        if top_rest != 0 or bot_rest != 0:
            return term
        max_order = max(len(top_map), len(bot_map))
        sym_vals = _subset_sym_values(root_poly, omitted, max_order)
        top_sub = {tag: sym_vals[pos + 1] for pos, (tag, _) in enumerate(top_map)}
        bot_sub = {tag: sym_vals[pos + 1] for pos, (tag, _) in enumerate(bot_map)}
        new_top = sp.cancel(top_sym.subs(top_sub))
        new_bottom = sp.cancel(bot_sym.subs(bot_sub))
        if new_bottom == 0:
            return term
        candidate = sp.cancel(new_top / new_bottom)
        if len(candidate.atoms(sp.CRootOf)) >= count:
            return term
        return candidate
    except EXACT_METHOD_ERRORS:
        return term


def simplify_root_sets(term: sp.Expr) -> sp.Expr:
    """Reduce symmetric expressions involving a majority of common roots.

    When more than half of a polynomial's roots occur, elementary symmetric
    identities rewrite the selected roots through the smaller complementary
    root set. Complete sets therefore collapse entirely to coefficient data.
    """
    if not within_budget(sp.sympify(term)):
        return sp.sympify(term)

    term = sp.sympify(term)
    current = term
    for roots in _root_groups(current).values():
        current = _rewrite_root_group(current, roots)
    return current


def _prime_power(base: sp.Rational, exponent: sp.Rational) -> sp.Expr:
    """Express a positive rational power on a canonical prime-exponent basis."""
    pieces: list[sp.Expr] = []
    powers = dict(sp.factorint(int(base.p)))
    for prime, count in sp.factorint(int(base.q)).items():
        powers[prime] = powers.get(prime, 0) - count
    for prime in sorted(powers):
        scaled = sp.Rational(powers[prime]) * exponent
        whole = int(sp.floor(scaled))
        frac = scaled - whole
        if whole:
            pieces.append(sp.Integer(prime) ** whole)
        if frac:
            pieces.append(sp.Pow(sp.Integer(prime), frac))
    return sp.Mul(*pieces) if pieces else sp.Integer(1)


def _quadratic_radical_parts(base: sp.Expr):
    """Return ``(a, b, c)`` for ``a + b*sqrt(c)`` with rational data."""
    if not base.is_Add or len(base.args) != 2:
        return None
    rational = None
    radical = None
    for part in base.args:
        if part.is_Rational:
            rational = sp.Rational(part)
            continue
        coeff, rest = part.as_coeff_Mul(rational=True)
        if (
            coeff.is_Rational
            and rest.is_Pow
            and rest.exp == sp.Rational(1, 2)
            and rest.base.is_Rational
            and rest.base.is_positive
        ):
            radical = (sp.Rational(coeff), sp.Rational(rest.base))
    if rational is None or radical is None:
        return None
    return rational, radical[0], radical[1]


def denest_quadratic_radical(term: sp.Expr) -> sp.Expr:
    """Cheap exact denesting for ``sqrt(a + b*sqrt(c))`` rational shapes."""
    term = sp.sympify(term)
    if int(term.count_ops()) > cfg.RADICAL_DENEST_MAX_OPS:
        return term
    if not (term.is_Pow and term.exp == sp.Rational(1, 2)):
        return term
    parts = _quadratic_radical_parts(term.base)
    if parts is None:
        return term
    a, b, c = parts
    disc = sp.cancel(a * a - b * b * c)
    if not (disc.is_Rational and disc >= 0):
        return term
    root_disc = sp.sqrt(disc)
    if root_disc.is_Rational is not True:
        return term
    u = sp.cancel((a + root_disc) / 2)
    v = sp.cancel((a - root_disc) / 2)
    if not (u.is_Rational and v.is_Rational and u >= 0 and v >= 0):
        return term
    sign = 1 if b > 0 else -1
    candidate = sp.sqrt(u) + sign * sp.sqrt(v)
    # The nonnegative principal root is selected because u >= v >= 0.
    return candidate


def normalize_radicals(term: sp.Expr) -> sp.Expr:
    """Canonicalize rational-base radicals onto shared prime generators."""
    if not within_budget(sp.sympify(term)):
        return sp.sympify(term)

    term = sp.sympify(term)
    if not term.args:
        return term
    rebuilt = term.func(*(normalize_radicals(part) for part in term.args))
    denested = denest_quadratic_radical(rebuilt)
    if denested != rebuilt:
        return normalize_radicals(denested)
    if (
        rebuilt.is_Pow
        and rebuilt.base.is_Rational
        and rebuilt.base.is_positive
        and rebuilt.exp.is_Rational
    ):
        return sp.cancel(_prime_power(rebuilt.base, rebuilt.exp))
    if rebuilt.is_number and rebuilt.is_algebraic is True:
        try:
            return sp.sqrtdenest(rebuilt)
        except EXACT_METHOD_ERRORS:
            return rebuilt
    return rebuilt
