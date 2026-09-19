"""Exact arithmetic for roots of unity and rational-angle trig constants."""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from math import gcd

import sympy as sp
from sympy.polys.numberfields import to_number_field

from . import _config as cfg
from ._cost import ExactBudget, stage_allowed
from ._errors import EXACT_METHOD_ERRORS
from ._memo import cyclo_for
from ._result import Verdict, ZeroClassification


@dataclass(frozen=True)
class CyclotomicForm:
    """Reduced element of ``QQ[z] / Phi_n(z)`` for a primitive root of unity."""

    order: int
    variable: sp.Symbol
    modulus: sp.Poly
    value: sp.Expr

    @property
    def is_zero(self) -> bool:
        return self.value == 0

    def as_expr(self) -> sp.Expr:
        """Return the exact value using ``exp(2*pi*I/order)`` as the root."""
        root = sp.exp(2 * sp.pi * sp.I / self.order)
        return sp.cancel(self.value.subs(self.variable, root))


def _lcm(a: int, b: int) -> int:
    return abs(a * b) // gcd(a, b) if a and b else 0


def _pi_ratio(arg: sp.Expr) -> sp.Rational | None:
    coeff = arg.coeff(sp.pi)
    if arg == coeff * sp.pi and coeff.is_Rational:
        return sp.Rational(coeff)
    return None


def _root_ratio(arg: sp.Expr) -> sp.Rational | None:
    coeff = arg.coeff(sp.pi)
    if arg != coeff * sp.pi:
        return None
    ratio = coeff / (2 * sp.I)
    return sp.Rational(ratio) if ratio.is_Rational else None


def _needed_order(term: sp.Expr) -> int | None:
    """Return a common root-of-unity order, or ``None`` if unsupported."""
    order = 4
    for node in sp.preorder_traversal(term):
        if node.func in (sp.sin, sp.cos, sp.tan) and len(node.args) == 1:
            ratio = _pi_ratio(node.args[0])
            if ratio is None:
                return None
            order = _lcm(order, 2 * int(ratio.q))
        elif node.func is sp.exp and len(node.args) == 1:
            ratio = _root_ratio(node.args[0])
            if ratio is None:
                return None
            order = _lcm(order, int(ratio.q))
    if order > cfg.CYCLOTOMIC_MAX_ORDER:
        return None
    return order


def _reduce(value: sp.Expr, var: sp.Symbol, mod: sp.Poly) -> sp.Expr:
    return sp.Poly(value, var, domain=sp.QQ).rem(mod).as_expr()


def _inverse(value: sp.Expr, var: sp.Symbol, mod: sp.Poly) -> sp.Expr | None:
    try:
        inv = sp.invert(sp.Poly(value, var, domain=sp.QQ), mod)
        return _reduce(inv.as_expr(), var, mod)
    except EXACT_METHOD_ERRORS:
        return None


def _power(value: sp.Expr, power: int, var: sp.Symbol, mod: sp.Poly) -> sp.Expr | None:
    if power < 0:
        value = _inverse(value, var, mod)
        if value is None:
            return None
        power = -power
    result = sp.Integer(1)
    base = _reduce(value, var, mod)
    while power:
        if power & 1:
            result = _reduce(result * base, var, mod)
        power >>= 1
        if power:
            base = _reduce(base * base, var, mod)
    return result


def _root_power(k: int, order: int, var: sp.Symbol, mod: sp.Poly) -> sp.Expr:
    return _power(var, k % order, var, mod) or sp.Integer(0)


def _trig_value(term: sp.Expr, order: int, var: sp.Symbol, mod: sp.Poly) -> sp.Expr | None:
    ratio = _pi_ratio(term.args[0])
    if ratio is None:
        return None
    scaled = sp.Rational(order) * ratio / 2
    if not scaled.is_Integer:
        return None
    k = int(scaled)
    pos = _root_power(k, order, var, mod)
    neg = _root_power(-k, order, var, mod)
    if term.func is sp.cos:
        return _reduce((pos + neg) / 2, var, mod)
    imag = _root_power(order // 4, order, var, mod)
    if term.func is sp.sin:
        inv_i = _inverse(2 * imag, var, mod)
        return None if inv_i is None else _reduce((pos - neg) * inv_i, var, mod)
    if term.func is sp.tan:
        den = _reduce(pos + neg, var, mod)
        inv_den = _inverse(den, var, mod)
        inv_i = _inverse(imag, var, mod)
        if inv_den is None or inv_i is None:
            return None
        return _reduce((pos - neg) * inv_i * inv_den, var, mod)
    return None


def _eval(term: sp.Expr, order: int, var: sp.Symbol, mod: sp.Poly) -> sp.Expr | None:
    if term.is_Rational:
        return term
    if term == sp.I:
        return _root_power(order // 4, order, var, mod)
    if term.func in (sp.sin, sp.cos, sp.tan):
        return _trig_value(term, order, var, mod)
    if term.func is sp.exp and len(term.args) == 1:
        ratio = _root_ratio(term.args[0])
        if ratio is None:
            return None
        scaled = sp.Rational(order) * ratio
        if not scaled.is_Integer:
            return None
        return _root_power(int(scaled), order, var, mod)
    if term.is_Add:
        acc = sp.Integer(0)
        for arg in term.args:
            part = _eval(arg, order, var, mod)
            if part is None:
                return None
            acc = _reduce(acc + part, var, mod)
        return acc
    if term.is_Mul:
        acc = sp.Integer(1)
        for arg in term.args:
            part = _eval(arg, order, var, mod)
            if part is None:
                return None
            acc = _reduce(acc * part, var, mod)
        return acc
    if term.is_Pow and term.exp.is_Integer:
        base = _eval(term.base, order, var, mod)
        if base is None:
            return None
        return _power(base, int(term.exp), var, mod)
    return None


@lru_cache(maxsize=1024)
def _form_cached(term: sp.Expr, limit: int) -> CyclotomicForm | None:
    order = _needed_order(term)
    if order is None or order > limit:
        return None
    var = sp.Dummy("z")
    try:
        mod = cyclo_for(order, var)
        if mod.degree() > cfg.CYCLOTOMIC_MAX_DEGREE:
            return None
        value = _eval(term, order, var, mod)
        if value is None:
            return None
        return CyclotomicForm(order, var, mod, _reduce(value, var, mod))
    except EXACT_METHOD_ERRORS:
        return None


def cyclotomic_form(term: sp.Expr, order_limit: int | None = None) -> CyclotomicForm | None:
    """Represent a supported exact constant in a cyclotomic quotient field."""
    budget = ExactBudget()
    term = sp.sympify(term)
    if not stage_allowed(term, "cyclotomic", budget):
        return None
    limit = min(cfg.CYCLOTOMIC_MAX_ORDER, budget.max_cyclo) if order_limit is None else min(int(order_limit), budget.max_cyclo)
    return _form_cached(term, limit)


def cyclotomic_zero_test(term: sp.Expr) -> ZeroClassification:
    """Prove zero or nonzero directly in a cyclotomic quotient when possible."""
    form = cyclotomic_form(term)
    if form is None:
        return ZeroClassification(Verdict.UNKNOWN, "cyclotomic", detail="expression is outside the cyclotomic subset")
    if form.is_zero:
        return ZeroClassification(
            Verdict.ZERO_PROVEN,
            "cyclotomic",
            detail=f"reduction modulo Phi_{form.order} is exactly zero",
            evidence="cyclotomic-polynomial-reduction",
        )
    return ZeroClassification(
        Verdict.NONZERO_PROVEN,
        "cyclotomic",
        detail=f"nonzero residue modulo Phi_{form.order}",
        evidence="cyclotomic-polynomial-reduction",
    )


def cyclotomic_sign(term: sp.Expr) -> int | None:
    """Return the exact sign of a supported real cyclotomic constant."""
    term = sp.sympify(term)
    if term.is_real is not True:
        return None
    form = cyclotomic_form(term)
    if form is None:
        return None
    if form.is_zero:
        return 0
    try:
        value = to_number_field(form.as_expr())
        exact = value.to_root(radicals=False) if isinstance(value, sp.AlgebraicNumber) else value
        pos = sp.ask(sp.Q.positive(exact))
        if pos is True:
            return 1
        neg = sp.ask(sp.Q.negative(exact))
        if neg is True:
            return -1
    except EXACT_METHOD_ERRORS:
        return None
    return None
