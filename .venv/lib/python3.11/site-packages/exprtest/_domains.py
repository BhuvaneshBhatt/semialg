"""Cheap number-domain facts and transcendence classification."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

import sympy as sp

from ._assumptions import assumption_facts, normalize_assumptions
from ._cost import ExactBudget, stage_allowed, within_budget
from ._errors import RECOVERABLE_POLY_ERRORS, RECOVERABLE_SYMPY_ERRORS
from ._function_facts import PROPERTY_FUNCS, function_property

_alg_var = sp.Dummy("a")
Algebraics = sp.ConditionSet(_alg_var, sp.Q.algebraic(_alg_var), sp.S.Complexes)


class NumberKind(Enum):
    """Proof status for the arithmetic nature of an exact numeric expression."""

    ALGEBRAIC = "algebraic"
    TRANSCENDENTAL = "transcendental-proven"
    UNKNOWN = "unknown"


@dataclass(frozen=True)
class DomainFacts:
    """Cheap theorem/property facts used by membership and branch decisions."""

    integer: bool | None = None
    rational: bool | None = None
    algebraic: bool | None = None
    real: bool | None = None
    complex: bool | None = None
    positive: bool | None = None
    negative: bool | None = None
    nonzero: bool | None = None
    finite: bool | None = None
    root_of_unity: bool | None = None


def _ask(prop, assumptions) -> bool | None:
    facts = assumption_facts(assumptions)
    try:
        if facts:
            with sp.assumptions.assuming(*facts):
                ans = sp.ask(prop)
        else:
            ans = sp.ask(prop)
        return None if ans is None else bool(ans)
    except RECOVERABLE_SYMPY_ERRORS:
        return None


def _alg_atom(term: sp.Expr) -> bool:
    if term.is_Rational or term in (sp.I, sp.GoldenRatio):
        return True
    return bool(isinstance(term, sp.AlgebraicNumber) or getattr(term, "is_CRootOf", False))


def _cheap_alg(term: sp.Expr, budget: ExactBudget) -> bool:
    if _alg_atom(term):
        return True
    if term.is_Add or term.is_Mul:
        return all(_cheap_alg(arg, budget) for arg in term.args)
    if term.is_Pow and term.exp.is_Rational:
        if term.exp.is_negative and term.base.is_zero is True:
            return False
        return _cheap_alg(term.base, budget)
    if term.func in (sp.sin, sp.cos, sp.tan) and len(term.args) == 1:
        arg = term.args[0]
        coeff = arg.coeff(sp.pi)
        return bool(arg == coeff * sp.pi and coeff.is_Rational)
    return False


def _known_nonzero_alg(term: sp.Expr, budget: ExactBudget) -> bool:
    if not _cheap_alg(term, budget):
        return False
    if term.is_zero is False:
        return True
    if term.is_Rational:
        return term != 0
    return False


def _gelfond_schneider_pow(term: sp.Expr, budget: ExactBudget) -> bool:
    """Recognize algebraic-base, algebraic-irrational powers as transcendental."""
    if not term.is_Pow:
        return False
    base, exponent = term.args
    if not (_cheap_alg(base, budget) and _cheap_alg(exponent, budget)):
        return False
    if exponent.is_rational is not False:
        return False
    if not _known_nonzero_alg(base, budget):
        return False
    return _known_nonzero_alg(base - 1, budget)


def _trans_atom(term: sp.Expr, budget: ExactBudget) -> bool:
    """Recognize single constants known transcendental by cheap theorems."""
    if term in (sp.pi, sp.E):
        return True
    if _gelfond_schneider_pow(term, budget):
        return True
    if term.func is sp.exp and len(term.args) == 1:
        return _known_nonzero_alg(term.args[0], budget)
    if term.func is sp.log and len(term.args) == 1:
        arg = term.args[0]
        return _cheap_alg(arg, budget) and arg.is_zero is False and arg != 1
    if term.func in (sp.sin, sp.cos, sp.tan) and len(term.args) == 1:
        arg = term.args[0]
        coeff = arg.coeff(sp.pi)
        if arg == coeff * sp.pi and coeff.is_Rational:
            return False
        return _known_nonzero_alg(arg, budget)
    return False


def _rat_parts(term: sp.Expr, atom: sp.Expr, var: sp.Symbol,
               budget: ExactBudget):
    """Build a rational function in one transcendental atom structurally."""
    if term == atom:
        return sp.Poly(var, var, domain=sp.QQ), sp.Poly(1, var, domain=sp.QQ)
    if _cheap_alg(term, budget):
        # Polynomial coefficients may live in an algebraic extension. Use EX;
        # no algebraic simplifier or primitive element is invoked here.
        return sp.Poly(term, var, domain="EX"), sp.Poly(1, var, domain="EX")
    if term.is_Add:
        num = sp.Poly(0, var, domain="EX")
        den = sp.Poly(1, var, domain="EX")
        for arg in term.args:
            part = _rat_parts(arg, atom, var, budget)
            if part is None:
                return None
            pn, pd = part
            num = num * pd + pn * den
            den = den * pd
            if len(num.terms()) > budget.max_poly_terms or len(den.terms()) > budget.max_poly_terms:
                return None
        return num, den
    if term.is_Mul:
        num = sp.Poly(1, var, domain="EX")
        den = sp.Poly(1, var, domain="EX")
        for arg in term.args:
            part = _rat_parts(arg, atom, var, budget)
            if part is None:
                return None
            pn, pd = part
            num *= pn
            den *= pd
            if len(num.terms()) > budget.max_poly_terms or len(den.terms()) > budget.max_poly_terms:
                return None
        return num, den
    if term.is_Pow and term.exp.is_Integer:
        part = _rat_parts(term.base, atom, var, budget)
        if part is None:
            return None
        num, den = part
        power = int(term.exp)
        if power < 0:
            num, den = den, num
            power = -power
        if power.bit_length() > budget.max_pow_bits:
            return None
        return num ** power, den ** power
    return None


def rational_transcendental_kind(term: sp.Expr,
                                  budget: ExactBudget | None = None) -> NumberKind:
    """Classify a rational function of one known transcendental constant.

    A nonzero polynomial with algebraic coefficients cannot vanish at a
    transcendental argument. After exact polynomial gcd cancellation, a
    nonconstant rational function is transcendental; a constant one is
    algebraic. Unsupported forms return ``UNKNOWN``.
    """
    budget = budget or ExactBudget()
    term = sp.sympify(term)
    if term.free_symbols or not stage_allowed(term, "minpoly", budget):
        return NumberKind.UNKNOWN
    atoms = []
    for node in sp.preorder_traversal(term):
        if (
            _trans_atom(node, budget)
            and not any(node.has(prior) and node != prior for prior in atoms)
            and node not in atoms
        ):
            atoms.append(node)
    # Keep this theorem path deliberately one-generator and cheap.
    maximal = [a for a in atoms if not any(a != b and b.has(a) for b in atoms)]
    if len(maximal) != 1:
        return NumberKind.UNKNOWN
    atom = maximal[0]
    var = sp.Dummy("t")
    try:
        parts = _rat_parts(term, atom, var, budget)
        if parts is None:
            return NumberKind.UNKNOWN
        num, den = parts
        if den.is_zero:
            return NumberKind.UNKNOWN
        if num.is_zero:
            # Zero is algebraic, but leave cancellation to the zero oracle;
            # this classifier avoids deriving algebraicity from two
            # transcendental terms that happen to cancel.
            return NumberKind.UNKNOWN
        common = sp.polys.polytools.gcd(num, den)
        if common.degree() > 0:
            num = num.exquo(common)
            den = den.exquo(common)
        if num.degree() == 0 and den.degree() == 0:
            return NumberKind.ALGEBRAIC
        return NumberKind.TRANSCENDENTAL
    except RECOVERABLE_POLY_ERRORS:
        return NumberKind.UNKNOWN


def number_kind(term: sp.Expr, budget: ExactBudget | None = None) -> NumberKind:
    """Cheap sufficient algebraic/transcendental classifier."""
    budget = budget or ExactBudget()
    term = sp.sympify(term)
    if term.free_symbols or not within_budget(term, budget):
        return NumberKind.UNKNOWN
    if term.has(sp.Float, sp.nan, sp.zoo, sp.oo, -sp.oo):
        return NumberKind.UNKNOWN
    if _cheap_alg(term, budget):
        return NumberKind.ALGEBRAIC
    if _trans_atom(term, budget):
        return NumberKind.TRANSCENDENTAL

    rat_kind = rational_transcendental_kind(term, budget)
    if rat_kind is not NumberKind.UNKNOWN:
        return rat_kind

    if term.is_Add:
        trans = 0
        for arg in term.args:
            kind = number_kind(arg, budget)
            if kind is NumberKind.TRANSCENDENTAL:
                trans += 1
            elif kind is not NumberKind.ALGEBRAIC:
                return NumberKind.UNKNOWN
        if trans == 1:
            return NumberKind.TRANSCENDENTAL

    if term.is_Mul:
        trans = 0
        for arg in term.args:
            kind = number_kind(arg, budget)
            if kind is NumberKind.TRANSCENDENTAL:
                trans += 1
            elif kind is NumberKind.ALGEBRAIC:
                if arg.is_zero is True:
                    return NumberKind.ALGEBRAIC
                if arg.is_zero is not False and not arg.is_Rational:
                    return NumberKind.UNKNOWN
            else:
                return NumberKind.UNKNOWN
        if trans == 1:
            return NumberKind.TRANSCENDENTAL
    return NumberKind.UNKNOWN


def _root_unity(term: sp.Expr) -> bool | None:
    if term in (sp.Integer(1), sp.Integer(-1), sp.I, -sp.I):
        return True
    if term.func is sp.exp and len(term.args) == 1:
        arg = term.args[0]
        coeff = arg.coeff(sp.pi)
        if arg == coeff * sp.pi:
            ratio = coeff / (2 * sp.I)
            if ratio.is_Rational:
                return True
    return None


def domain_facts(term: sp.Expr, assumptions=True,
                 budget: ExactBudget | None = None) -> DomainFacts:
    """Collect inexpensive exact domain/sign facts without simplification."""
    budget = budget or ExactBudget()
    assumptions = normalize_assumptions(assumptions)
    term = sp.sympify(term)
    kind = number_kind(term, budget) if not term.free_symbols else NumberKind.UNKNOWN
    integer = bool(term.is_integer) if term.is_integer is not None else _ask(sp.Q.integer(term), assumptions)
    rational = bool(term.is_rational) if term.is_rational is not None else _ask(sp.Q.rational(term), assumptions)
    real = bool(term.is_real) if term.is_real is not None else _ask(sp.Q.real(term), assumptions)
    complex_ = bool(term.is_complex) if term.is_complex is not None else _ask(sp.Q.complex(term), assumptions)
    positive = bool(term.is_positive) if term.is_positive is not None else _ask(sp.Q.positive(term), assumptions)
    negative = bool(term.is_negative) if term.is_negative is not None else _ask(sp.Q.negative(term), assumptions)
    if term.func in PROPERTY_FUNCS:
        if integer is None:
            integer = function_property(term, "integer", assumptions)
        if real is None:
            real = function_property(term, "real", assumptions)
        if positive is None:
            positive = function_property(term, "positive", assumptions)
        if negative is None:
            negative = function_property(term, "negative", assumptions)
    finite = bool(term.is_finite) if term.is_finite is not None else _ask(sp.Q.finite(term), assumptions)
    zero = term.is_zero
    nonzero = (not bool(zero)) if zero is not None else None
    if nonzero is None and (positive is True or negative is True):
        nonzero = True
    if kind is NumberKind.TRANSCENDENTAL:
        nonzero = True
    algebraic = True if kind is NumberKind.ALGEBRAIC else (False if kind is NumberKind.TRANSCENDENTAL else None)
    return DomainFacts(integer, rational, algebraic, real, complex_, positive,
                       negative, nonzero, finite, _root_unity(term))



def is_integer(term: sp.Expr, assumptions=True) -> bool | None:
    """Return whether *term* is provably an integer.

    The result is ``True`` or ``False`` when the inexpensive domain engine can
    prove membership or nonmembership, and ``None`` otherwise.
    """
    return domain_facts(term, assumptions).integer


def is_rational(term: sp.Expr, assumptions=True) -> bool | None:
    """Return whether *term* is provably rational, or ``None`` if unknown."""
    return domain_facts(term, assumptions).rational


def is_real(term: sp.Expr, assumptions=True) -> bool | None:
    """Return whether *term* is provably real, or ``None`` if unknown."""
    return domain_facts(term, assumptions).real


def is_algebraic(term: sp.Expr, assumptions=True) -> bool | None:
    """Return whether *term* is provably algebraic, or ``None`` if unknown."""
    return domain_facts(term, assumptions).algebraic


def is_prime(term: sp.Expr, assumptions=True) -> bool | None:
    """Return whether *term* is provably a positive prime integer.

    Exact integer inputs use SymPy's deterministic integer primality test.
    Symbolic inputs use inexpensive assumptions only; this function does not
    factor general symbolic expressions.
    """
    term = sp.sympify(term)
    assumptions = normalize_assumptions(assumptions)
    if term.is_Integer:
        return bool(sp.ntheory.primetest.isprime(int(term)))
    facts = domain_facts(term, assumptions)
    if facts.integer is False or facts.rational is False or facts.algebraic is False:
        return False
    ans = _ask(sp.Q.prime(term), assumptions)
    return ans


def _direct_member(term: sp.Expr, set_, assumptions) -> bool | None:
    facts = domain_facts(term, assumptions)
    if set_ == sp.S.Integers:
        return facts.integer
    if set_ == sp.S.Rationals:
        return facts.rational
    if set_ == sp.S.Reals:
        return facts.real
    if set_ == sp.S.Complexes:
        return facts.complex
    if set_ == Algebraics:
        return facts.algebraic
    return None


def _strip_member_parts(term: sp.Expr, set_, assumptions) -> sp.Expr:
    if not (term.is_Add or term.is_Mul):
        return term
    kept = []
    removed = []
    for arg in term.args:
        if _direct_member(arg, set_, assumptions) is True:
            removed.append(arg)
        else:
            kept.append(arg)
    if not removed or not kept:
        return term
    if term.is_Mul:
        factor = sp.Mul(*removed)
        if set_ == sp.S.Integers:
            if factor not in (sp.Integer(1), sp.Integer(-1)):
                return term
        elif set_ in (sp.S.Rationals, sp.S.Reals) or set_ == Algebraics:
            if domain_facts(factor, assumptions).nonzero is not True:
                return term
        else:
            return term
    return term.func(*kept)


def ElementOf(elem: sp.Expr, set_, assumptions=True):
    """Simplify exact membership in a standard numeric domain."""
    assumptions = normalize_assumptions(assumptions)
    term = sp.sympify(elem)
    ans = _direct_member(term, set_, assumptions)
    if ans is not None:
        return ans
    if set_ == Algebraics and term.is_Pow and term.exp.is_Rational:
        facts = domain_facts(term.base, assumptions)
        if term.exp.is_positive or facts.nonzero is True:
            return ElementOf(term.base, set_, assumptions)
    reduced = _strip_member_parts(term, set_, assumptions)
    if reduced != term:
        ans = _direct_member(reduced, set_, assumptions)
        if ans is not None:
            return ans
        term = reduced
    return sp.Contains(term, set_, evaluate=False)
