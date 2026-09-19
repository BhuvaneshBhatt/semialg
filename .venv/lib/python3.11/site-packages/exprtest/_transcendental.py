"""Exact structural rules for selected transcendental constants.

All transformations are branch-safe. Logarithms are expanded only across
positive factors, rational powers are extracted only from positive bases,
and closed exp/log inversions are used only where the principal branch is
unambiguous.
"""

from __future__ import annotations

from math import lcm

import sympy as sp

from . import _config as cfg
from ._assumptions import assumption_facts, normalize_assumptions
from ._cost import ExactBudget, stage_allowed, within_budget
from ._domains import domain_facts
from ._errors import EXACT_METHOD_ERRORS
from ._fast import quick_reduce
from ._negative_cache import log_rel_inapplicable
from ._result import Verdict, ZeroClassification


def _ask(prop, assumptions) -> bool | None:
    facts = assumption_facts(assumptions)
    try:
        if facts:
            with sp.assumptions.assuming(*facts):
                value = sp.ask(prop)
        else:
            value = sp.ask(prop)
        return None if value is None else bool(value)
    except EXACT_METHOD_ERRORS:
        return None


def _positive(term: sp.Expr, assumptions) -> bool:
    facts = domain_facts(term, assumptions)
    if facts.positive is True:
        return True
    if (not term.free_symbols and facts.algebraic is True
            and facts.real is True):
        from ._exact_constants import algebraic_sign

        return algebraic_sign(term) == 1
    return False


def _integer(term: sp.Expr, assumptions) -> bool:
    if term.is_integer is True:
        return True
    return _ask(sp.Q.integer(term), assumptions) is True


def _safe_log_once(term: sp.Expr, assumptions) -> sp.Expr:
    if not term.args:
        return term
    args = tuple(_safe_log_once(arg, assumptions) for arg in term.args)
    try:
        rebuilt = term.func(*args)
    except EXACT_METHOD_ERRORS:
        rebuilt = term
    if rebuilt.func is not sp.log or len(rebuilt.args) != 1:
        return rebuilt
    arg = rebuilt.args[0]
    if arg.is_Mul and all(_positive(part, assumptions) for part in arg.args):
        return sp.Add(*(sp.log(part) for part in arg.args))
    if arg.is_Pow and arg.exp.is_Rational and _positive(arg.base, assumptions):
        return arg.exp * sp.log(arg.base)
    return rebuilt


def normalize_logs(term: sp.Expr, assumptions=True) -> sp.Expr:
    """Apply only branch-safe logarithm expansion rules."""
    assumptions = normalize_assumptions(assumptions)
    current = sp.sympify(term)
    for _ in range(cfg.LOG_NORMALIZE_PASSES):
        updated = _safe_log_once(current, assumptions)
        updated = quick_reduce(updated)
        if updated == current:
            break
        current = updated
    return current


def _log_term(term: sp.Expr, assumptions) -> tuple[sp.Rational, sp.Expr] | None:
    """Extract ``q*log(a)`` with rational q and positive exact algebraic a."""
    coeff, rest = term.as_coeff_Mul(rational=True)
    if rest.func is not sp.log or len(rest.args) != 1 or not coeff.is_Rational:
        return None
    base = rest.args[0]
    if base.free_symbols or base.is_algebraic is not True or not _positive(base, assumptions):
        return None
    return sp.Rational(coeff), base


def _rational_log_vector(term: sp.Expr, assumptions, budget: ExactBudget):
    """Return the prime-exponent vector for a rational logarithmic relation."""
    if not stage_allowed(term, "factor-rational", budget):
        return None
    parts = term.args if term.is_Add else (term,)
    vector: dict[int, sp.Rational] = {}
    saw = False
    for part in parts:
        item = _log_term(part, assumptions)
        if item is None:
            return None
        coeff, base = item
        if not base.is_Rational or base <= 0:
            return None
        saw = True
        rat = sp.Rational(base)
        for value, sign in ((int(rat.p), 1), (int(rat.q), -1)):
            if value == 1:
                continue
            try:
                factors = sp.factorint(value)
            except EXACT_METHOD_ERRORS:
                return None
            for prime, power in factors.items():
                vector[int(prime)] = vector.get(int(prime), sp.Rational(0)) + coeff * sign * int(power)
    return {prime: power for prime, power in vector.items() if power != 0} if saw else None


def rational_log_test(term: sp.Expr, assumptions=True) -> ZeroClassification:
    """Decide rational log relations by prime-exponent vectors."""
    assumptions = normalize_assumptions(assumptions)
    term = normalize_logs(sp.sympify(term), assumptions)
    budget = ExactBudget()
    vector = _rational_log_vector(term, assumptions, budget)
    if vector is None:
        return ZeroClassification(Verdict.UNKNOWN, "rational-log", detail="not a bounded positive-rational log relation")
    if not vector:
        return ZeroClassification(Verdict.ZERO_PROVEN, "rational-log",
                                  detail="prime-exponent vector cancels exactly",
                                  evidence="rational-prime-factorization")
    return ZeroClassification(Verdict.NONZERO_PROVEN, "rational-log",
                              detail="prime-exponent vector is nonzero",
                              evidence="rational-prime-factorization")


def log_dependence_test(term: sp.Expr, assumptions=True) -> ZeroClassification:
    """Decide rational linear relations among logs of positive algebraic values.

    A relation ``sum(q_i*log(a_i))`` is multiplied by a common denominator,
    then reduced to the exact algebraic product ``prod(a_i**n_i)``. Since all
    logarithms are real on positive inputs, equality of that product to one is
    equivalent to the original logarithmic relation being zero.
    """
    budget = ExactBudget()
    if not within_budget(sp.sympify(term), budget):
        return ZeroClassification(Verdict.UNKNOWN, "log-dependence", detail="exact-method budget exceeded")

    assumptions = normalize_assumptions(assumptions)
    term = normalize_logs(sp.sympify(term), assumptions)
    if not log_rel_inapplicable(term):
        rational = rational_log_test(term, assumptions)
        if rational.verdict is not Verdict.UNKNOWN:
            return rational
    parts = term.args if term.is_Add else (term,)
    if len(parts) > budget.max_logs:
        return ZeroClassification(Verdict.UNKNOWN, "log-dependence", detail="log-term budget exceeded")
    items = [_log_term(part, assumptions) for part in parts]
    if not items or any(item is None for item in items):
        return ZeroClassification(Verdict.UNKNOWN, "log-dependence", detail="expression is not a supported rational log relation")

    pairs = [item for item in items if item is not None]

    # Exact structural cancellation is much cheaper than constructing a number
    # field and catches repeated logarithms after branch-safe normalization.
    coeffs: dict[sp.Expr, sp.Rational] = {}
    for coeff, base in pairs:
        coeffs[base] = coeffs.get(base, sp.Rational(0)) + coeff
    coeffs = {base: coeff for base, coeff in coeffs.items() if coeff != 0}
    if not coeffs:
        return ZeroClassification(
            Verdict.ZERO_PROVEN, "log-dependence",
            detail="branch-safe logarithm normalization cancels all exact bases",
            evidence="exact-log-base-vector",
        )

    scale = 1
    for coeff, _ in pairs:
        scale = lcm(scale, int(coeff.q))
    product = sp.Integer(1)
    for coeff, base in pairs:
        power = int(coeff * scale)
        product *= sp.Pow(base, power)

    relation = quick_reduce(product - 1)
    if relation == 0 or relation.is_zero is True:
        return ZeroClassification(
            Verdict.ZERO_PROVEN, "log-dependence",
            detail="exact multiplicative relation reduced structurally to one",
            evidence="exact-log-multiplicative-relation",
        )
    try:
        from ._algebraic import tower_algebraic_test
        tower = tower_algebraic_test(relation)
        if tower.verdict is Verdict.ZERO_PROVEN:
            return ZeroClassification(
                Verdict.ZERO_PROVEN, "log-dependence",
                detail="bounded algebraic-tower arithmetic proves the log arguments multiply to one",
                evidence="exact-log-tower-relation",
            )
        if tower.verdict is Verdict.NONZERO_PROVEN:
            return ZeroClassification(
                Verdict.NONZERO_PROVEN, "log-dependence",
                detail="bounded algebraic-tower arithmetic proves the log arguments do not multiply to one",
                evidence="exact-log-tower-relation",
            )
    except EXACT_METHOD_ERRORS:
        pass

    try:
        from ._exact_constants import canonical_algebraic, compare_algebraic

        reduced = canonical_algebraic(quick_reduce(product))
        cmp = compare_algebraic(reduced, sp.Integer(1))
        if cmp is not None:
            if cmp.order == 0:
                return ZeroClassification(
                    Verdict.ZERO_PROVEN,
                    "log-dependence",
                    detail="exact multiplicative relation among positive algebraic log arguments equals one",
                    evidence="exact-log-multiplicative-relation",
                )
            return ZeroClassification(
                Verdict.NONZERO_PROVEN,
                "log-dependence",
                detail="exact multiplicative relation among positive algebraic log arguments is not one",
                evidence="exact-log-multiplicative-relation",
            )
    except EXACT_METHOD_ERRORS:
        pass
    return ZeroClassification(Verdict.UNKNOWN, "log-dependence", detail="algebraic product relation was inconclusive")


def _normalize_exp_log_once(term: sp.Expr, assumptions) -> sp.Expr:
    if not term.args:
        return term
    args = tuple(_normalize_exp_log_once(arg, assumptions) for arg in term.args)

    if term.func is sp.exp and len(args) == 1:
        arg = args[0]
        # Remove exact integer multiples of 2*pi*I from an exponent.
        if arg.is_Add:
            kept = []
            changed = False
            for part in arg.args:
                q = sp.cancel(part / (2 * sp.pi * sp.I))
                if q * (2 * sp.pi * sp.I) == part and _integer(q, assumptions):
                    changed = True
                else:
                    kept.append(part)
            if changed:
                arg = sp.Add(*kept)
        if arg.func is sp.log and len(arg.args) == 1:
            base = arg.args[0]
            if (not base.free_symbols and base.is_finite is True
                    and base.is_zero is False):
                return base
            if _positive(base, assumptions):
                return base
        coeff, rest = arg.as_coeff_Mul(rational=True)
        if rest.func is sp.log and coeff.is_Rational:
            base = rest.args[0]
            if _positive(base, assumptions):
                return sp.Pow(base, coeff)
        if arg.is_Add:
            items = [_log_term(part, assumptions) for part in arg.args]
            if items and all(item is not None for item in items):
                product = sp.Integer(1)
                for coeff, base in (item for item in items if item is not None):
                    product *= sp.Pow(base, coeff)
                return quick_reduce(product)
        try:
            return sp.exp(arg, evaluate=False)
        except EXACT_METHOD_ERRORS:
            return term

    if term.func is sp.log and len(args) == 1:
        arg = args[0]
        if arg.func is sp.exp and len(arg.args) == 1:
            exponent = arg.args[0]
            # Principal log inverts exp exactly on the strip -pi < Im(z) <= pi.
            if exponent.is_real is True:
                return exponent
            imag = sp.im(exponent)
            coeff = imag.coeff(sp.pi)
            if imag == coeff * sp.pi and coeff.is_Rational and -1 < coeff <= 1:
                return exponent
        if arg.is_Pow and arg.exp.is_Rational and _positive(arg.base, assumptions):
            return arg.exp * sp.log(arg.base)
        if arg.is_Mul and all(_positive(part, assumptions) for part in arg.args):
            return sp.Add(*(sp.log(part) for part in arg.args))
        try:
            return sp.log(arg, evaluate=False)
        except EXACT_METHOD_ERRORS:
            return term

    if term.is_Mul:
        exp_args = [arg.args[0] for arg in args if arg.func is sp.exp]
        if len(exp_args) >= 2:
            others = [arg for arg in args if arg.func is not sp.exp]
            return sp.Mul(*others, sp.exp(sp.Add(*exp_args)))

    try:
        return term.func(*args)
    except EXACT_METHOD_ERRORS:
        return term


def normalize_exp_log(term: sp.Expr, assumptions=True) -> sp.Expr:
    """Normalize closed exp/log expressions using principal-branch-safe rules."""
    if not within_budget(sp.sympify(term)):
        return sp.sympify(term)

    assumptions = normalize_assumptions(assumptions)
    current = sp.sympify(term)
    for _ in range(cfg.LOG_NORMALIZE_PASSES):
        updated = _normalize_exp_log_once(current, assumptions)
        updated = normalize_logs(updated, assumptions)
        updated = quick_reduce(updated)
        if updated == current:
            break
        current = updated
    return current


def _periodic_zero(term: sp.Expr, assumptions) -> bool | None:
    if len(term.args) != 1:
        return None
    arg = term.args[0]
    coeff = arg.coeff(sp.pi)
    if arg != coeff * sp.pi:
        return None
    if term.func in (sp.sin, sp.tan):
        return _integer(coeff, assumptions)
    if term.func is sp.cos:
        return _integer(coeff - sp.Rational(1, 2), assumptions)
    if term.func in (sp.sinh, sp.tanh):
        return _integer(coeff / sp.I, assumptions)
    if term.func is sp.cosh:
        return _integer(coeff / sp.I - sp.Rational(1, 2), assumptions)
    return None


def transcendental_zero_test(term: sp.Expr, assumptions=True) -> ZeroClassification:
    """Apply exact zero-set, logarithmic, and branch-safe exp/log rules."""
    assumptions = normalize_assumptions(assumptions)
    term = sp.sympify(term)

    if not log_rel_inapplicable(term):
        rational = rational_log_test(term, assumptions)
        if rational.verdict is not Verdict.UNKNOWN:
            return rational

    normalized = normalize_exp_log(term, assumptions) if term.has(sp.log, sp.exp) else term
    if normalized != term:
        if normalized == 0:
            return ZeroClassification(
                Verdict.ZERO_PROVEN,
                "transcendental-normalize",
                detail="branch-safe exp/log normalization reduced expression to zero",
                evidence="exact-exp-log-identity",
            )
        if normalized.is_number and normalized.is_zero is False:
            # Only use exact SymPy knowledge here; unresolved exact constants
            # continue through the dedicated algebraic/cyclotomic stages.
            return ZeroClassification(
                Verdict.NONZERO_PROVEN,
                "transcendental-normalize",
                detail=f"branch-safe exp/log normalization reduced expression to {normalized}",
                evidence="exact-exp-log-identity",
            )

    if not log_rel_inapplicable(normalized):
        log_rel = log_dependence_test(normalized, assumptions)
        if log_rel.verdict is not Verdict.UNKNOWN:
            return log_rel

    periodic = _periodic_zero(normalized, assumptions)
    if periodic is True:
        return ZeroClassification(
            Verdict.ZERO_PROVEN,
            "transcendental-zero-set",
            detail="argument lies exactly in the function's zero lattice",
            evidence="exact-zero-set",
        )

    if len(normalized.args) == 1:
        arg = normalized.args[0]
        if normalized.func in (sp.asin, sp.atan, sp.asinh, sp.atanh, sp.Abs, sp.sign):
            if arg == 0:
                return ZeroClassification(
                    Verdict.ZERO_PROVEN,
                    "transcendental-zero-set",
                    detail=f"{normalized.func.__name__} has zero argument",
                    evidence="exact-zero-set",
                )
            if arg.is_number and arg.is_zero is False and normalized.is_finite is True:
                return ZeroClassification(
                    Verdict.NONZERO_PROVEN,
                    "transcendental-zero-set",
                    detail=f"finite {normalized.func.__name__} values can vanish only at zero",
                    evidence="exact-zero-set",
                )
        if normalized.func is sp.log and arg == 1:
            return ZeroClassification(
                Verdict.ZERO_PROVEN,
                "transcendental-zero-set",
                detail="principal logarithm vanishes at one",
                evidence="exact-zero-set",
            )
        if normalized.func in (sp.exp, sp.gamma) and normalized.is_finite is True:
            return ZeroClassification(
                Verdict.NONZERO_PROVEN,
                "transcendental-zero-set",
                detail=f"finite {normalized.func.__name__} values do not vanish",
                evidence="exact-zero-set",
            )

    return ZeroClassification(
        Verdict.UNKNOWN,
        "transcendental",
        detail="no exact branch-safe transcendental zero rule applied",
    )
