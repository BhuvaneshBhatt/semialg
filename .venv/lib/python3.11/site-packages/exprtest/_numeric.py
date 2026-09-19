"""Rigorous adaptive-precision numerics and assumption-aware witnesses."""

from __future__ import annotations

import random
from contextlib import contextmanager

try:
    import flint
except ImportError:
    flint = None
import sympy as sp

from . import _config as cfg
from ._assumptions import (
    assumptions_hold,
    equality_substitutions,
    normalize_assumptions,
)
from ._errors import RECOVERABLE_SYMPY_ERRORS
from ._features import expression_features
from ._result import Verdict, ZeroClassification


class UnsupportedNode(Exception):
    """Raised when an expression cannot be evaluated safely with Arb/Acb."""


_ACB_CALLABLE_UNARY = {
    sp.exp: "exp", sp.log: "log",
    sp.sin: "sin", sp.cos: "cos", sp.tan: "tan",
    sp.cot: "cot", sp.sec: "sec", sp.csc: "csc",
    sp.sinh: "sinh", sp.cosh: "cosh", sp.tanh: "tanh",
    sp.coth: "coth", sp.sech: "sech", sp.csch: "csch",
    sp.asin: "asin", sp.acos: "acos", sp.atan: "atan",
    sp.asinh: "asinh", sp.acosh: "acosh", sp.atanh: "atanh",
    sp.sign: "sgn", sp.conjugate: "conjugate", sp.arg: "arg",
}
_ACB_PROPERTY_UNARY = {sp.re: "real", sp.im: "imag"}


@contextmanager
def _flint_precision(prec_bits: int):
    old_prec = flint.ctx.prec
    flint.ctx.prec = prec_bits
    try:
        yield
    finally:
        flint.ctx.prec = old_prec


def _sympy_to_acb_ball(expr: sp.Expr, prec_bits: int, subs: dict):
    if expr in subs:
        return subs[expr]
    if expr.is_Symbol:
        raise UnsupportedNode(f"free symbol {expr} with no numeric value")
    if expr.is_Integer or expr.is_Rational:
        value = flint.fmpq(int(expr.p), int(expr.q)) if expr.is_Rational and not expr.is_Integer else int(expr)
        return flint.acb(value)
    if expr.is_Float:
        return flint.acb(str(expr))
    if expr is sp.pi:
        return flint.acb.pi()
    if expr is sp.E:
        return flint.acb(1).exp()
    if expr is sp.I:
        return flint.acb(0, 1)
    if expr is sp.EulerGamma:
        return flint.acb(flint.arb.const_euler())
    if expr is sp.Catalan:
        return flint.acb(flint.arb.const_catalan())
    if expr.is_Add:
        acc = flint.acb(0)
        for arg in expr.args:
            acc += _sympy_to_acb_ball(arg, prec_bits, subs)
        return acc
    if expr.is_Mul:
        acc = flint.acb(1)
        for arg in expr.args:
            acc *= _sympy_to_acb_ball(arg, prec_bits, subs)
        return acc
    if expr.is_Pow:
        base, exponent = expr.args
        bval = _sympy_to_acb_ball(base, prec_bits, subs)
        if exponent.is_Integer:
            return bval ** int(exponent)
        return bval ** _sympy_to_acb_ball(exponent, prec_bits, subs)
    if isinstance(expr, sp.Function):
        head = type(expr)
        if head is sp.Abs and len(expr.args) == 1:
            return flint.acb(abs(_sympy_to_acb_ball(expr.args[0], prec_bits, subs)))
        if head in _ACB_CALLABLE_UNARY and len(expr.args) == 1:
            arg = _sympy_to_acb_ball(expr.args[0], prec_bits, subs)
            return getattr(arg, _ACB_CALLABLE_UNARY[head])()
        if head in _ACB_PROPERTY_UNARY and len(expr.args) == 1:
            arg = _sympy_to_acb_ball(expr.args[0], prec_bits, subs)
            return flint.acb(getattr(arg, _ACB_PROPERTY_UNARY[head]))
        if head is sp.atan2:
            raise UnsupportedNode("atan2 requires quadrant-aware real-ball evaluation")
    raise UnsupportedNode(f"no Arb rule for {type(expr).__name__}")


def _log2_radius(rad: flint.arb) -> float:
    mantissa, exponent = rad.man_exp()
    mantissa = int(mantissa)
    exponent = int(exponent)
    if mantissa == 0:
        return float("-inf")
    return float(exponent + mantissa.bit_length())



def _precision_from_lower(lower) -> int:
    """Choose a conservative Arb starting precision from an exact lower bound."""
    try:
        value = sp.Rational(lower)
    except RECOVERABLE_SYMPY_ERRORS:
        return cfg.INITIAL_PREC_BITS
    if value <= 0 or value >= 1:
        return cfg.INITIAL_PREC_BITS
    num = abs(int(value.p))
    den = int(value.q)
    # ceil(log2(den/num)) using integer arithmetic, plus guard bits.
    bits = max(0, den.bit_length() - num.bit_length() + 1)
    target = bits + max(cfg.GUARD_BITS, cfg.ARB_HINT_MIN_GUARD_BITS)
    return min(cfg.ARB_HINT_MAX_START_BITS, max(cfg.INITIAL_PREC_BITS, target))



def _precision_from_expr(expr: sp.Expr) -> int:
    """Predict a modest starting precision from structural cancellation risk."""
    features = expression_features(sp.sympify(expr))
    extra = 0
    if sp.sympify(expr).is_Add:
        extra += min(96, max(0, features.max_int_bits // 2))
        extra += min(48, 4 * max(0, len(expr.args) - 2))
    if features.has_log_exp:
        extra += min(32, 4 * (features.logs + int(features.has_exp)))
    if features.depth > 10:
        extra += min(32, 2 * (features.depth - 10))
    return min(cfg.ARB_HINT_MAX_START_BITS, cfg.INITIAL_PREC_BITS + extra)


def _next_precision(prec: int, cur_log2_rad: float,
                    prev_prec: int | None, prev_log2_rad: float | None,
                    lower=None) -> int:
    """Predict the next useful Arb precision from observed enclosure shrinkage."""
    default = min(cfg.MAX_PREC_BITS + 1, max(prec + cfg.ARB_MIN_STEP_BITS, prec * 2))
    if prev_prec is None or prev_log2_rad is None:
        return default
    gain = prec - prev_prec
    shrink = prev_log2_rad - cur_log2_rad
    if gain <= 0 or shrink <= 0:
        return min(cfg.MAX_PREC_BITS + 1, prec * cfg.ARB_MAX_GROWTH_FACTOR)
    slope = shrink / gain
    if lower is not None:
        try:
            bound = sp.Rational(lower)
            if 0 < bound < 1:
                lower_bits = max(0, int(bound.q).bit_length() - abs(int(bound.p)).bit_length() + 1)
                # Radius log2 <= -lower_bits-guard should separate a nonzero value.
                target_rad = -lower_bits - cfg.ARB_TARGET_GUARD_BITS
                need_shrink = max(0.0, cur_log2_rad - target_rad)
                predicted = prec + int(need_shrink / max(slope, 0.1)) + cfg.ARB_MIN_STEP_BITS
                return min(cfg.MAX_PREC_BITS + 1,
                           max(prec + cfg.ARB_MIN_STEP_BITS,
                               min(predicted, prec * cfg.ARB_MAX_GROWTH_FACTOR)))
        except RECOVERABLE_SYMPY_ERRORS:
            pass
    # Good convergence needs only a measured step; poor convergence gets a
    # larger jump so we do not spend the budget on many nearly identical balls.
    factor = 2 if slope >= cfg.MIN_SHRINK_FRACTION else cfg.ARB_MAX_GROWTH_FACTOR
    return min(cfg.MAX_PREC_BITS + 1, max(prec + cfg.ARB_MIN_STEP_BITS, prec * factor))

def adaptive_precision_ball_test(expr: sp.Expr, subs: dict | None = None,
                                 separation_bound=None) -> ZeroClassification:
    """Escalate Arb precision and record quantitative zero-like evidence.

    When an exact positive lower bound for a nonzero value is already known,
    it is used only to choose a better starting precision. The enclosure still
    supplies the proof; the hint never changes the verdict by itself.
    """
    if flint is None:
        return ZeroClassification(Verdict.UNKNOWN, "numerical", detail="python-flint is unavailable")
    expr = sp.sympify(expr)
    if expr.has(sp.Float):
        return ZeroClassification(
            Verdict.UNKNOWN, "numerical",
            detail="inexact SymPy Float input is not promoted to a rigorous exact enclosure",
            evidence="inexact-input",
        )
    subs = subs or {}
    structural_prec = _precision_from_expr(expr)
    lower_prec = _precision_from_lower(separation_bound) if separation_bound is not None else cfg.INITIAL_PREC_BITS
    prec = min(cfg.ARB_HINT_MAX_START_BITS, max(structural_prec, lower_prec))
    prev_log2_rad: float | None = None
    prev_prec: int | None = None
    consecutive_shrinks = 0
    history = []

    while prec <= cfg.MAX_PREC_BITS:
        try:
            with _flint_precision(prec):
                val = _sympy_to_acb_ball(expr, prec, subs)
        except UnsupportedNode as exc:
            return ZeroClassification(Verdict.UNKNOWN, "numerical", detail=str(exc))
        if not val.is_finite:
            return ZeroClassification(Verdict.UNKNOWN, "numerical", detail="non-finite enclosure during evaluation")

        history.append(f"{prec} bits: {val}")
        if 0 not in val:
            return ZeroClassification(
                Verdict.NONZERO_PROVEN, "numerical", precision_bits=prec,
                detail=f"Arb ball {val} rigorously excludes zero",
                evidence="rigorous-enclosure", enclosure_history=tuple(history),
            )
        if val.rad() == 0:
            return ZeroClassification(
                Verdict.ZERO_PROVEN, "numerical", precision_bits=prec,
                detail=f"Arb ball {val} is exact and equals zero",
                evidence="exact-enclosure", enclosure_history=tuple(history),
            )

        cur_log2_rad = _log2_radius(val.rad())
        if prev_log2_rad is not None and prev_prec is not None:
            precision_gain = prec - prev_prec
            actual_shrink = prev_log2_rad - cur_log2_rad
            if precision_gain > 0 and actual_shrink >= cfg.MIN_SHRINK_FRACTION * precision_gain:
                consecutive_shrinks += 1
            else:
                consecutive_shrinks = 0
            if consecutive_shrinks >= cfg.REQUIRED_SHRINKS:
                return ZeroClassification(
                    Verdict.ZERO_UNPROVEN, "numerical", precision_bits=prec,
                    detail=(f"zero remains inside the enclosure while its radius tracks increasing precision "
                            f"for {consecutive_shrinks} consecutive escalations"),
                    evidence="precision-convergence", enclosure_history=tuple(history),
                )
        next_prec = _next_precision(prec, cur_log2_rad, prev_prec, prev_log2_rad, separation_bound)
        prev_log2_rad = cur_log2_rad
        prev_prec = prec
        prec = next_prec

    return ZeroClassification(
        Verdict.UNKNOWN, "numerical",
        detail=f"undetermined after escalating to {cfg.MAX_PREC_BITS} bits",
        evidence="inconclusive-enclosures", enclosure_history=tuple(history),
    )


def _random_exact_sample(symbol: sp.Symbol, rng: random.Random):
    """Generate an exact SymPy point, respecting intrinsic symbol assumptions."""
    if symbol.is_integer:
        if symbol.is_positive:
            return sp.Integer(rng.randint(1, 100))
        if symbol.is_nonnegative:
            return sp.Integer(rng.randint(0, 100))
        if symbol.is_negative:
            return sp.Integer(-rng.randint(1, 100))
        return sp.Integer(rng.randint(-100, 100))
    if symbol.is_real or symbol.is_positive or symbol.is_nonnegative or symbol.is_negative:
        num = rng.randint(1, 10_000)
        den = rng.randint(1, 97)
        value = sp.Rational(num, den)
        if symbol.is_negative:
            return -value
        if symbol.is_nonnegative and rng.random() < 0.05:
            return sp.Integer(0)
        if symbol.is_positive or symbol.is_nonnegative:
            return value
        return value if rng.random() < 0.5 else -value
    re = sp.Rational(rng.randint(-10_000, 10_000), rng.randint(1, 97))
    im = sp.Rational(rng.randint(-10_000, 10_000), rng.randint(1, 97))
    if rng.random() < 0.5:
        im = sp.Integer(0)
    return re + sp.I * im


def _sample_satisfying_assumptions(free_symbols, assumptions, rng: random.Random, max_attempts: int = 64):
    assumptions = normalize_assumptions(assumptions)
    symbols = tuple(sorted(free_symbols, key=str))
    equalities = equality_substitutions(assumptions)
    for _ in range(max_attempts):
        exact_subs = {symbol: _random_exact_sample(symbol, rng) for symbol in symbols}
        for symbol, value in equalities.items():
            if symbol in exact_subs:
                exact_subs[symbol] = sp.sympify(value).subs(exact_subs)
        if assumptions_hold(assumptions, exact_subs) is True:
            return exact_subs
    return None


def random_witness_nonzero_check(
    expr: sp.Expr,
    free_symbols,
    prec_bits: int = cfg.WITNESS_PREC_BITS,
    num_points: int = cfg.NUM_WITNESS_POINTS,
    assumptions=True,
    rng: random.Random | None = None,
) -> ZeroClassification | None:
    """Seek rigorous nonzero witness points that satisfy caller assumptions."""
    if flint is None:
        return None

    rng = rng or random.Random()
    for attempt in range(1, num_points + 1):
        exact_subs = _sample_satisfying_assumptions(free_symbols, assumptions, rng)
        if exact_subs is None:
            return None
        try:
            with _flint_precision(prec_bits):
                ball_subs = {
                    symbol: _sympy_to_acb_ball(value, prec_bits, {})
                    for symbol, value in exact_subs.items()
                }
                val = _sympy_to_acb_ball(expr, prec_bits, ball_subs)
        except UnsupportedNode:
            return None
        if val.is_finite and 0 not in val:
            return ZeroClassification(
                Verdict.NONZERO_LIKELY, "numeric-probe", precision_bits=prec_bits,
                detail=(f"expression is rigorously nonzero at assumption-satisfying witness point "
                        f"{attempt}/{num_points}: {exact_subs}"),
                evidence="assumption-satisfying-witness",
            )
    return None
