"""Expansion-free finite-field identity testing with explicit error bounds."""

from __future__ import annotations

import math
import random
from dataclasses import dataclass

try:
    import flint
except ImportError:
    flint = None
import sympy as sp

from . import _config as cfg
from ._result import Verdict, ZeroClassification


class NotPolynomialRing(Exception):
    """Raised when an expression is not a rational function over Q."""


@dataclass(frozen=True)
class RationalBounds:
    """Conservative bounds for N/D without constructing N or D."""

    num_degree: int
    den_degree: int
    num_log2_l1: int
    den_log2_l1: int


def _log2_l1_integer(value: int) -> int:
    value = abs(int(value))
    return 0 if value == 0 else max(1, value.bit_length())


def _logadd2(a: int, b: int) -> int:
    if a == 0:
        return b
    if b == 0:
        return a
    return max(a, b) + 1


def _mul_bounds(a: RationalBounds, b: RationalBounds) -> RationalBounds:
    return RationalBounds(
        a.num_degree + b.num_degree,
        a.den_degree + b.den_degree,
        a.num_log2_l1 + b.num_log2_l1,
        a.den_log2_l1 + b.den_log2_l1,
    )


def _add_bounds(a: RationalBounds, b: RationalBounds) -> RationalBounds:
    # a/b + c/d = (ad + cb)/(bd)
    left = a.num_log2_l1 + b.den_log2_l1
    right = b.num_log2_l1 + a.den_log2_l1
    return RationalBounds(
        max(a.num_degree + b.den_degree, b.num_degree + a.den_degree),
        a.den_degree + b.den_degree,
        _logadd2(left, right),
        a.den_log2_l1 + b.den_log2_l1,
    )


def rational_function_bounds(expr: sp.Expr, symbols) -> RationalBounds:
    """Bound rational numerator/denominator directly on the expression DAG."""
    symbol_set = set(symbols)
    if expr in symbol_set:
        return RationalBounds(1, 0, 1, 1)
    if expr.is_Integer:
        return RationalBounds(0, 0, _log2_l1_integer(int(expr)), 1)
    if expr.is_Rational:
        return RationalBounds(
            0,
            0,
            _log2_l1_integer(int(expr.p)),
            _log2_l1_integer(int(expr.q)),
        )
    if expr.is_Symbol:
        raise NotPolynomialRing(f"unrecognized symbol {expr}")
    if expr.is_Add:
        result = RationalBounds(0, 0, 0, 1)
        for arg in expr.args:
            result = _add_bounds(result, rational_function_bounds(arg, symbols))
        return result
    if expr.is_Mul:
        result = RationalBounds(0, 0, 1, 1)
        for arg in expr.args:
            result = _mul_bounds(result, rational_function_bounds(arg, symbols))
        return result
    if expr.is_Pow:
        base, exponent = expr.args
        if not exponent.is_Integer:
            raise NotPolynomialRing(f"non-integer exponent {exponent}")
        power = int(exponent)
        bounds = rational_function_bounds(base, symbols)
        if power == 0:
            return RationalBounds(0, 0, 1, 1)
        if power < 0:
            bounds = RationalBounds(
                bounds.den_degree,
                bounds.num_degree,
                bounds.den_log2_l1,
                bounds.num_log2_l1,
            )
            power = -power
        return RationalBounds(
            bounds.num_degree * power,
            bounds.den_degree * power,
            bounds.num_log2_l1 * power,
            bounds.den_log2_l1 * power,
        )
    if not (expr.free_symbols & symbol_set) and expr.is_Rational:
        return rational_function_bounds(sp.Rational(expr), symbols)
    raise NotPolynomialRing(f"unsupported rational-function node {type(expr).__name__}")


def _evaluate_tree_modulo_prime(expr: sp.Expr, point: dict, ctx: flint.fmpz_mod_ctx, p: int):
    if expr in point:
        return ctx(point[expr])
    if expr.is_Integer:
        return ctx(int(expr) % p)
    if expr.is_Rational:
        num, den = int(expr.p), int(expr.q)
        if den % p == 0:
            raise NotPolynomialRing("rational denominator divisible by p")
        return ctx(num % p) * (ctx(den % p) ** -1)
    if expr.is_Add:
        total = ctx(0)
        for arg in expr.args:
            total += _evaluate_tree_modulo_prime(arg, point, ctx, p)
        return total
    if expr.is_Mul:
        total = ctx(1)
        for arg in expr.args:
            total *= _evaluate_tree_modulo_prime(arg, point, ctx, p)
        return total
    if expr.is_Pow:
        base, exponent = expr.args
        if not exponent.is_Integer:
            raise NotPolynomialRing(f"non-integer exponent {exponent}")
        value = _evaluate_tree_modulo_prime(base, point, ctx, p)
        power = int(exponent)
        if power < 0 and int(value) % p == 0:
            raise ZeroDivisionError("sample lies on a rational-function pole modulo p")
        return value ** power
    if expr.is_Symbol:
        raise NotPolynomialRing(f"unassigned symbol {expr}")
    raise NotPolynomialRing(f"unsupported node {type(expr).__name__}")


def _prime_for_bits(bits: int, rng: random.Random) -> int:
    candidate = rng.getrandbits(bits) | (1 << (bits - 1)) | 1
    return int(sp.nextprime(candidate))


def _plan_field(degree_num: int, degree_den: int, coeff_bits: int, target_error: float):
    if not (0.0 < target_error < 1.0):
        raise ValueError("target_error must lie strictly between 0 and 1")
    best = None
    for trials in range(1, cfg.SZ_MAX_TRIALS + 1):
        allowed = target_error ** (1.0 / trials)
        required_p = degree_den + degree_num / allowed
        error_bits = max(2, math.ceil(math.log2(max(required_p, 2))) + 1)
        bits = max(cfg.SZ_FIELD_MIN_BITS, coeff_bits + 1, error_bits)
        candidate = (bits * trials, bits, trials)
        if best is None or candidate < best:
            best = candidate
    _, bits, trials = best
    return bits, trials


def finite_field_identity_test(
    expr: sp.Expr,
    symbols,
    target_error: float = cfg.SZ_TARGET_FALSE_POSITIVE,
    rng: random.Random | None = None,
) -> ZeroClassification:
    """Test a rational-function identity without symbolic common denominators.

    The original DAG is evaluated directly.  Separate recursive bounds track
    the degree and coefficient size of the implicit numerator/denominator,
    which is enough to choose good primes and a Schwartz-Zippel error budget.
    """
    if flint is None:
        return ZeroClassification(Verdict.UNKNOWN, "finite-field", detail="python-flint is unavailable")

    rng = rng or random.Random()
    try:
        bounds = rational_function_bounds(expr, symbols)
        degree_num = max(bounds.num_degree, 1)
        degree_den = bounds.den_degree
        field_bits, trials_needed = _plan_field(
            degree_num,
            degree_den,
            max(bounds.num_log2_l1, bounds.den_log2_l1, 1),
            target_error,
        )
    except (NotPolynomialRing, ValueError) as exc:
        return ZeroClassification(Verdict.UNKNOWN, "finite-field", detail=str(exc))

    trials = 0
    compounded = 1.0
    for _ in range(trials_needed):
        p = _prime_for_bits(field_bits, rng)
        ctx = flint.fmpz_mod_ctx(p)
        if p <= degree_den:
            return ZeroClassification(Verdict.UNKNOWN, "finite-field", detail="chosen field is too small")

        value = None
        point = None
        for _attempt in range(cfg.SZ_DENOM_RETRIES):
            point = {symbol: rng.randrange(0, p) for symbol in symbols}
            try:
                value = _evaluate_tree_modulo_prime(expr, point, ctx, p)
                break
            except ZeroDivisionError:
                continue
            except NotPolynomialRing as exc:
                return ZeroClassification(Verdict.UNKNOWN, "finite-field", detail=str(exc))
        if value is None:
            return ZeroClassification(
                Verdict.UNKNOWN,
                "finite-field",
                detail="could not sample away from rational-function poles",
            )

        trials += 1
        if int(value) % p != 0:
            return ZeroClassification(
                Verdict.NONZERO_PROVEN,
                "finite-field",
                trials=trials,
                detail=(f"expression is nonzero modulo a {p.bit_length()}-bit good prime at {point}; "
                        "the prime exceeds the implicit numerator coefficient bound"),
                evidence="certified-good-prime-witness",
            )
        per_trial = min(1.0, degree_num / (p - degree_den))
        compounded *= per_trial

    return ZeroClassification(
        Verdict.ZERO_UNPROVEN,
        "finite-field",
        trials=trials,
        detail=(f"expression vanished in {trials} independent pole-free samples over "
                f"~{field_bits}-bit prime fields; implicit numerator degree bound {degree_num}, "
                f"denominator degree bound {degree_den}"),
        evidence="schwartz-zippel",
        error_bound=compounded,
        requested_error=target_error,
    )
