from __future__ import annotations

import random
from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import sympy as sp

from ..exact_arithmetic import compare_exact_reals
from .modular_sampling import sample_modular_points


@dataclass(frozen=True)
class VariableDomainSpec:
    variable: sp.Symbol
    domains: tuple[object, ...]


def safe_integer_floor_bound(value, *, strict: bool = False):
    """Greatest integer ``n`` with ``n <= value`` (or ``n < value`` if strict).

    This routine is used in certified witness generation and therefore never
    converts an exact endpoint through binary floating point.
    """
    value = sp.sympify(value)
    if value in (sp.oo, -sp.oo):
        return value
    bound = sp.ceiling(value) - 1 if strict else sp.floor(value)
    if bound.is_Integer:
        return int(bound)
    raise ValueError(f"could not determine exact integer upper bound for {value!s}")


def safe_int_ceiling_bound(value, *, strict: bool = False):
    """Least integer ``n`` with ``n >= value`` (or ``n > value`` if strict)."""
    value = sp.sympify(value)
    if value in (sp.oo, -sp.oo):
        return value
    bound = sp.floor(value) + 1 if strict else sp.ceiling(value)
    if bound.is_Integer:
        return int(bound)
    raise ValueError(f"could not determine exact integer lower bound for {value!s}")


def _dyadic_between(lower: object, upper: object) -> sp.Rational:
    """Return a certified dyadic rational strictly between two exact reals."""

    lower = sp.sympify(lower)
    upper = sp.sympify(upper)
    if compare_exact_reals(lower, upper) >= 0:
        raise ValueError("lower bound must be strictly smaller than upper bound")
    if lower.is_Rational and upper.is_Rational:
        return sp.Rational(lower + upper) / 2

    # Numerical evaluation is used only to propose a dyadic.  The returned
    # candidate is accepted exclusively after exact comparison.
    for bits in (8, 16, 32, 64, 128, 256, 512):
        scale = sp.Integer(2) ** bits
        lo_approx = sp.N(lower, max(30, bits // 3 + 12))
        hi_approx = sp.N(upper, max(30, bits // 3 + 12))
        lo_int = sp.floor(lo_approx * scale)
        hi_int = sp.ceiling(hi_approx * scale)
        start = int(lo_int) - 2
        stop = int(hi_int) + 2
        if stop - start > 64:
            mid = (start + stop) // 2
            candidates = range(mid - 4, mid + 5)
        else:
            candidates = range(start, stop + 1)
        for numerator in candidates:
            candidate = sp.Rational(numerator, scale)
            if (
                compare_exact_reals(lower, candidate) < 0
                and compare_exact_reals(candidate, upper) < 0
            ):
                return candidate
    raise ValueError(f"could not construct a certified rational between {lower} and {upper}")


def _exact_interval_samples(
    lower: object, upper: object, *, strict: bool, sample_count: int, rng: random.Random
) -> list[sp.Expr]:
    lower = sp.sympify(lower)
    upper = sp.sympify(upper)
    comparison = compare_exact_reals(lower, upper)
    if comparison > 0:
        return []
    if comparison == 0:
        return [] if strict else [lower]

    first = _dyadic_between(lower, upper)
    if sample_count <= 1:
        return [first]

    # Recursively split exact intervals.  This yields distinct certified
    # rational samples without relying on binary floating point or epsilons.
    intervals: list[tuple[sp.Expr, sp.Expr]] = [(lower, first), (first, upper)]
    samples: list[sp.Expr] = [first]
    while intervals and len(samples) < sample_count:
        index = rng.randrange(len(intervals))
        lo, hi = intervals.pop(index)
        try:
            candidate = _dyadic_between(lo, hi)
        except ValueError:
            continue
        if candidate not in samples:
            samples.append(candidate)
        intervals.extend(((lo, candidate), (candidate, hi)))
    return samples


def random_sample_from_intv(
    lower,
    upper,
    *,
    strict: bool,
    sample_count: int,
    rng: random.Random,
    integral: bool,
    exact: bool = True,
) -> list[object]:
    """Generate bounded exact or numerical samples from one interval."""
    if integral or lower in (sp.oo, -sp.oo) or upper in (sp.oo, -sp.oo):
        lo = safe_int_ceiling_bound(lower, strict=strict)
        hi = safe_integer_floor_bound(upper, strict=strict)
        if lo in (sp.oo, -sp.oo) or hi in (sp.oo, -sp.oo):
            lo = -50 if lo == -sp.oo else lo
            hi = 50 if hi == sp.oo else hi
        if lo > hi:
            return []
        width = hi - lo + 1
        if width <= sample_count:
            return list(range(lo, hi + 1))
        chosen = set()
        while len(chosen) < sample_count:
            chosen.add(rng.randint(lo, hi))
        return sorted(chosen)
    if exact:
        return _exact_interval_samples(
            lower, upper, strict=strict, sample_count=sample_count, rng=rng
        )

    lo = float(sp.N(lower, 40))
    hi = float(sp.N(upper, 40))
    if strict:
        eps = max(1e-9, (hi - lo) * 1e-9)
        lo += eps
        hi -= eps
    if hi < lo:
        return []
    if hi == lo:
        return [sp.Float(lo, 30)]
    return list(dict.fromkeys(sp.Float(rng.uniform(lo, hi), 30) for _ in range(sample_count)))


def sample_free_assignments(
    variables: Sequence[sp.Symbol],
    *,
    domain_rules: Mapping[sp.Symbol, Sequence[object]] | None = None,
    sample_count: int = 5,
    modulus: int | None = None,
    seed: int | None = None,
) -> list[dict[sp.Symbol, object]]:
    """Generate reproducible sample assignments that respect variable domain restrictions."""
    if sample_count <= 0:
        return []
    variables = tuple(variables)
    if modulus is not None:
        points = sample_modular_points(len(variables), modulus, sample_count, seed=seed)
        return [{var: value for var, value in zip(variables, pt, strict=True)} for pt in points]
    rng = random.Random(seed)
    per_var: list[list[object]] = []
    for var in variables:
        doms = tuple((domain_rules or {}).get(var, ()))
        if any(d == sp.Integers for d in doms) or bool(var.is_integer):
            per_var.append(
                random_sample_from_intv(
                    -50, 50, strict=False, sample_count=sample_count, rng=rng, integral=True
                )
            )
        elif any(d == sp.Reals for d in doms) or bool(var.is_real):
            per_var.append(
                random_sample_from_intv(
                    -20, 20, strict=False, sample_count=sample_count, rng=rng, integral=False
                )
            )
        else:
            reals = random_sample_from_intv(
                -10, 10, strict=False, sample_count=sample_count, rng=rng, integral=False
            )
            imags = random_sample_from_intv(
                -10, 10, strict=False, sample_count=sample_count, rng=rng, integral=False
            )
            per_var.append([r + sp.I * i for r, i in zip(reals, imags, strict=True)])
    assignments = []
    for idx in range(sample_count):
        assn = {}
        for var, pool in zip(variables, per_var, strict=True):
            assn[var] = pool[idx % len(pool)]
        assignments.append(assn)
    return assignments


__all__ = [
    "VariableDomainSpec",
    "safe_integer_floor_bound",
    "safe_int_ceiling_bound",
    "sample_free_assignments",
]
