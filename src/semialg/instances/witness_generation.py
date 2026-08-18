from __future__ import annotations

import math
import random
from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import sympy as sp

from .modular_sampling import sample_modular_points


@dataclass(frozen=True)
class VariableDomainSpec:
    variable: sp.Symbol
    domains: tuple[object, ...]


def safe_integer_floor_bound(value, *, strict: bool = False):
    value = sp.sympify(value)
    if value in (sp.oo, -sp.oo):
        return value
    if bool(value.is_integer):
        return int(value) - 1 if strict else int(value)
    numeric = sp.N(value, 80)
    floor = math.floor(float(numeric))
    if strict and abs(float(numeric) - floor) < 1e-14:
        floor -= 1
    return floor


def safe_int_ceiling_bound(value, *, strict: bool = False):
    value = sp.sympify(value)
    if value in (sp.oo, -sp.oo):
        return value
    if bool(value.is_integer):
        return int(value) + 1 if strict else int(value)
    numeric = sp.N(value, 80)
    ceiling = math.ceil(float(numeric))
    if strict and abs(ceiling - float(numeric)) < 1e-14:
        ceiling += 1
    return ceiling


def random_sample_from_intv(
    lower, upper, *, strict: bool, sample_count: int, rng: random.Random, integral: bool
) -> list[object]:
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
    lo = float(sp.N(lower, 40))
    hi = float(sp.N(upper, 40))
    if strict:
        eps = max(1e-9, (hi - lo) * 1e-9)
        lo += eps
        hi -= eps
    if hi < lo:
        return []
    if hi == lo:
        return [sp.nsimplify(lo)]
    return list(dict.fromkeys(sp.nsimplify(rng.uniform(lo, hi)) for _ in range(sample_count)))


def sample_free_assignments(
    variables: Sequence[sp.Symbol],
    *,
    domain_rules: Mapping[sp.Symbol, Sequence[object]] | None = None,
    sample_count: int = 5,
    modulus: int | None = None,
    seed: int | None = None,
) -> list[dict[sp.Symbol, object]]:
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
