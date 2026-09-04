import random

import sympy as sp

from semialg.exact_arithmetic import compare_exact_reals
from semialg.instances.real_fallbacks import rational_bound
from semialg.instances.witness_generation import random_sample_from_intv


def test_exact_sampling_preserves_sub_float_width_interval():
    lo = sp.Integer(1)
    hi = lo + sp.Rational(1, 10**30)
    samples = random_sample_from_intv(
        lo, hi, strict=True, sample_count=4, rng=random.Random(7), integral=False, exact=True
    )
    assert samples
    assert all(value.is_Rational for value in samples)
    assert all(compare_exact_reals(lo, value) < 0 for value in samples)
    assert all(compare_exact_reals(value, hi) < 0 for value in samples)


def test_exact_sampling_handles_algebraic_endpoints():
    lo = sp.sqrt(2)
    hi = sp.sqrt(2) + sp.Rational(1, 10**20)
    samples = random_sample_from_intv(
        lo, hi, strict=True, sample_count=3, rng=random.Random(11), integral=False, exact=True
    )
    assert samples
    assert all(compare_exact_reals(lo, value) < 0 for value in samples)
    assert all(compare_exact_reals(value, hi) < 0 for value in samples)


def test_rational_bound_is_certified_exactly():
    value = sp.sqrt(2)
    lower = rational_bound(value, -1)
    upper = rational_bound(value, 1)
    assert lower.is_Rational and upper.is_Rational
    assert compare_exact_reals(lower, value) <= 0
    assert compare_exact_reals(value, upper) <= 0
