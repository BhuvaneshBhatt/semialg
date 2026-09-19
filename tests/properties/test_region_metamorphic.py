from __future__ import annotations

import pytest
import sympy as sp

hypothesis = pytest.importorskip("hypothesis")
from hypothesis import given, settings
from hypothesis import strategies as st

from semialg import SemialgebraicRegion

pytestmark = pytest.mark.slow

small = st.integers(min_value=-3, max_value=3)


@settings(max_examples=12, deadline=None)
@given(a=small, b=small, c=small)
def test_interval_boolean_operations_are_commutative(a, b, c):
    lo1, hi1 = sorted((a, b))
    lo2, hi2 = sorted((b, c))
    x = sp.symbols("x", real=True)
    r1 = SemialgebraicRegion(sp.And(x >= lo1, x <= hi1), (x,))
    r2 = SemialgebraicRegion(sp.And(x >= lo2, x <= hi2), (x,))

    assert r1.union(r2).equals_region(r2.union(r1))
    assert r1.intersection(r2).equals_region(r2.intersection(r1))
    assert r1.difference(r2).equals_region(r1.intersection(r2.complement()))


@settings(max_examples=10, deadline=None)
@given(a=small, b=small)
def test_interval_euler_inclusion_exclusion(a, b):
    lo, hi = sorted((a, b))
    x = sp.symbols("x", real=True)
    a_region = SemialgebraicRegion(sp.And(x >= lo, x <= hi), (x,))
    b_region = SemialgebraicRegion(sp.And(x >= 0, x <= 2), (x,))
    union = a_region.union(b_region)
    inter = a_region.intersection(b_region)

    assert union.euler_characteristic() == (
        a_region.euler_characteristic()
        + b_region.euler_characteristic()
        - inter.euler_characteristic()
    )


@settings(max_examples=10, deadline=None)
@given(scale=st.integers(min_value=1, max_value=4), shift=small)
def test_affine_image_preimage_roundtrip(scale, shift):
    x, y = sp.symbols("x y", real=True)
    source = SemialgebraicRegion(sp.And(x >= -1, x <= 2), (x,))
    image = source.image(scale * x + shift, variables=(y,))
    assert image.preimage(scale * x + shift, (x,)).equals_region(source)
