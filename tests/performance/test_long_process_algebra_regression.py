"""Regression for cross-subsystem algebraic process degradation."""

from __future__ import annotations

import gc
import resource
from time import perf_counter

import pytest
import sympy as sp

import semialg
from semialg.algebraic_decomposition import primary_decomposition
from semialg.algebraic_function_fields import (
    MonogenicFunctionField,
    RationalFunctionField,
    compress_primitive_element,
    verify_primitive_element_compression,
)
from semialg.cache_control import clear_caches


def _rss_mb() -> float:
    # Linux reports KiB; macOS reports bytes. The CI/performance environment is Linux.
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0


def _primitive_workload():
    x, a, b = sp.symbols("x a b")
    base = RationalFunctionField((x,))
    first = MonogenicFunctionField(base, a, (-x, 0, 1))
    second = MonogenicFunctionField(first, b, (-first.alpha, first.zero, first.one))
    result = compress_primitive_element(second)
    assert verify_primitive_element_compression(result)

    c = sp.Symbol("c")
    # Degree-one third layer exercises arbitrary-depth transport without
    # changing the total finite extension degree.
    third = MonogenicFunctionField(second, c, (-second.alpha, second.one))
    deep = compress_primitive_element(third)
    assert verify_primitive_element_compression(deep)


@pytest.mark.performance
@pytest.mark.slow
def test_gtz_then_primitive_compression_has_bounded_process_degradation():
    """A GTZ-heavy prefix must not make later primitive compression pathological."""
    x, y = sp.symbols("x y")
    clear_caches(include_sympy=True, collect=False)
    before_objects = len(gc.get_objects())
    before_rss = _rss_mb()

    # Populate the same exact-algebra families implicated by the original
    # long-process reproducer without depending on test execution order.
    ideals = (
        ((x + y) ** 2, y * (x + y)),
        (x**2, x * y),
        (x * y,),
        ((x - y) ** 2 * (x + y),),
    )
    for generators in ideals:
        result = primary_decomposition(generators, (x, y))
        assert result.complete
        assert semialg.replay_certificate(result).verified

    start = perf_counter()
    _primitive_workload()
    elapsed = perf_counter() - start
    after_objects = len(gc.get_objects())
    rss_growth = _rss_mb() - before_rss

    # These generous tripwires catch the known pathological long-process behavior.
    assert elapsed < 30.0
    assert rss_growth < 512.0
    assert after_objects - before_objects < 1_000_000
