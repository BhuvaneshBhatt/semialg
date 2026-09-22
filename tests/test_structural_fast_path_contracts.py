"""Structural performance contracts guard algorithm selection without timing noise."""

from __future__ import annotations

import inspect

import sympy as sp

from semialg.algebraic import certify_polynomial_root_interval
from semialg.algebraic.sample_points import choose_sector_sample
from semialg.algebraic.samples import RationalSample


def test_descartes_decisive_intervals_report_the_fast_exact_method():
    x = sp.Symbol("x")
    unique = certify_polynomial_root_interval(x**2 - 2, 0, 2, var=x)
    empty = certify_polynomial_root_interval(x**2 - 2, 2, 3, var=x)

    assert unique.root_count == 1 and unique.method == "descartes"
    assert empty.root_count == 0 and empty.method == "descartes"


def test_ambiguous_descartes_interval_reports_exact_sturm_fallback():
    x = sp.Symbol("x")
    certificate = certify_polynomial_root_interval(
        (x + 3) * (x + 1) * (x - 1) * (x - 2), -4, 4, var=x
    )

    assert certificate.descartes_variations is not None
    assert certificate.descartes_variations > 1
    assert certificate.method == "descartes+sturm"
    assert certificate.root_count == 4


def test_sector_sampling_uses_shared_separator():
    source = inspect.getsource(choose_sector_sample)
    assert "rational_between_algebraic_reals" in source
    assert "range(32)" not in source

    assert choose_sector_sample(RationalSample(-1), RationalSample(1)).as_expr() == 0
