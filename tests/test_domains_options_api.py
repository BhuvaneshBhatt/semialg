import sympy as sp

from semialg.domains import SemialgOptions, SolveDomain, apply_assumptions, normalize_domain


def test_domains_options_01():
    assert normalize_domain("reals") is SolveDomain.REALS
    assert normalize_domain("Q") is SolveDomain.RATIONALS
    assert normalize_domain("Z") is SolveDomain.INTEGERS
    assert normalize_domain("booleans") is SolveDomain.BOOLEANS


def test_domains_options_02():
    x = sp.Symbol("x", real=True)
    opts = SemialgOptions.from_values(domain="reals", assumptions=[x > 0], count=3)
    assert opts.domain is SolveDomain.REALS
    assert opts.count == 3
    assert opts.assumptions == (x > 0,)


def test_domains_options_03():
    x = sp.Symbol("x", real=True)
    expr = apply_assumptions(x < 2, [x > 0])
    assert expr == sp.And(x > 0, x < 2)
