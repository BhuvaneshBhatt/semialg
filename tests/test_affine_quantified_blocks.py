"""Safe existential substitutions precede expensive graph elimination."""

from unittest.mock import patch

import sympy as sp

from semialg import SemialgebraicRegion
from semialg.quantifiers import Exists, ForAll


def test_affine_block_avoids_vs():
    x, t = sp.symbols("x t", real=True)
    region = SemialgebraicRegion(Exists(t, sp.And(sp.Eq(2 * t, x + 1), t >= 0)), (x,))
    with patch(
        "semialg.solve.reduce.try_quadratic_vs_qe", side_effect=AssertionError("unneeded VS")
    ):
        result = region.quantifier_free_formula()
    for value in (-2, -1, 0, 1):
        assert bool(result.subs(x, value)) == (value >= -1)


def test_universal_equation_is_not_substitution():
    x, t = sp.symbols("x t", real=True)
    region = SemialgebraicRegion(ForAll(t, sp.Eq(t, x)), (x,))
    assert region.quantifier_free_formula() is sp.false


def test_zero_coefficient_stratum():
    a, x, t = sp.symbols("a x t", real=True)
    region = SemialgebraicRegion(Exists(t, sp.Eq(a * t, x)), (a, x))
    result = region.quantifier_free_formula()
    for av in (-1, 0, 1):
        for xv in (-1, 0, 1):
            assert bool(result.subs({a: av, x: xv})) == (av != 0 or xv == 0)


def test_mixed_prefix():
    u, v = sp.symbols("u v", real=True)
    region = SemialgebraicRegion(Exists(v, ForAll(u, sp.Eq(u, v))), ())
    assert region.quantifier_free_formula() is sp.false
