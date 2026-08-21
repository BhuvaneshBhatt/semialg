from __future__ import annotations

import sympy as sp

from semialg import cad, classify_real_roots, root_of
from semialg.algebraic.roots import isolate_real_roots
from semialg.algebraic.sample_points import choose_sector_sample
from semialg.algebraic.samples import AlgebraicRoot
from semialg.simplify.intervals import Interval1D, merge_intervals


def test_root_classification_string_names_preserve_input_symbol_identity():
    x = sp.Symbol("x")
    a = sp.Symbol("a")
    result = classify_real_roots(x**2 + a, "x", parameters=["a"])
    assert result.variable == x
    assert result.parameters == (a,)


def test_root_of_treats_fiber_variable_as_bound_under_substitution():
    x, y = sp.symbols("x y", real=True)
    expression = root_of(x**2 + y**2 - 1, y, 1)
    assert sp.simplify(expression.subs(x, 1)) == 0
    # Substituting the ambient y coordinate must not rewrite the binder inside
    # an otherwise unspecialized root selector.
    assert expression.subs(y, 17) == expression


def test_exact_root_isolation_fallback_never_uses_nroots(monkeypatch):
    x = sp.Symbol("x", real=True)
    original_real_roots = sp.real_roots

    def fail_real_roots(*args, **kwargs):
        raise NotImplementedError("force RootOf fallback")

    def forbidden_nroots(*args, **kwargs):
        raise AssertionError("nroots must not be used by exact root isolation")

    monkeypatch.setattr(sp, "real_roots", fail_real_roots)
    monkeypatch.setattr(sp, "nroots", forbidden_nroots)
    roots = isolate_real_roots(sp.Poly(x**5 - x - 1, x))
    assert len(roots) == 1
    assert roots[0].as_expr().is_real is True
    # Restore explicitly before comparing the RootOf to avoid depending on the
    # monkeypatch for helper internals.
    monkeypatch.setattr(sp, "real_roots", original_real_roots)
    assert sp.simplify((x**5 - x - 1).subs(x, roots[0].as_expr())) == 0


def test_sector_sample_refines_overlapping_algebraic_intervals_exactly():
    x = sp.Symbol("x", real=True)
    roots = isolate_real_roots(sp.Poly((x**2 - 2) * (x**2 - 3), x))
    left = next(root for root in roots if root.as_expr() == sp.sqrt(2))
    right = next(root for root in roots if root.as_expr() == sp.sqrt(3))
    # Deliberately widen the intervals so they overlap; choose_sector_sample
    # must refine them rather than compare decimal approximations.
    widened_left = AlgebraicRoot(
        left.polynomial, left.interval.expand(1), left.root_index, left.multiplicity, left.root_expr
    )
    widened_right = AlgebraicRoot(
        right.polynomial,
        right.interval.expand(1),
        right.root_index,
        right.multiplicity,
        right.root_expr,
    )
    sample = choose_sector_sample(widened_left, widened_right)
    assert sample.value > sp.sqrt(2)
    assert sample.value < sp.sqrt(3)


def test_interval_simplification_orders_extremely_close_algebraics_exactly():
    a = sp.sqrt(2)
    b = a + sp.Rational(1, 10**120)
    result = merge_intervals(
        (
            Interval1D(b, b + 1, True, True),
            Interval1D(a, a + sp.Rational(1, 10**121), True, True),
        )
    )
    assert result[0].left == a
    assert result[1].left == b


def test_disk_boundary_excludes_interior_and_root_of_substitutes_at_endpoints():
    x, y = sp.symbols("x y", real=True)
    boundary = cad(x**2 + y**2 <= 1, [x, y], operation="boundary")
    assert bool(boundary.formula.subs({x: 1, y: 0}))
    assert bool(boundary.formula.subs({x: 0, y: 1}))
    assert not bool(boundary.formula.subs({x: 0, y: 0}))
