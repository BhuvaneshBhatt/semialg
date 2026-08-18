from __future__ import annotations

import sympy as sp

from semialg import GenericSplit, generic_cad


def test_generic_cad_nonparametric_disk_splits_interior_and_boundary():
    x, y = sp.symbols("x y", real=True)
    result = generic_cad(x**2 + y**2 <= 1, [x, y])

    assert result.status == "complete"
    assert isinstance(result.generic_split, GenericSplit)
    assert result.generic_split.generic_cells
    assert result.generic_split.exceptional_cells
    assert bool(result.generic_formula.subs({x: 0, y: 0}))
    assert not bool(result.generic_formula.subs({x: 1, y: 0}))
    assert bool(result.exceptional_formula.subs({x: 1, y: 0}))
    assert any(
        sp.expand(poly.as_expr() - (x**2 + y**2 - 1)) == 0
        for poly in result.generic_split.exceptional_polys
    )


def test_generic_cad_nonparametric_strict_region_keeps_deleted_boundary_exceptional():
    x = sp.symbols("x", real=True)
    result = generic_cad(x**2 < 1, [x])

    assert bool(result.generic_formula.subs({x: 0}))
    assert not bool(result.generic_formula.subs({x: 1}))
    assert bool(result.exceptional_formula.subs({x: 1}))
    assert bool(result.exceptional_formula.subs({x: -1}))


def test_parameterized_generic_cad_keeps_parameter_exceptional_cases():
    a, x = sp.symbols("a x", real=True)
    result = generic_cad(sp.Eq(a * x - 1, 0), variables=[x], parameters=[a])

    assert result.generic_split is None
    assert result.generic_cases
    assert result.exceptional_cases
    assert bool(result.generic_formula.subs({a: 2}))
    assert bool(result.exceptional_formula.subs({a: 0}))
    assert result.diagnostics["exceptional_polynomial_count"] >= 1
