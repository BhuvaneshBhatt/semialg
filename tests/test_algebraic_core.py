import pytest
import sympy as sp

from semialg.algebraic import (
    AlgebraicRoot,
    RationalInterval,
    RationalSample,
    compare_samples,
    isolate_real_roots,
    refine_isol_intv,
    sign_at_sample,
)
from semialg.cad import decomp_collins_complete


def test_algebraic_core_01():
    with pytest.raises(TypeError):
        compare_samples(RationalSample(1), sp.Integer(2))


def test_algebraic_core_02():
    x = sp.symbols("x")
    roots = isolate_real_roots(sp.Poly(x**2 - 2, x))
    assert len(roots) == 2
    assert all(isinstance(root, AlgebraicRoot) for root in roots)
    assert compare_samples(roots[0], RationalSample(0)) < 0
    assert compare_samples(roots[1], RationalSample(0)) > 0
    assert roots[0].multiplicity == 1


def test_algebraic_core_03():
    x = sp.symbols("x")
    roots = isolate_real_roots(sp.Poly((x - 1) ** 2, x))
    assert len(roots) == 1
    assert roots[0].multiplicity == 2
    assert roots[0].interval == RationalInterval(1, 1)


def test_algebraic_core_04():
    x = sp.symbols("x")
    roots = isolate_real_roots(sp.Poly(x**2 - 2, x))
    assert sign_at_sample(sp.Poly(x**2 - 2, x), [roots[0]]) == 0
    assert sign_at_sample(sp.Poly(x, x), [roots[0]]) == -1
    assert sign_at_sample(sp.Poly(x, x), [roots[1]]) == 1


def test_algebraic_core_05():
    x = sp.symbols("x")
    root = isolate_real_roots(sp.Poly(x**2 - 2, x))[1]
    refined = refine_isol_intv(root, steps=3)
    assert refined.polynomial == root.polynomial
    assert refined.root_index == root.root_index
    assert refined.interval.width <= root.interval.width


def test_algebraic_core_06():
    x = sp.symbols("x")
    cad = decomp_collins_complete([x**2 - 1], [x])
    assert len(cad.cells) == 5
    assert all(isinstance(cell.sample[0], (RationalSample, AlgebraicRoot)) for cell in cad.cells)
    assert [cell.signs[sp.sstr(x**2 - 1)] for cell in cad.cells] == [1, 0, -1, 0, 1]
