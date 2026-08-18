import sympy as sp

from semialg.algebraic import RationalSample, compare_samples, sign_at_sample
from semialg.cad import build_collins_proj_set, decomp_collins_complete
from semialg.cad.lifting.sign_invariance import verify_recorded_signs


def test_collins_core_01():
    x, y = sp.symbols("x y")
    tower = build_collins_proj_set([y**2 - x, y - 1], [x, y])
    level1 = {sp.expand(poly.as_expr()) for poly in tower.level(1).polynomials}
    assert x in level1 or x - 1 in level1 or 1 - x in level1
    assert tower.metadata["complete"] is True


def test_collins_core_02():
    x = sp.symbols("x")
    sample = RationalSample(sp.Rational(3, 2))
    assert compare_samples(sample, RationalSample(1)) == 1
    assert sign_at_sample(sp.Poly(x**2 - 2, x), [sample]) == 1


def test_collins_core_03():
    x = sp.symbols("x")
    cad = decomp_collins_complete([x**2 - 1], [x])
    assert cad.complete is True
    assert len(cad.cells) == 5
    check = verify_recorded_signs(cad.cells, cad.tower.level(1).polynomials)
    assert check.ok, check.failures
    signs = [cell.signs[sp.sstr(x**2 - 1)] for cell in cad.cells]
    assert signs == [1, 0, -1, 0, 1]
