import sympy as sp

from semialg.solve.equality_ideal import EqualityIdealContext


def test_exact_dimension_and_quotient_dimension():
    x, y = sp.symbols("x y", real=True)
    zero = EqualityIdealContext.build([x**2 - 1, y - x], (x, y))
    assert zero.dimension == 0
    assert zero.quotient_dimension == 2
    positive = EqualityIdealContext.build([x * y], (x, y))
    assert positive.dimension == 1
    assert positive.quotient_dimension is None


def test_fglm_coordinate_polynomials():
    x, y = sp.symbols("x y", real=True)
    ctx = EqualityIdealContext.build([x**2 - 2, y - x], (x, y))
    px, py = ctx.fglm_coordinate_polynomials()
    assert sp.Poly(px, x).monic().as_expr() == x**2 - 2
    assert sp.Poly(py, y).monic().as_expr() == y**2 - 2


def test_gbconvert_style_linear_elimination():
    x, y = sp.symbols("x y", real=True)
    ctx = EqualityIdealContext.build([x - y**2, y**2 - 2], (x, y))
    elim = ctx.eliminate_linear_variable(x)
    assert elim is not None
    assert ctx.normal_form(x - elim.replacement) == 0
    assert x not in sp.Tuple(*elim.remaining_basis).free_symbols


def test_radical_simplification_collapses_inequalities():
    x = sp.symbols("x", real=True)
    ctx = EqualityIdealContext.build([x**2], (x,))
    assert ctx.radical_contains(x)
    simplified, meta = ctx.simplify_relations(sp.And(sp.Eq(x**2, 0), x >= 0))
    assert simplified == sp.Eq(x**2, 0)
    assert meta["radical_collapses"] == 1
    impossible, _ = ctx.simplify_relations(sp.And(sp.Eq(x**2, 0), x > 0))
    assert impossible is sp.false


def test_zero_dimensional_inequality_filtering():
    x, y = sp.symbols("x y", real=True)
    ctx = EqualityIdealContext.build([x**2 - 1, y - x], (x, y))
    result = ctx.filter_zero_dimensional(x > 0)
    assert result.points == ((sp.Integer(1), sp.Integer(1)),)
    assert result.quotient_dimension == 2
    assert len(result.coordinate_polynomials) == 2


def test_unit_ideal_and_nonvanishing_analysis():
    x, y = sp.symbols("x y")
    unit = EqualityIdealContext.build([sp.Integer(1)], (x, y))
    assert unit.unit_ideal
    assert unit.dimension == -1
    assert unit.quotient_dimension == 0

    ctx = EqualityIdealContext.build([x**2 - 2], (x,))
    assert ctx.has_common_zero_with(x) is False
    assert ctx.simplify_constraints(sp.Ne(x, 0)) is sp.true


def test_all_linear_eliminations_are_certified():
    x, y = sp.symbols("x y")
    ctx = EqualityIdealContext.build((x - y, y**2 - 2), (x, y))
    eliminations = ctx.linear_eliminations()
    assert eliminations
    assert all(ctx.normal_form(item.certificate) == 0 for item in eliminations)
