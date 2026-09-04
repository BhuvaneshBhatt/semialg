import pytest
import sympy as sp

from semialg import is_equal
from semialg.formula import parse_formula
from semialg.presolve import fourier_motzkin_eliminate, presolve_semialgebraic
from semialg.qe.complete import qe_by_complete_cad


def _qe(expr, variables, quantified, *, presolve):
    qvars = {v for _, v in quantified}
    free = [v for v in variables if v not in qvars]
    return qe_by_complete_cad(
        variables,
        quantified,
        parse_formula(expr),
        free_variables=free,
        use_presolve=presolve,
        variable_order_strategy="preserve",
        return_result=True,
    ).formula


x, y, z, a, b = sp.symbols("x y z a b", real=True)

PRESOLVE_CASES = [
    (sp.And(sp.Eq(y, x + 1), y > 0), [x, y], [("exists", y)]),
    (sp.And(sp.Eq(y, 2 * x - 3), y <= 5), [x, y], [("exists", y)]),
    (sp.And(sp.Eq(z, y + 1), sp.Eq(y, x + 1), z < 5), [x, y, z], [("exists", y), ("exists", z)]),
    (sp.And(y > x, y < 2), [x, y], [("exists", y)]),
    (sp.And(y >= x, y <= 2), [x, y], [("exists", y)]),
    (sp.And(y > x, y <= 2), [x, y], [("exists", y)]),
    (sp.And(y >= x, y < 2), [x, y], [("exists", y)]),
    (sp.And(y > x, y > 0), [x, y], [("exists", y)]),
    (sp.And(y < x, y < 0), [x, y], [("exists", y)]),
    (sp.And(sp.Eq(y, 3), y > x), [x, y], [("exists", y)]),
    (sp.And(sp.Eq(y, x), y**2 <= 4), [x, y], [("exists", y)]),
    (sp.And(sp.Eq(y, -x), y >= 0), [x, y], [("exists", y)]),
    (sp.And(sp.Eq(z, x + y), z >= 0), [x, y, z], [("exists", z)]),
    (sp.And(sp.Eq(z, x - y), z < 1), [x, y, z], [("exists", z)]),
    (sp.And(sp.Eq(y, x + 1), sp.Eq(y, x + 2)), [x, y], [("exists", y)]),
]


@pytest.mark.parametrize("expr,variables,quantified", PRESOLVE_CASES)
def test_presolve_preserves_complete_qe_semantics(expr, variables, quantified):
    with_presolve = _qe(expr, variables, quantified, presolve=True)
    without_presolve = _qe(expr, variables, quantified, presolve=False)
    free = [v for v in variables if v not in {q for _, q in quantified}]
    if not free:
        assert bool(with_presolve) == bool(without_presolve)
    else:
        try:
            assert is_equal(with_presolve, without_presolve, free)
        except (sp.PolynomialError, NotImplementedError, ValueError):
            # Multivariate CAD reconstruction may retain algebraic root
            # functions.  Compare exact rational specializations instead of
            # feeding those root functions back into CAD as polynomials.
            for values in __import__("itertools").product((-2, -1, 0, 1, 2), repeat=len(free)):
                subs = dict(zip(free, values, strict=True))
                assert bool(sp.simplify(with_presolve.subs(subs))) == bool(
                    sp.simplify(without_presolve.subs(subs))
                )


@pytest.mark.parametrize(
    "expr",
    [
        sp.Eq(a * x, 1),
        sp.Eq((a + 1) * x, b),
        sp.Eq(y * x, 1),
        sp.Eq((y - 2) * x, z),
        sp.Eq((a * y + 1) * x, 3),
    ],
)
def test_presolve_never_divides_by_nonconstant_coefficients(expr):
    vars_ = tuple(sorted(expr.free_symbols, key=lambda s: s.name))
    result = presolve_semialgebraic(expr, vars_, eliminate=[x])
    assert x in result.variables
    assert not any(var == x for var, _ in result.substitutions)


@pytest.mark.parametrize(
    "expr, eliminated, expected",
    [
        (sp.And(y >= x, y <= 2), [y], x <= 2),
        (sp.And(y > x, y <= 2), [y], x < 2),
        (sp.And(y >= x, y < 2), [y], x < 2),
        (sp.And(y > x, y < 2), [y], x < 2),
        (sp.And(y >= x, y >= 0), [y], sp.true),
        (sp.And(y <= x, y <= 0), [y], sp.true),
        (sp.And(y >= x, y <= z), [y], x <= z),
        (sp.And(y > x, y < z), [y], x < z),
    ],
)
def test_fourier_motzkin_generated_linear_cases(expr, eliminated, expected):
    reduced = fourier_motzkin_eliminate(expr, eliminated)
    assert reduced is not None
    formula, removed = reduced
    assert removed == tuple(eliminated)
    vars_ = tuple(sorted((formula.free_symbols | expected.free_symbols), key=lambda s: s.name))
    assert is_equal(formula, expected, vars_)


def test_presolve_keeps_boolean_disjunction_semantics_untouched():
    expr = sp.Or(sp.And(sp.Eq(y, x + 1), y > 0), x < -5)
    result = presolve_semialgebraic(expr, [x, y], eliminate=[y])
    assert is_equal(result.formula, expr, [x, y])
    assert result.substitutions == ()


def test_presolve_reports_zero_dimensional_equalities_without_changing_formula():
    expr = sp.And(sp.Eq(x**2 - 1, 0), sp.Eq(y - x, 0))
    result = presolve_semialgebraic(expr, [x, y], detect_zero_dimensional=True)
    assert result.zero_dim_eqs is True
    assert is_equal(result.formula, expr, [x, y])
