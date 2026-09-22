import pytest
import sympy as sp

from semialg import Ball, Box, Polygon, Simplex, Sphere
from semialg.region_integrate import integrate_over_region
from semialg.standard_region_integrate import integrate_over_standard_region


def _formula_integral(region, expression, variables, *, dimension):
    formula = region.as_formula(variables, eliminate=True)
    return sp.sympify(
        integrate_over_region(
            expression,
            formula,
            variables,
            measure_dimension=dimension,
            method="symbolic",
        )
    )


@pytest.mark.parametrize(
    ("region", "variables", "expression"),
    [
        (
            Box(((0, 2), (-1, 1))),
            sp.symbols("x y", real=True),
            sp.Symbol("x", real=True) ** 2,
        ),
        (
            Simplex(((0, 0), (2, 0), (0, 2))),
            sp.symbols("x y", real=True),
            sp.Symbol("x", real=True) + sp.Symbol("y", real=True),
        ),
        (
            Ball((0, 0), 2),
            sp.symbols("x y", real=True),
            sp.Symbol("x", real=True) ** 2 + sp.Symbol("y", real=True) ** 2,
        ),
        (
            Polygon(((0, 0), (2, 0), (2, 1), (0, 1))),
            sp.symbols("x y", real=True),
            sp.Symbol("x", real=True),
        ),
    ],
    ids=["box", "simplex", "ball", "polygon"],
)
def test_canonical_formula_integrals_agree(region, variables, expression):
    canonical = sp.sympify(integrate_over_standard_region(expression, region, variables))
    formula = _formula_integral(region, expression, variables, dimension=region.dimension())
    assert sp.simplify(canonical - formula) == 0


def test_canonical_and_cad_integrals_agree():
    x, y = sp.symbols("x y", real=True)
    region = Polygon(((0, 0), (2, 0), (2, 1), (0, 1)))
    canonical = integrate_over_standard_region(x + y, region, (x, y))
    result = integrate_over_region(
        x + y,
        region.as_formula((x, y), eliminate=True),
        (x, y),
        return_result=True,
    )
    assert "cad" in result.method
    assert sp.simplify(result.value - canonical) == 0


@pytest.mark.parametrize(
    ("region", "variables"),
    [
        (Simplex(((0, 0), (3, 4))), sp.symbols("x y", real=True)),
        (Sphere((0, 0), 2), sp.symbols("x y", real=True)),
    ],
    ids=["segment", "circle"],
)
def test_canonical_formula_intrinsic_measures_agree(region, variables):
    formula = region.as_formula(variables, eliminate=True)
    from_formula = sp.sympify(
        integrate_over_region(
            1,
            formula,
            variables,
            measure_dimension=region.dimension(),
            method="symbolic",
        )
    )
    assert sp.simplify(region.measure() - from_formula) == 0


def test_intrinsic_chart_and_formula_agree_for_line_graph():
    x, y = sp.symbols("x y", real=True)
    formula = sp.And(sp.Eq(y, x), x >= 0, x <= 1)
    from semialg.parametric_geometry import intrinsic_parametric_cover

    cover = intrinsic_parametric_cover(formula, (x, y), dimension=1)
    assert len(cover.charts) == 1
    chart = cover.charts[0]
    chart_value = sp.integrate(chart.metric_factor(), *chart.bounds)
    formula_value = integrate_over_region(
        1, formula, (x, y), measure_dimension=1, method="symbolic"
    )
    assert sp.simplify(chart_value - formula_value) == 0
    assert sp.simplify(chart_value - sp.sqrt(2)) == 0


@pytest.mark.slow
def test_embedded_simplex_cad_measure_agrees():
    x, y, z = sp.symbols("x y z", real=True)
    region = Simplex(((0, 0, 1), (1, 0, 1), (0, 1, 1)))
    formula = region.as_formula((x, y, z), eliminate=True)
    result = integrate_over_region(
        1,
        formula,
        (x, y, z),
        measure_dimension=2,
        return_result=True,
    )
    assert "intrinsic" in result.method or "cad" in result.method
    assert sp.simplify(result.value - region.measure()) == 0


def test_embedded_simplex_intrinsic_measure_precedes_ambient_cad():
    x, y, z = sp.symbols("x y z", real=True)
    region = Simplex(((0, 0, 1), (1, 0, 1), (0, 1, 1)))
    formula = region.as_formula((x, y, z), eliminate=True)

    result = integrate_over_region(1, formula, (x, y, z), measure_dimension=2, return_result=True)
    assert result.method == "affine_intrinsic_parametric_measure"
    assert sp.simplify(result.value - sp.Rational(1, 2)) == 0
