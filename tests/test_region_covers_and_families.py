"""Structural optimization geometry and exact parameter-fiber transitions."""

import pytest
import sympy as sp

from semialg import (
    Ball,
    BooleanRegion,
    Box,
    Ellipsoid,
    ParametricRegion,
    Simplex,
    Sphere,
    TransformedRegion,
)
from semialg.region_covers import bounding_cover, optimization_cover
from semialg.region_families import region_family

x, y, t, a = sp.symbols("x y t a", real=True)


@pytest.mark.parametrize(
    "region,expected",
    [
        (Box(((-1, 2), (3, 4))), {x: (-1, 2), y: (3, 4)}),
        (Ball((1, 2), 3), {x: (-2, 4), y: (-1, 5)}),
        (Sphere((1, 2), 3), {x: (-2, 4), y: (-1, 5)}),
        (Ellipsoid((0, 0), [[4, 1], [1, 9]]), {x: (-2, 2), y: (-3, 3)}),
        (Simplex(((0, 0), (1, 1))), {x: (0, 1), y: (0, 1)}),
    ],
)
def test_structural_bounds(region, expected):
    assert bounding_cover(region, (x, y)) == expected


def test_polynomial_bounds():
    region = ParametricRegion((t,), ((t, -1, 2),), (t * t, t**3))
    assert bounding_cover(region, (x, y)) == {x: (0, 4), y: (-1, 8)}


def test_transformed_cover():
    region = TransformedRegion(Box(((0, 1),)), (t, t * t), (t,))
    cover = optimization_cover(region, (x, y))
    assert cover.exact and len(cover.charts) == 1
    chart = cover.charts[0]
    parameter = chart.parameters[0]
    assert tuple(expr.subs(parameter, sp.Rational(1, 2)) for expr in chart.mapping) == (
        sp.Rational(1, 2),
        sp.Rational(1, 4),
    )
    assert bounding_cover(region, (x, y)) == {x: (0, 1), y: (0, 1)}


def test_union_retains_points():
    from semialg import Point

    region = BooleanRegion("union", (Simplex(((0, 0), (1, 1))), Point((-1, 0))))
    cover = optimization_cover(region, (x, y))
    assert len(cover.charts) == 2
    assert cover.charts[1].parameters == ()
    assert cover.charts[1].mapping == (-1, 0)
    assert bounding_cover(region, (x, y)) == {x: (-1, 1), y: (0, 1)}


def test_family_root_collision():
    family = region_family(sp.Eq(x * x, a), (x,), (a,))
    expected = {-1: (-1, 0), 0: (0, 1), 1: (0, 2)}
    for parameter, result in expected.items():
        selected = [
            s for s in family.strata if sp.simplify(s.condition.subs(a, parameter)) is sp.true
        ]
        assert len(selected) == 1
        assert (selected[0].fiber_dimension, selected[0].component_count) == result
    assert family.hardt.certified
    assert sp.simplify(family.empty_condition.subs(a, -1)) is sp.true


def test_family_dimension_change():
    family = region_family(sp.And(sp.Eq(a * x, 0), sp.Eq(y, 0)), (x, y), (a,))
    for parameter, dimension in [(0, 1), (1, 0), (-1, 0)]:
        strata = [
            s for s in family.strata if sp.simplify(s.condition.subs(a, parameter)) is sp.true
        ]
        assert len(strata) == 1 and strata[0].fiber_dimension == dimension
        assert strata[0].component_count is None


def test_family_validation():
    with pytest.raises(ValueError):
        region_family(sp.Eq(x, a), (x,), (x,))
    with pytest.raises(ValueError):
        region_family(sp.Eq(x, a), (x,), (t,))


def test_documentation_session():
    import re
    from pathlib import Path

    page = Path(__file__).resolve().parents[1] / "docs/optimization_geometry.md"
    namespace = {"__name__": "__geometry_documentation__"}
    for block in re.findall(r"```python\n(.*?)\n```", page.read_text(encoding="utf-8"), re.DOTALL):
        exec(compile(block, str(page), "exec"), namespace)


def test_nested_limit_bounds():
    region = ParametricRegion((t, a), ((t, 0, a), (a, 0, 1)), (t, a))
    assert bounding_cover(region, (x, y)) == {x: (0, 1), y: (0, 1)}
