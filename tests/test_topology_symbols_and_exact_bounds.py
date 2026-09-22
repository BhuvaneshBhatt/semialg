from __future__ import annotations

import pytest
import sympy as sp

from semialg.algebraic.rational_univariate.signs import sign_of_algebraic_expression
from semialg.algebraic_geometry import tangent_cone, tangent_space
from semialg.geometry_queries import (
    distance_to_region,
    fiber,
    path_between,
    semialgebraic_projection,
)
from semialg.instances.witness_generation import safe_int_ceiling_bound, safe_integer_floor_bound
from semialg.region_transformations import region_image
from semialg.regions.operations import region_closure, region_dimension, region_interior
from semialg.simplify.implication import _is_unsatisfiable_cached


def test_topology_ne_is_open_and_dense():
    x = sp.symbols("x", real=True)
    assert region_closure(sp.Ne(x, 0), [x]) is sp.true
    interior = region_interior(sp.Ne(x, 0), [x])
    assert sp.simplify(interior.subs(x, 1)) is sp.true
    assert sp.simplify(interior.subs(x, 0)) is sp.false


def test_strict_polynomial_closure_handles_even_factor_content():
    x, y = sp.symbols("x y", real=True)
    closure = region_closure(x**2 * y > 0, [x, y])
    # The closure is y >= 0 and includes the x=0 fibre.
    assert sp.simplify(closure.subs({x: 0, y: 1})) is sp.true
    assert sp.simplify(closure.subs({x: 0, y: -1})) is sp.false
    assert not any(
        getattr(node, "func", None).__name__ == "root_of" and node.args[0] == 0
        for node in sp.preorder_traversal(closure)
        if isinstance(node, sp.Basic)
    )


def test_lower_dimensional_3d_region_has_empty_interior_and_exact_dimension():
    x, y, z = sp.symbols("x y z", real=True)
    region = sp.And(sp.Eq(x * y, 0), sp.Eq(x * z, 0))
    assert region_interior(region, [x, y, z]) is sp.false
    assert region_dimension(region, [x, y, z]) == 2


def test_integer_endpoint_bounds_never_round_through_float():
    n = 10**20
    above = sp.Integer(n) + sp.sqrt(2) / 10
    below = sp.Integer(n) - sp.sqrt(2) / 10
    assert safe_integer_floor_bound(above, strict=True) == n
    assert safe_integer_floor_bound(above) == n
    assert safe_int_ceiling_bound(below, strict=True) == n
    assert safe_int_ceiling_bound(below) == n
    assert safe_integer_floor_bound(sp.Integer(n), strict=True) == n - 1
    assert safe_int_ceiling_bound(sp.Integer(n), strict=True) == n + 1


def test_projection_resolves_string_variable_against_contextual_symbol():
    x = sp.Symbol("x", positive=True)
    y = sp.Symbol("y", real=True)
    result = semialgebraic_projection(sp.And(x > 1, y > x), ["x"], variables=[x, y])
    assert sp.simplify(sp.Equivalent(result, y > 1)) is sp.true


def test_image_preserves_unspecified_free_parameter():
    x, a, y = sp.symbols("x a y", real=True)
    result = region_image(
        sp.And(x >= 0, x <= a),
        x,
        variables=[x],
        image_variables=[y],
    )
    assert a in result.free_symbols
    assert x not in result.free_symbols
    for av, yv, expected in [(2, 1, True), (2, 3, False), (-1, 0, False)]:
        value = sp.simplify(result.subs({a: av, y: yv}))
        assert bool(value) is expected


def test_fiber_rejects_ambiguous_same_name_symbol():
    xr = sp.Symbol("x", real=True)
    xp = sp.Symbol("x", positive=True)
    formula = sp.And(xr < 4, xp < 3)
    with pytest.raises(ValueError, match="ambiguous"):
        fiber(formula, {"x": 1})
    assert fiber(formula, {xr: 1}) == (xp < 3)


def test_tangent_directions_do_not_collide_for_same_name_symbols():
    xr = sp.Symbol("x", real=True)
    xp = sp.Symbol("x", positive=True)
    equations = (xr, xp)
    point = {xr: 0, xp: 0}
    cone = tangent_cone(equations, point, variables=[xr, xp])
    space = tangent_space(equations, point, variables=[xr, xp])
    assert len(set(cone.direction_variables)) == 2
    assert len(set(space.equations[0].free_symbols | space.equations[1].free_symbols)) == 2


def test_implication_cache_preserves_assumption_distinct_symbols():
    xr = sp.Symbol("x", real=True)
    xp = sp.Symbol("x", positive=True)
    expr = sp.And(xr < 0, xp > 0)
    variables = (xr, xp)
    # They are two independent CAD variables, not one reconstructed real x.
    assert _is_unsatisfiable_cached(expr, variables) is False


def test_algebraic_sign_handles_exact_algebraic_expression():
    expression = sp.sqrt(2) - 1
    assert sign_of_algebraic_expression(expression) == 1


def test_shared_point_normalization_accepts_contextual_string_keys():
    x = sp.Symbol("x", positive=True)
    region = x >= 1
    assert distance_to_region({"x": 0}, region, [x]) == 1
    path = path_between(sp.And(x >= 1, x <= 3), {"x": 1}, {"x": 2}, [x])
    assert path.connected
