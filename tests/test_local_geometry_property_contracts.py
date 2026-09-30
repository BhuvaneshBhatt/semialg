import itertools

import pytest
import sympy as sp

from semialg import (
    find_instance,
    local_components,
    local_germ,
    local_preimage,
    local_sign_strata,
    parameter_strata,
    point_in_closure,
)


@pytest.mark.parametrize(
    ("symbol", "constraint"),
    [
        (sp.Symbol("p", positive=True), lambda x: x > 0),
        (sp.Symbol("n", negative=True), lambda x: x < 0),
        (sp.Symbol("nn", nonnegative=True), lambda x: x >= 0),
        (sp.Symbol("np", nonpositive=True), lambda x: x <= 0),
        (sp.Symbol("nz", nonzero=True, real=True), lambda x: sp.Ne(x, 0)),
    ],
)
def test_witness_satisfies_lowered_symbol_assumption(symbol, constraint):
    point = find_instance(sp.true, (symbol,))
    assert point is not None
    assert sp.simplify(constraint(point[symbol])) is sp.true


def test_multiple_witnesses_respect_assumptions_and_formula():
    x = sp.Symbol("x", positive=True)
    points = find_instance(x < 4, (x,), count=3)
    assert len(points) == 3
    for point in points:
        assert point[x] > 0
        assert point[x] < 4


def test_incompatible_explicit_constraint_has_no_witness():
    x = sp.Symbol("x", positive=True)
    # Unevaluated keeps the contradictory public expression available to Semialg.
    contradiction = sp.Lt(x, 0, evaluate=False)
    assert find_instance(contradiction, (x,)) is None


@pytest.mark.slow
@pytest.mark.parametrize(
    "region, expected",
    [
        (lambda x, y: sp.Eq(y**2, x**2), 4),
        (lambda x, y: sp.Eq(y**2, x**3), 2),
        (lambda x, y: sp.Eq(y**2, x**4), 4),
        (lambda x, y: sp.And(sp.Eq(y, 0), sp.Ne(x, 0)), 2),
    ],
)
def test_local_component_singularity_families(region, expected):
    x, y = sp.symbols("x y", real=True)
    assert len(local_components(region(x, y), (0, 0), (x, y))) == expected


@pytest.mark.slow
def test_local_components_are_translation_equivariant():
    x, y = sp.symbols("x y", real=True)
    base = sp.Eq(y**2, x**2)
    moved = base.xreplace({x: x - 3, y: y + 2})
    assert len(local_components(base, (0, 0), (x, y))) == len(
        local_components(moved, (3, -2), (x, y))
    )


def test_local_preimage_unique_algebraic_fiber():
    x, y = sp.symbols("x y", real=True)
    target = local_germ(sp.Eq(y, 2), (2,), (y,))
    preimage = local_preimage(target, (x**3,), (x,))
    # Restricting the source map to x >= 0 is represented in the map's target
    # formula, not in the fiber equation, so the unrestricted map is nonunique.
    assert preimage.point[x] ** 3 == 2


def test_local_preimage_empty_and_positive_dimensional_fibers():
    x, z, y = sp.symbols("x z y", real=True)
    impossible = local_germ(sp.Eq(y, -1), (-1,), (y,))
    with pytest.raises(ValueError, match="no exact source"):
        local_preimage(impossible, (x**2,), (x,))

    target = local_germ(sp.Eq(y, 0), (0,), (y,))
    with pytest.raises(ValueError, match="unique"):
        local_preimage(target, (x * 0,), (x, z))


def _exhaustive_signs(expressions, region, point, variables):
    result = set()
    for signs in itertools.product((-1, 0, 1), repeat=len(expressions)):
        atoms = [
            expr < 0 if sign < 0 else sp.Eq(expr, 0) if sign == 0 else expr > 0
            for expr, sign in zip(expressions, signs, strict=True)
        ]
        formula = sp.And(region, *atoms)
        if point_in_closure(formula, dict(zip(variables, point, strict=True)), variables):
            result.add(signs)
    return result


@pytest.mark.parametrize(
    "expressions",
    [
        lambda x, y: (x, y),
        lambda x, y: (x, x**2),
        lambda x, y: (x + y, x - y),
        lambda x, y: (x, 2 * x, x**2 + y**2),
    ],
)
def test_incremental_sign_strata_match_exhaustive_reference(expressions):
    x, y = sp.symbols("x y", real=True)
    exprs = expressions(x, y)
    actual = {s.signs for s in local_sign_strata(exprs, sp.true, (0, 0), (x, y))}
    assert actual == _exhaustive_signs(exprs, sp.true, (0, 0), (x, y))


def test_sign_strata_are_permutation_equivariant():
    x, y = sp.symbols("x y", real=True)
    forward = {s.signs for s in local_sign_strata((x, y), sp.true, (0, 0), (x, y))}
    reverse = {tuple(reversed(s.signs)) for s in local_sign_strata((y, x), sp.true, (0, 0), (x, y))}
    assert forward == reverse


def test_parameter_strata_are_disjoint_and_cover_declared_cases():
    x, a = sp.symbols("x a", real=True)
    strata = parameter_strata(
        sp.Eq(x**2, a),
        (x,),
        (a,),
        conditions_by_value={
            "zero": sp.Eq(x, 0),
            "positive": x > 0,
            "negative": x < 0,
        },
    )
    for left, right in itertools.combinations(strata, 2):
        assert find_instance(sp.And(left.condition, right.condition), (a,)) is None
    assert any(sp.simplify(s.condition.subs(a, 0)) is sp.true for s in strata)
    assert any(sp.simplify(s.condition.subs(a, 4)) is sp.true for s in strata)
