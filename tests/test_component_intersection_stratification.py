import sympy as sp

from semialg.region_analysis import (
    algebraic_singularity_stratification,
    region_singular_locus,
)


def _truth(formula, substitutions):
    return sp.simplify(formula.subs(substitutions))


def test_crossing_lines_are_union_induced_not_intrinsic_singularities():
    x, y = sp.symbols("x y", real=True)
    result = algebraic_singularity_stratification((x * y,), (x, y))

    assert result.complete
    assert len(result.components) == 2
    assert all(component.intrinsic_singular is sp.false for component in result.components)
    assert len(result.intersections) == 1
    stratum = result.intersections[0]
    assert stratum.component_indices == (0, 1)
    assert stratum.dimension == 0
    assert _truth(result.union_singular, {x: 0, y: 0}) is sp.true
    assert _truth(result.union_singular, {x: 1, y: 0}) is sp.false


def test_tangent_smooth_components_are_still_singular_as_a_union():
    x, y = sp.symbols("x y", real=True)
    result = algebraic_singularity_stratification((y * (y - x**2),), (x, y))

    assert all(component.intrinsic_singular is sp.false for component in result.components)
    assert len(result.intersections) == 1
    assert result.intersections[0].dimension == 0
    assert _truth(result.singular_formula, {x: 0, y: 0}) is sp.true
    assert _truth(result.singular_formula, {x: 1, y: 0}) is sp.false
    assert _truth(result.singular_formula, {x: 1, y: 1}) is sp.false


def test_triple_intersection_has_exact_incidence_support():
    x, y = sp.symbols("x y", real=True)
    result = algebraic_singularity_stratification((x * y * (x - y),), (x, y))

    assert len(result.components) == 3
    assert len(result.intersections) == 1
    stratum = result.intersections[0]
    assert stratum.component_indices == (0, 1, 2)
    assert stratum.order == 3
    assert stratum.dimension == 0
    assert _truth(stratum.formula, {x: 0, y: 0}) is sp.true


def test_disjoint_components_have_no_union_singular_strata():
    x, y = sp.symbols("x y", real=True)
    result = algebraic_singularity_stratification((x * (x - 1),), (x, y))

    assert len(result.components) == 2
    assert result.intersections == ()
    assert sp.simplify(result.union_singular) is sp.false
    assert sp.simplify(result.singular_formula) is sp.false


def test_contained_redundant_branch_does_not_create_union_singularity():
    x, y = sp.symbols("x y", real=True)
    result = algebraic_singularity_stratification((x, x * y), (x, y))

    assert len(result.components) == 1
    assert result.intersections == ()
    assert result.singular_formula is sp.false


def test_intrinsic_and_union_singularities_are_recorded_separately():
    x, y, z = sp.symbols("x y z", real=True)
    equations = (x * (x - 1), x * (z**2 - y**3))
    result = algebraic_singularity_stratification(equations, (x, y, z))

    assert len(result.components) == 2
    assert _truth(result.intrinsic_singular, {x: 1, y: 0, z: 0}) is sp.true
    assert _truth(result.union_singular, {x: 1, y: 0, z: 0}) is sp.false
    assert _truth(result.singular_formula, {x: 1, y: 0, z: 0}) is sp.true
    assert _truth(result.singular_formula, {x: 0, y: 2, z: 3}) is sp.false


def test_ordinary_regular_component_locus_excludes_intersections():
    x, y = sp.symbols("x y", real=True)
    result = algebraic_singularity_stratification((x * y,), (x, y))

    x_axis = next(component for component in result.components if component.equations == (y,))
    assert _truth(x_axis.ordinary_regular_formula, {x: 1, y: 0}) is sp.true
    assert _truth(x_axis.ordinary_regular_formula, {x: 0, y: 0}) is sp.false


def test_region_singular_locus_uses_stratified_union_semantics():
    x, y = sp.symbols("x y", real=True)
    singular = region_singular_locus(sp.Eq(y * (y - x**2), 0), (x, y))

    assert _truth(singular, {x: 0, y: 0}) is sp.true
    assert _truth(singular, {x: 2, y: 0}) is sp.false
    assert _truth(singular, {x: 2, y: 4}) is sp.false


def test_incidence_metadata_cap_does_not_change_exact_singular_formula():
    x, y = sp.symbols("x y", real=True)
    equation = x * y * (x - y)
    result = algebraic_singularity_stratification((equation,), (x, y), max_incidence_subsets=1)

    assert not result.complete
    assert result.intersections == ()
    assert _truth(result.singular_formula, {x: 0, y: 0}) is sp.true
    assert _truth(result.singular_formula, {x: 1, y: 0}) is sp.false


def test_pair_and_triple_incidence_strata_are_disjoint_and_dimensioned():
    x, y, z = sp.symbols("x y z", real=True)
    result = algebraic_singularity_stratification((x * y * z,), (x, y, z))

    pair_strata = [stratum for stratum in result.intersections if stratum.order == 2]
    triple_strata = [stratum for stratum in result.intersections if stratum.order == 3]
    assert len(pair_strata) == 3
    assert all(stratum.dimension == 1 for stratum in pair_strata)
    assert len(triple_strata) == 1
    assert triple_strata[0].dimension == 0

    assert all(_truth(stratum.formula, {x: 0, y: 0, z: 0}) is sp.false for stratum in pair_strata)
    assert _truth(triple_strata[0].formula, {x: 0, y: 0, z: 0}) is sp.true
