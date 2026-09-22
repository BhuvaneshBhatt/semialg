import sympy as sp

from semialg import (
    local_dimension_strata,
    minimal_prime_intersections,
    stratified_singular_geometry,
)


def test_crossing_lines_are_component_intersection_not_intrinsic_singularities():
    x, y = sp.symbols("x y")
    result = stratified_singular_geometry((x * y,), (x, y))

    assert len(result.components) == 2
    assert all(item.formula is sp.false for item in result.component_singular_loci)
    assert len(result.minimal_prime_intersections) == 1
    crossing = result.minimal_prime_intersections[0]
    assert crossing.dimension == 0
    assert set(crossing.equations) == {x, y}

    crossing_strata = [item for item in result.singular_strata if item.component_intersection]
    assert len(crossing_strata) == 1
    assert crossing_strata[0].intrinsically_singular is False
    assert sp.simplify_logic(crossing_strata[0].formula ^ (sp.Eq(x, 0) & sp.Eq(y, 0))) is sp.false


def test_local_dimension_strata_handle_mixed_dimensional_components():
    x, y, z = sp.symbols("x y z")
    # V(x*y, x*z) = plane x=0 union line y=z=0.
    strata = local_dimension_strata((x * y, x * z), (x, y, z))
    assert [item.dimension for item in strata] == [2, 1]
    assert sp.simplify_logic(strata[0].formula ^ sp.Eq(x, 0)) is sp.false
    expected_line_only = sp.Eq(y, 0) & sp.Eq(z, 0) & sp.Ne(x, 0)
    assert sp.simplify_logic(strata[1].formula ^ expected_line_only) is sp.false


def test_minimal_prime_intersections_are_deterministic_groebner_sums():
    x, y, z = sp.symbols("x y z")
    intersections = minimal_prime_intersections((x * y, x * z), (x, y, z))
    assert len(intersections) == 1
    assert intersections[0].component_indices == (0, 1)
    assert intersections[0].equations == (x, y, z)
    assert intersections[0].dimension == 0


def test_nonradical_input_is_reduced_before_stratification():
    x = sp.symbols("x")
    result = stratified_singular_geometry((x**2,), (x,))
    assert len(result.components) == 1
    assert result.components[0].equations == (x,)
    assert result.component_singular_loci[0].formula is sp.false
    assert result.minimal_prime_intersections == ()
    assert result.singular_strata == ()
