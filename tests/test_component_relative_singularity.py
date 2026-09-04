import sympy as sp

from semialg.algebraic_decomposition import equidimensional_decomposition
from semialg.region_analysis import region_singular_locus


def test_lower_dimensional_cusp_uses_its_own_expected_jacobian_rank():
    x, y, z = sp.symbols("x y z", real=True)
    # V(x) union {x=1, z^2=y^3}.  The plane has codimension one while the
    # cusp curve has codimension two.  A single global threshold of one misses
    # the singular tip of the lower-dimensional component.
    equations = (x * (x - 1), x * (z**2 - y**3))
    decomposition = equidimensional_decomposition(equations, (x, y, z))

    assert decomposition.complete
    assert [(piece.dimension, piece.codimension) for piece in decomposition.pieces] == [
        (2, 1),
        (1, 2),
    ]

    region = sp.And(*(sp.Eq(equation, 0) for equation in equations))
    singular = region_singular_locus(region, (x, y, z))

    assert singular.subs({x: 1, y: 0, z: 0}) is sp.true
    assert singular.subs({x: 1, y: 1, z: 1}) is sp.false
    assert singular.subs({x: 0, y: 0, z: 0}) is sp.false


def test_regular_isolated_component_is_not_singular_from_dimension_gap():
    x, y, z = sp.symbols("x y z", real=True)
    # z=0 is a plane; the other component is the isolated point (0, 0, 1).
    equations = (z * x, z * y, z * (z - 1))
    decomposition = equidimensional_decomposition(equations, (x, y, z))

    assert decomposition.complete
    assert [(piece.dimension, piece.codimension) for piece in decomposition.pieces] == [
        (2, 1),
        (0, 3),
    ]

    region = sp.And(*(sp.Eq(equation, 0) for equation in equations))
    singular = region_singular_locus(region, (x, y, z))

    assert singular.subs({x: 0, y: 0, z: 1}) is sp.false
    assert singular.subs({x: 1, y: 1, z: 0}) is sp.false


def test_reduced_repeated_equations_do_not_create_singularity():
    x = sp.symbols("x", real=True)
    singular = region_singular_locus(sp.Eq(x**4, 0), (x,))
    assert singular is sp.false


def test_real_radical_dimension_drop_prevents_false_hypersurface_singularity():
    x, y, z = sp.symbols("x y z", real=True)
    # Over C this is a hypersurface and its defining gradient vanishes along
    # the real locus. Over R, however, x^2+y^2=0 is exactly the smooth line
    # x=y=0. Component regularity must use the certified real radical/codim 2.
    singular = region_singular_locus(sp.Eq(x**2 + y**2, 0), (x, y, z))
    assert singular is sp.false
