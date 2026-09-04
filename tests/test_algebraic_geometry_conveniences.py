import sympy as sp

from semialg import is_singular, is_smooth, tangent_dimension


def test_smoothness_and_point_singularity_conveniences():
    x, y = sp.symbols("x y", real=True)
    assert is_smooth([x**2 + y**2 - 1], [x, y])
    cusp = y**2 - x**3
    assert not is_smooth([cusp], [x, y])
    assert is_singular([cusp], {x: 0, y: 0}, [x, y])
    assert not is_singular([cusp], {x: 1, y: 1}, [x, y])
    assert tangent_dimension([cusp], {x: 0, y: 0}, [x, y]) == 2
