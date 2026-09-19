import sympy as sp

import semialg.algebraic_decomposition as decomposition


def test_algebraic_decomposition_retains_nonreal_components():
    """Algebraic decomposition must not discard components with no real locus."""
    x, y = sp.symbols("x y")
    equations = (x**2 + y**2 + 1,)

    geometric = decomposition.equidimensional_decomposition(equations, (x, y))
    algebraic = decomposition._equidimensional_decomposition_impl(
        equations, (x, y), _algebraic_semantics=True
    )

    assert geometric.complete
    assert geometric.pieces == ()
    assert algebraic.complete
    assert len(algebraic.pieces) == 1
    assert algebraic.pieces[0].dimension == 1


def test_fiber_specialization_keeps_opaque_parent_tower_native():
    from semialg.algebraic.samples import AlgebraicRoot, FiberRootContext, RationalInterval
    from semialg.cad_algorithms.decomposition import (
        TowerFiberSpecialization,
        _specialize_fiber_polynomial,
    )

    r, a, b = sp.symbols("r a b")
    r_root = AlgebraicRoot(
        sp.Poly(3 * r**2 - 6 * r + 2, r), RationalInterval(0, sp.Rational(3, 2)), 0
    )
    a_context = FiberRootContext(
        sp.Poly(
            9 * r**2 * a**2 - 6 * r**2 * a - 15 * r * a**2 + 12 * r * a - 2 * r + 6 * a**2 - 4 * a,
            r,
            a,
        ),
        (r, a),
        (r_root,),
        a,
    )
    a_root = AlgebraicRoot(
        sp.Poly(
            (9 * r**2 - 15 * r + 6) * a**2 + (-6 * r**2 + 12 * r - 4) * a - 2 * r,
            a,
            domain=sp.QQ.frac_field(r),
        ),
        RationalInterval(-1, 0),
        0,
        fiber_context=a_context,
    )
    poly = sp.Poly(
        3 * r * a**2 - 4 * r * a + 3 * r * b**2 + r - 2 * a**2 + 2 * a - 2 * b**2, r, a, b
    )
    result = _specialize_fiber_polynomial(poly, (r, a, b), (r_root, a_root))
    assert isinstance(result, TowerFiberSpecialization)
    assert result.variables == (r, a, b)
    assert result.parent_samples == (r_root, a_root)
    assert not result.expression.has(sp.CRootOf)
