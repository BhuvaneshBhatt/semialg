import sympy as sp

from semialg.algebraic.roots import isolate_real_roots, refine_isol_intvl
from semialg.algebraic.samples import AlgebraicRoot


def test_mignotte_like_cluster_isolated_exactly():
    x = sp.Symbol("x", real=True)
    polynomial = sp.Poly(x**8 - (32 * x - 1) ** 2, x)

    roots = isolate_real_roots(polynomial)

    assert len(roots) == 4
    assert all(
        int(polynomial.count_roots(root.interval.left, root.interval.right)) == 1 for root in roots
    )
    assert all(
        roots[index].interval.right <= roots[index + 1].interval.left
        for index in range(len(roots) - 1)
    )


def test_newton_refinement_contracts_cluster_interval():
    x = sp.Symbol("x", real=True)
    polynomial = sp.Poly(x**8 - (32 * x - 1) ** 2, x)
    root = isolate_real_roots(polynomial)[2]
    initial_width = root.interval.width

    refined = refine_isol_intvl(root, steps=2)

    assert int(polynomial.count_roots(refined.interval.left, refined.interval.right)) == 1
    assert refined.interval.width <= initial_width / 4


def test_adaptive_proposal_remains_exact_for_algebraic_coefficients():
    x = sp.Symbol("x", real=True)
    polynomial = sp.Poly(x**2 - sp.sqrt(2), x, extension=True)
    root = isolate_real_roots(polynomial)[1]
    opaque_root = AlgebraicRoot(
        root.polynomial, root.interval, root.root_index, root.multiplicity, None
    )

    refined = refine_isol_intvl(opaque_root, steps=3)

    assert int(polynomial.count_roots(refined.interval.left, refined.interval.right)) == 1
    assert refined.interval.width < opaque_root.interval.width
