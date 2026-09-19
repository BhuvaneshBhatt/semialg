import sympy as sp

from semialg.algebraic_decomposition import (
    equidimensional_decomposition,
    verify_decomposition_certificate,
)


def _zero_set_formula(piece):
    return sp.And(*(sp.Eq(eq, 0) for eq in piece.equations)) if piece.equations else sp.true


def test_repeated_generator_has_reduced_set_semantics():
    x = sp.symbols("x", real=True)
    result = equidimensional_decomposition((x**4,), (x,))

    assert result.complete
    assert result.dimensions == (0,)
    assert len(result.pieces) == 1
    assert sp.simplify(result.pieces[0].equations[0] / x) in (1, -1)


def test_crossing_axes_split_into_pure_dimensional_linear_components():
    x, y = sp.symbols("x y", real=True)
    result = equidimensional_decomposition((x * y,), (x, y))

    assert result.complete
    assert result.dimensions == (1,)
    assert len(result.pieces) == 2
    assert all(piece.dimension == 1 and piece.equidimensional for piece in result.pieces)
    formulas = {_zero_set_formula(piece) for piece in result.pieces}
    assert formulas == {sp.Eq(x, 0), sp.Eq(y, 0)}


def test_non_equidimensional_plane_union_line_is_separated_by_dimension():
    x, y, z = sp.symbols("x y z", real=True)
    # V(x*y, x*z) = V(x) union V(y, z).
    result = equidimensional_decomposition((x * y, x * z), (x, y, z))

    assert result.complete
    assert result.dimensions == (2, 1)
    assert [piece.dimension for piece in result.pieces] == [2, 1]
    assert all(piece.equidimensional for piece in result.pieces)


def test_surface_union_isolated_point_is_not_misclassified_by_global_dimension():
    x, y, z = sp.symbols("x y z", real=True)
    # V(z*x, z*y, z*(z-1)) = V(z) union V(x, y, z-1).
    result = equidimensional_decomposition(
        (z * x, z * y, z * (z - 1)),
        (x, y, z),
    )

    assert result.complete
    assert result.dimensions == (2, 0)
    assert [piece.dimension for piece in result.pieces] == [2, 0]


def test_disjoint_parallel_planes_remain_two_same_dimension_components():
    x, y, z = sp.symbols("x y z", real=True)
    result = equidimensional_decomposition((z * (z - 1),), (x, y, z))

    assert result.complete
    assert result.dimensions == (2,)
    assert len(result.pieces) == 2
    assert all(piece.dimension == 2 for piece in result.pieces)


def test_empty_locus_has_no_pieces():
    x = sp.symbols("x", real=True)
    result = equidimensional_decomposition((x, x - 1), (x,))

    assert result.complete
    assert result.pieces == ()
    assert result.dimensions == ()


def test_piece_budget_preserves_exact_cover_but_marks_unresolved_piece():
    x, y = sp.symbols("x y", real=True)
    result = equidimensional_decomposition((x * y * (x - 1),), (x, y), max_pieces=1)

    assert len(result.pieces) == 1
    # The unsplit principal hypersurface is still equidimensional as a reduced
    # set, so a tight branch budget does not make this particular case unsafe.
    assert result.complete
    assert result.pieces[0].dimension == 1


def test_nonradical_zero_dimensional_generators_reduce_to_one_point_piece():
    x, y = sp.symbols("x y", real=True)
    result = equidimensional_decomposition((x**2, y**3), (x, y))

    assert result.complete
    assert result.dimensions == (0,)
    assert len(result.pieces) == 1
    assert set(result.pieces[0].equations) == {x, y}


def test_budget_limited_non_equidimensional_piece_is_not_certified_complete():
    x, y, z = sp.symbols("x y z", real=True)
    result = equidimensional_decomposition((x * y, x * z), (x, y, z), max_pieces=1)

    assert not result.complete
    assert len(result.pieces) == 1
    assert result.pieces[0].dimension == 2
    assert not result.pieces[0].equidimensional


def test_real_reduction_collapses_sum_of_squares_hypersurface_to_regular_point_piece():
    x, y = sp.symbols("x y", real=True)
    result = equidimensional_decomposition((x**2 + y**2,), (x, y))

    assert result.complete
    assert result.dimensions == (0,)
    assert len(result.pieces) == 1
    piece = result.pieces[0]
    assert piece.dimension == 0
    assert piece.algebraic_dimension == 1
    assert set(piece.equations) == {x, y}


def test_monomial_minimal_primes_avoid_branch_budget_false_incompleteness():
    x, y, z = sp.symbols("x y z", real=True)
    # The edge ideal of a triangle has three minimal coordinate primes:
    # (x, y), (x, z), and (y, z). A naive factor-branch traversal creates a
    # duplicate intermediate branch and can exhaust a budget of three before
    # proving completeness. The monomial minimal-prime certificate does not.
    result = equidimensional_decomposition(
        (x * y, x * z, y * z),
        (x, y, z),
        max_pieces=3,
    )

    assert result.complete
    assert result.dimensions == (1,)
    assert len(result.pieces) == 3
    assert {frozenset(piece.equations) for piece in result.pieces} == {
        frozenset((x, y)),
        frozenset((x, z)),
        frozenset((y, z)),
    }


def test_monomial_minimal_primes_remove_embedded_component_under_reduced_semantics():
    x, y = sp.symbols("x y", real=True)
    # I=(x^2, x*y) has radical (x). The associated coordinate prime (x, y)
    # is embedded and must not survive reduced-set decomposition.
    result = equidimensional_decomposition((x**2, x * y), (x, y))

    assert result.complete
    assert result.dimensions == (1,)
    assert len(result.pieces) == 1
    assert set(result.pieces[0].equations) == {x}


def test_nonmonomial_factor_incidence_decomposes_surface_union_line_exactly():
    x, y, z = sp.symbols("x y z", real=True)
    circle_cylinder = x**2 + y**2 - 1
    # V(f*z, f*x) = V(f) union V(x, z).  The shared nonlinear factor is not
    # monomial, so this exercises the generalized factor-incidence certificate.
    result = equidimensional_decomposition(
        (circle_cylinder * z, circle_cylinder * x),
        (x, y, z),
        max_pieces=2,
    )

    assert result.complete
    assert result.dimensions == (2, 1)
    assert len(result.pieces) == 2
    assert {frozenset(piece.equations) for piece in result.pieces} == {
        frozenset((circle_cylinder,)),
        frozenset((x, z)),
    }


def test_nonmonomial_factor_incidence_removes_embedded_reduced_set_branch():
    x, y, z = sp.symbols("x y z", real=True)
    f = x**2 + y**2 - 1
    # V(f^2, f*z) = V(f).  The candidate V(f, z) is contained in V(f) and
    # must not survive reduced-set decomposition.
    result = equidimensional_decomposition((f**2, f * z), (x, y, z))

    assert result.complete
    assert result.dimensions == (2,)
    assert len(result.pieces) == 1
    assert set(result.pieces[0].equations) == {f}


def test_factor_incidence_budget_counts_minimal_candidates_not_naive_branches():
    x, y, z = sp.symbols("x y z", real=True)
    f = x**2 + y**2 - 1
    # Minimal factor covers are {f} and {x,z}; a budget of two is sufficient.
    result = equidimensional_decomposition((f * x, f * z), (x, y, z), max_pieces=2)

    assert result.complete
    assert len(result.pieces) == 2


def test_rabinowitsch_saturation_exactly_removes_component_in_splitter_hyperplane():
    from semialg.algebraic.equality_ideal import EqualityIdealContext
    from semialg.algebraic_decomposition import _saturate_generators

    x, y, z = sp.symbols("x y z")
    # V(xy, xz) = V(x) union V(y,z). Saturating by x removes V(x).
    saturated = _saturate_generators((x * y, x * z), x, (x, y, z))
    assert saturated is not None
    context = EqualityIdealContext.build(saturated, (x, y, z))
    target = EqualityIdealContext.build((y, z), (x, y, z))
    assert all(context.in_radical(g) for g in target.generators)
    assert all(target.in_radical(g) for g in context.generators)


def test_saturation_split_is_exact_and_both_branches_are_strict():
    from semialg.algebraic.equality_ideal import EqualityIdealContext
    from semialg.algebraic_decomposition import _saturation_split_generators

    x, y, z = sp.symbols("x y z")
    context = EqualityIdealContext.build((x * (x - 1), y * (x - 1), x * z, y * z), (x, y, z))
    split = _saturation_split_generators(context)
    assert split is not None
    branches = [EqualityIdealContext.build(gens, (x, y, z)) for gens in split]

    expected = [
        EqualityIdealContext.build((x, y), (x, y, z)),
        EqualityIdealContext.build((x - 1, z), (x, y, z)),
    ]
    assert all(
        any(
            all(branch.in_radical(g) for g in target.generators)
            and all(target.in_radical(g) for g in branch.generators)
            for branch in branches
        )
        for target in expected
    )
    assert all(
        not (
            all(context.in_radical(g) for g in branch.generators)
            and all(branch.in_radical(g) for g in context.generators)
        )
        for branch in branches
    )


def test_saturation_splitter_does_not_split_an_irreducible_hypersurface():
    from semialg.algebraic.equality_ideal import EqualityIdealContext
    from semialg.algebraic_decomposition import _saturation_split_generators

    x, y = sp.symbols("x y")
    context = EqualityIdealContext.build((x**2 + y**2 - 1,), (x, y))
    assert _saturation_split_generators(context) is None


def test_real_radical_dimension_drop_to_positive_dimensional_linear_locus():
    x, y, z = sp.symbols("x y z", real=True)
    # Complex V(x^2+y^2) is a hypersurface of dimension 2.  Over R the same
    # equation is exactly x=y=0, a line. The decomposition must strengthen
    # the real ideal rather than relabel the complex hypersurface.
    result = equidimensional_decomposition((x**2 + y**2,), (x, y, z))

    assert result.complete
    assert result.dimensions == (1,)
    assert len(result.pieces) == 1
    piece = result.pieces[0]
    assert set(piece.equations) == {x, y}
    assert piece.dimension == 1
    assert piece.algebraic_dimension == 2
    assert piece.equidimensional


def test_real_radical_even_monomial_sum_can_create_multiple_real_components():
    x, y, z = sp.symbols("x y z", real=True)
    # x^2*y^2 + z^2 = 0 iff xy=0 and z=0 over R, hence the union of the two
    # coordinate lines V(x,z) and V(y,z).  The complex hypersurface itself has
    # dimension two, so this exercises positive-dimensional real dimension drop
    # followed by ordinary reduced component decomposition.
    result = equidimensional_decomposition((x**2 * y**2 + z**2,), (x, y, z))

    assert result.complete
    assert result.dimensions == (1,)
    assert len(result.pieces) == 2
    assert {frozenset(piece.equations) for piece in result.pieces} == {
        frozenset((x, z)),
        frozenset((y, z)),
    }
    assert all(piece.algebraic_dimension == 2 for piece in result.pieces)


def test_real_radical_certifies_empty_definite_component():
    x, y, z = sp.symbols("x y z", real=True)
    result = equidimensional_decomposition((x**2 + y**2 + 1,), (x, y, z))

    assert result.complete
    assert result.pieces == ()


def test_psd_quadratic_real_radical_handles_translated_nondiagonal_zero_set():
    x, y, z = sp.symbols("x y z", real=True)
    q = sp.expand((x + y - 1) ** 2 + (x - y) ** 2)
    result = equidimensional_decomposition((q,), (x, y, z))

    assert result.complete
    assert result.dimensions == (1,)
    assert len(result.pieces) == 1
    piece = result.pieces[0]
    # Compare the returned affine ideal to the expected line x=y=1/2.
    from semialg.algebraic.equality_ideal import EqualityIdealContext

    got = EqualityIdealContext.build(piece.equations, (x, y, z))
    expected = EqualityIdealContext.build((2 * x - 1, 2 * y - 1), (x, y, z))
    assert all(got.in_radical(g) for g in expected.generators)
    assert all(expected.in_radical(g) for g in got.generators)
    assert piece.dimension == 1
    assert piece.algebraic_dimension == 2


def test_real_radical_rewrite_preserves_other_equalities():
    x, y, z = sp.symbols("x y z", real=True)
    result = equidimensional_decomposition((x**2 + y**2, z - 3), (x, y, z))

    assert result.complete
    assert result.dimensions == (0,)
    assert len(result.pieces) == 1
    piece = result.pieces[0]
    from semialg.algebraic.equality_ideal import EqualityIdealContext

    got = EqualityIdealContext.build(piece.equations, (x, y, z))
    expected = EqualityIdealContext.build((x, y, z - 3), (x, y, z))
    assert all(got.in_radical(g) for g in expected.generators)
    assert all(expected.in_radical(g) for g in got.generators)


def test_initial_split_strategy_is_recorded_in_replayable_certificate():
    x, y, z = sp.symbols("x y z")
    decomposition = equidimensional_decomposition((x + y * z, x * y - z**2), (x, y, z))
    assert decomposition.certificate is not None
    assert "initial_split" in decomposition.certificate.methods
    assert verify_decomposition_certificate(decomposition.certificate)


def test_squarefree_regular_chain_certifies_monomial_curve_branch():
    x, y, z = sp.symbols("x y z", real=True)
    result = equidimensional_decomposition((x**2 - y, x * y - z), (x, y, z))

    assert result.complete
    assert result.dimensions == (1,)
    assert len(result.pieces) == 1
    assert result.pieces[0].equidimensional
    assert result.certificate is not None
    assert "squarefree_regular_chain" in result.certificate.methods
    assert verify_decomposition_certificate(result.certificate)


def test_squarefree_regular_chain_does_not_trust_triangular_shape_without_same_ideal():
    from semialg.algebraic.equality_ideal import EqualityIdealContext
    from semialg.algebraic_decomposition import _squarefree_regular_chain_certified

    x, y, z = sp.symbols("x y z", real=True)
    # Triangular shape alone is insufficient: y**2 has vanishing separant on
    # its reduced component.  Exact saturation by the separant must therefore
    # reject this as a squarefree regular-chain certificate.
    context = EqualityIdealContext.build((y**2, z - x), (x, y, z))
    assert not _squarefree_regular_chain_certified(context)
