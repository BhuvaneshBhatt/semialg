import sympy as sp

from semialg import local_branch_geometry


def test_transverse_crossing_has_two_smooth_branches():
    x, y = sp.symbols("x y")
    result = local_branch_geometry((x * y,), {x: 0, y: 0}, (x, y))
    assert result.component_indices == (0, 1)
    assert result.singular
    assert all(branch.tangent_dimension_excess == 0 for branch in result.branches)
    assert all(branch.multiplicity == 1 for branch in result.branches)
    assert result.intersection is not None
    assert result.intersection.transverse is True
    assert result.intersection.dimension == 0


def test_tangent_branches_are_nontransverse():
    x, y = sp.symbols("x y")
    # Two smooth branches y=0 and y=x**2 have the same tangent line at the origin.
    result = local_branch_geometry((y * (y - x**2),), {x: 0, y: 0}, (x, y))
    assert len(result.branches) == 2
    assert all(branch.tangent_dimension_excess == 0 for branch in result.branches)
    assert result.intersection is not None
    assert result.intersection.transverse is False


def test_cusp_has_tangent_excess_and_multiplicity_two():
    x, y = sp.symbols("x y")
    result = local_branch_geometry((y**2 - x**3,), {x: 0, y: 0}, (x, y))
    assert len(result.branches) == 1
    branch = result.branches[0]
    assert branch.tangent_space.dimension == 2
    assert branch.component.dimension == 1
    assert branch.tangent_dimension_excess == 1
    assert branch.multiplicity == 2
    assert branch.tangent_cone_degree == 2
    assert result.intersection is None
    assert result.singular


def test_nonradical_presentation_does_not_inflate_local_geometry():
    x, y = sp.symbols("x y")
    result = local_branch_geometry((y**2,), {x: 0, y: 0}, (x, y))
    assert len(result.branches) == 1
    branch = result.branches[0]
    assert branch.component.equations == (y,)
    assert branch.tangent_dimension_excess == 0
    assert branch.multiplicity == 1
    assert result.singular is False
