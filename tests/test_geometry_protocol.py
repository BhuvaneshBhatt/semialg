import pytest
import sympy as sp

from semialg import BallRegion, BoxRegion, Geometry, PointRegion, SemialgebraicRegion, SimplexRegion


def test_standard_regions_share_geometry_protocol():
    regions = (
        PointRegion((1, 2)),
        BoxRegion(((0, 1), (-1, 2))),
        BallRegion((0, 0), 2),
        SimplexRegion(((0, 0), (1, 0), (0, 1))),
    )
    assert all(isinstance(region, Geometry) for region in regions)
    assert [region.intrinsic_dimension for region in regions] == [0, 2, 2, 2]
    assert [region.ambient_dimension() for region in regions] == [2, 2, 2, 2]


def test_as_formula_uses_existing_semialgebraic_lowering():
    x, y = sp.symbols("x y", real=True)
    ball = BallRegion((1, -1), 3)
    formula = ball.as_formula((x, y))
    assert sp.simplify(formula ^ ((x - 1) ** 2 + (y + 1) ** 2 <= 9)) is sp.false


def test_as_formula_can_eliminate_structural_quantifiers():
    x, y = sp.symbols("x y", real=True)
    simplex = SimplexRegion(((0, 0), (1, 0), (0, 1)))
    raw = simplex.as_formula((x, y))
    assert raw.has(sp.Symbol) and "Exists" in type(raw).__name__
    quantifier_free = simplex.as_formula((x, y), eliminate=True)
    assert "Exists" not in sp.sstr(quantifier_free)
    region = SemialgebraicRegion(quantifier_free, (x, y))
    assert region.contains((sp.Rational(1, 4), sp.Rational(1, 4))) is True
    assert region.contains((1, 1)) is False


def test_contains_and_boundary_delegate_to_symbolic_geometry():
    x = sp.symbols("x", real=True)
    box = BoxRegion(((0, 1),))
    assert box.contains((sp.Rational(1, 2),)) is True
    assert box.contains((2,)) is False
    boundary = box.boundary((x,), strategy="syntactic")
    assert isinstance(boundary, SemialgebraicRegion)
    assert sp.simplify(boundary.formula ^ sp.Or(sp.Eq(x, 0), sp.Eq(x, 1))) is sp.false


def test_affine_hull_is_optional_capability():
    with pytest.raises(NotImplementedError, match="structural affine hull"):
        BallRegion((0, 0), 1).affine_hull()
