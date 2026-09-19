import pytest
import sympy as sp

from semialg import AffineSpace, Cone, ConeRegion, Cylinder, CylinderRegion, Line


def test_cylinder_canonical_geometry_and_formula():
    x, y, z = sp.symbols("x y z", real=True)
    region = Cylinder((0, 0, 0), (0, 0, 3), 2)
    assert region.dimension() == 3
    assert region.ambient_dimension() == 3
    assert region.axis == (0, 0, 3)
    assert region.height == 3
    assert region.affine_hull() == AffineSpace((0, 0, 0), ((1, 0, 0), (0, 1, 0), (0, 0, 1)))
    formula = region.as_formula((x, y, z), eliminate=True)
    assert bool(formula.subs({x: 2, y: 0, z: 1}))
    assert not bool(formula.subs({x: 3, y: 0, z: 1}))
    assert not bool(formula.subs({x: 0, y: 0, z: 4}))


def test_zero_radius_cylinder_is_axis_segment_with_line_affine_hull():
    region = Cylinder((0, 0), (2, 0), 0)
    assert region.dimension() == 1
    assert region.affine_hull() == Line((0, 0), (2, 0))
    assert region.contains((1, 0))
    assert not region.contains((3, 0))


def test_cone_uses_base_at_start_and_apex_at_end():
    x, y, z = sp.symbols("x y z", real=True)
    region = Cone((0, 0, 0), (0, 0, 2), 2)
    formula = region.as_formula((x, y, z), eliminate=True)
    assert bool(formula.subs({x: 2, y: 0, z: 0}))
    assert not bool(formula.subs({x: 2, y: 0, z: 2}))
    assert bool(formula.subs({x: 0, y: 0, z: 2}))
    assert not bool(formula.subs({x: 1, y: 0, z: sp.Rational(3, 2)}))


def test_canonical_axis_validation_and_symbolic_condition():
    a = sp.symbols("a", real=True)
    with pytest.raises(ValueError, match="axis endpoints must be distinct"):
        Cylinder((0, 0), (0, 0), 1)
    with pytest.raises(ValueError, match="axis endpoints must be distinct"):
        Cone((0, 0), (0, 0), 1)
    with pytest.raises(ValueError, match="ambient dimension"):
        Cylinder((0,), (1,), 1)
    region = Cylinder((0, 0), (a, 0), 1)
    assert region.construction_conditions == (a**2 > 0,)


def test_region_conversion_bridges_preserve_canonical_semantics():
    cylinder = Cylinder.from_region(CylinderRegion((0, 0, 0), (0, 0, 1), 3))
    cone = Cone.from_region(ConeRegion((0, 0, 0), (0, 0, 1), 3))
    assert isinstance(cylinder, Cylinder)
    assert isinstance(cone, Cone)
    assert cylinder.radius == cone.radius == 3
    with pytest.raises(TypeError):
        Cylinder.from_region(ConeRegion((0, 0, 0), (0, 0, 1), 1))


def test_existing_cone_region_formula_matches_specialized_orientation():
    x, y, z = sp.symbols("x y z", real=True)
    region = ConeRegion((0, 0, 0), (0, 0, 2), 2)
    formula = region.as_formula((x, y, z), eliminate=True)
    assert bool(formula.subs({x: 2, y: 0, z: 0}))
    assert not bool(formula.subs({x: 2, y: 0, z: 2}))
