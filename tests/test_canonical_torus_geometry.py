import pytest
import sympy as sp

from semialg import FilledTorus, Torus


def test_torus_structure_and_dimensions():
    torus = Torus((1, 2, 3), 3, 1)
    assert torus.dimension() == 2
    assert torus.ambient_dimension() == 3
    assert torus.inner_radius == 2
    assert torus.outer_radius == 4
    assert FilledTorus((1, 2, 3), 3, 1).dimension() == 3


def test_torus_formula_hits_surface_and_rejects_center():
    x, y, z = sp.symbols("x y z", real=True)
    torus = Torus((0, 0, 0), 3, 1)
    formula = torus.as_formula((x, y, z))
    assert formula.subs({x: 4, y: 0, z: 0}) is sp.true
    assert formula.subs({x: 3, y: 0, z: 1}) is sp.true
    assert formula.subs({x: 0, y: 0, z: 0}) is sp.false


def test_filled_torus_formula_contains_tube_interior():
    x, y, z = sp.symbols("x y z", real=True)
    solid = FilledTorus((0, 0, 0), 3, 1)
    formula = solid.as_formula((x, y, z))
    assert formula.subs({x: 3, y: 0, z: 0}) is sp.true
    assert formula.subs({x: 0, y: 0, z: 0}) is sp.false
    assert formula.subs({x: 5, y: 0, z: 0}) is sp.false


def test_horn_torus_is_supported():
    assert Torus((0, 0, 0), 1, 1).inner_radius == 0


def test_spindle_torus_is_rejected():
    with pytest.raises(ValueError):
        Torus((0, 0, 0), 1, 2)


def test_symbolic_radii_keep_validity_conditions():
    R, r = sp.symbols("R r", positive=True)
    torus = Torus((0, 0, 0), R, r)
    assert torus.construction_conditions == (R >= r,)
