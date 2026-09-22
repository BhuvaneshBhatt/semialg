from __future__ import annotations

import sympy as sp

from semialg.reconstruct.root_functions import root_of
from semialg.topology.incidence import _truth_condition_at_assignments


def test_incidence_specializes_nested_cylindrical_root_samples_exactly():
    x, y, z = sp.symbols("x y z", real=True)
    assignments = {
        x: sp.sqrt(2),
        y: root_of(y**2 - x, y, 1),
        z: root_of(z**2 - y, z, 1),
    }
    boundary = root_of(z**2 - y, z, 1)

    assert _truth_condition_at_assignments(sp.Eq(z, boundary), assignments)
    assert _truth_condition_at_assignments(z > 0, assignments)


def test_incidence_root_specialization_closes_base_samples_before_isolation():
    x, y, z = sp.symbols("x y z", real=True)
    assignments = {
        x: sp.sqrt(2),
        y: root_of(y**2 - x, y, 1),
        z: 0,
    }
    upper = root_of(z**2 - y, z, 1)

    assert _truth_condition_at_assignments(upper > 0, assignments)


def test_incidence_lower_root_specialization_ignores_higher_sample_coordinates():
    x, y, z = sp.symbols("x y z", real=True)
    assignments = {
        x: -2,
        y: sp.sqrt(5),
        z: root_of(z**2 - y, z, 1),
    }
    lower_level_boundary = root_of(y**2 - x**2 - 1, y, 1)

    assert _truth_condition_at_assignments(sp.Eq(y, lower_level_boundary), assignments)
