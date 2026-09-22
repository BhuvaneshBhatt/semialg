import pytest
import sympy as sp

from semialg import connected_components

x, y = sp.symbols("x y", real=True)


@pytest.mark.parametrize(
    ("formula", "count"),
    [
        (sp.Or((x + 2) ** 2 + y**2 < 1, (x - 2) ** 2 + y**2 < 1), 2),
        (sp.Or((x + 2) ** 2 + y**2 <= 1, (x - 2) ** 2 + y**2 <= 1), 2),
        (sp.Or((x + 1) ** 2 + y**2 <= 1, (x - 1) ** 2 + y**2 <= 1), 1),
        (sp.Or((x + 1) ** 2 + y**2 < 1, (x - 1) ** 2 + y**2 < 1), 2),
        (sp.Eq(x**2 + y**2, 1), 1),
        (sp.Or(x**2 + y**2 <= 1, sp.And(sp.Eq(x, 3), sp.Eq(y, 0))), 2),
        (sp.Or(sp.Eq(y, 0), sp.And(sp.Eq(x, 0), y >= 0)), 1),
    ],
    ids=[
        "open-disks",
        "closed-disks",
        "touching-closed-disks",
        "touching-open-disks",
        "equality-curve",
        "mixed-dimensional-disconnected",
        "mixed-dimensional-attached",
    ],
)
def test_connected_component_regression_matrix(formula, count):
    assert len(connected_components(formula, (x, y))) == count


@pytest.mark.slow
def test_connected_annulus_stress():
    annulus = sp.And(x**2 + y**2 >= 1, x**2 + y**2 <= 4)
    assert len(connected_components(annulus, (x, y))) == 1


@pytest.mark.slow
def test_punctured_planar_set_stress():
    punctured = sp.Or(sp.Ne(x, 0), sp.Ne(y, 0))
    assert len(connected_components(punctured, (x, y))) == 1
