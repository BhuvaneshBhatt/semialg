from pathlib import Path

import sympy as sp

from semialg._range_special_cases import _relation_from_scalar
from semialg._zero_testing import certified_pointwise_zero, certified_zero
from semialg.implicit_geometry import VerticalBoundCell2D


def test_scalar_unequality_preserves_unknown_zero_status():
    a = sp.Symbol("a", real=True)
    assert certified_zero(a) is False
    assert certified_pointwise_zero(a) is None
    assert _relation_from_scalar(a, "==") is None
    assert _relation_from_scalar(a, "!=") is None


def test_vertical_cell_does_not_assume_symbolic_interval_has_positive_width():
    x, y, a = sp.symbols("x y a", real=True)
    cell = VerticalBoundCell2D(x, y, (sp.Integer(0), a), ((sp.Integer(0), sp.Integer(1)),))
    assert cell.is_full_dimensional is False


def test_no_simplify_not_equal_zero_proof_decisions_remain():
    source = Path(__file__).parents[1] / "src" / "semialg"
    offenders = []
    for path in source.rglob("*.py"):
        text = path.read_text(encoding="utf-8")
        for line_number, line in enumerate(text.splitlines(), 1):
            if "simplify(" in line and "!= 0" in line:
                offenders.append(f"{path.relative_to(source)}:{line_number}")
    assert offenders == []


def test_circle_coefficient_requires_nonzero_proof():
    from semialg._region_integrate_intrinsic import _circle_radius_squared

    x, y, a = sp.symbols("x y a", real=True)
    assert _circle_radius_squared(sp.Eq(a * x**2 + a * y**2, 1), x, y) is None


def test_shared_constant_sign_boundary_preserves_parameter_uncertainty():
    from semialg._zero_testing import certified_constant_sign

    a = sp.Symbol("a", real=True)
    p = sp.Symbol("p", positive=True)
    assert certified_constant_sign(sp.Integer(0)) == 0
    assert certified_constant_sign(sp.Rational(-3, 7)) == -1
    assert certified_constant_sign(sp.sqrt(2)) == 1
    assert certified_constant_sign(a) is None
    assert certified_constant_sign(p) == 1
