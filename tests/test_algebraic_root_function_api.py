from __future__ import annotations

import pytest
import sympy as sp

from semialg import AlgebraicRootFunction, CertifiedRootComparison, DelineabilityCertificate
from semialg.cad.bounds import AlgebraicNumberBound

_DEFAULT_BASE_SYMBOL = sp.Symbol("x", real=True)


def _cert(poly, fiber, root_index, *, sample_value, regular=True, base=(_DEFAULT_BASE_SYMBOL,)):
    return DelineabilityCertificate(
        polynomial=poly,
        fiber_variable=fiber,
        root_index=root_index,
        base_variables=tuple(base),
        sign_invariant=True,
        stack_order_verified=True,
        sample_root_verified=True,
        sample_root_value=sample_value,
        regular_section_verified=regular,
    )


def test_certified_root_comparison_verify_contract():
    assert CertifiedRootComparison(-1, True).verify()
    assert not CertifiedRootComparison(None, False).verify()


def test_delineability_regularity_api_is_cell_wide_certificate():
    x, y = sp.symbols("x y", real=True)
    cert = _cert(y**2 - x, y, 1, sample_value=1, regular=True, base=(x,))
    assert cert.verify()
    assert cert.regular
    assert cert.verify_regularity()


def test_root_function_certified_comparison_on_same_stack():
    x, y = sp.symbols("x y", real=True)
    poly = y**2 - x
    lower = AlgebraicRootFunction(
        poly, y, 0, (x,), (1,), _cert(poly, y, 0, sample_value=-1, base=(x,))
    )
    upper = AlgebraicRootFunction(
        poly, y, 1, (x,), (1,), _cert(poly, y, 1, sample_value=1, base=(x,))
    )
    result = lower.compare_certified(upper)
    assert result.verify()
    assert result.relation == -1
    assert result.scope == "base-cell"


def test_root_function_comparison_declines_unrelated_global_roots():
    x, y = sp.symbols("x y", real=True)
    a = AlgebraicRootFunction(
        y**2 - x, y, 1, (x,), certificate=_cert(y**2 - x, y, 1, sample_value=1, base=(x,))
    )
    b = AlgebraicRootFunction(
        y**2 - x - 1,
        y,
        1,
        (x,),
        certificate=_cert(y**2 - x - 1, y, 1, sample_value=sp.sqrt(2), base=(x,)),
    )
    result = a.compare_certified(b)
    assert not result.verify()
    assert result.relation is None


def test_pointwise_root_comparison_is_exact():
    x, y = sp.symbols("x y", real=True)
    root = AlgebraicRootFunction(y**2 - x, y, 1, (x,))
    result = root.compare_certified(2, base_point={x: 4})
    assert result.verify()
    assert result.relation == 0
    assert result.scope == "point"


def test_full_specialization_returns_algebraic_number_bound():
    x, y = sp.symbols("x y", real=True)
    root = AlgebraicRootFunction(y**2 - x, y, 1, (x,), closed=True)
    specialized = root.specialize({x: 4})
    assert isinstance(specialized, AlgebraicNumberBound)
    assert specialized.closed
    assert sp.simplify(specialized.as_expr() - 2) == 0


def test_partial_specialization_preserves_root_identity_but_drops_certificate():
    x, z, y = sp.symbols("x z y", real=True)
    poly = y**2 - x - z
    cert = _cert(poly, y, 1, sample_value=sp.sqrt(2), base=(x, z))
    root = AlgebraicRootFunction(poly, y, 1, (x, z), certificate=cert)
    specialized = root.specialize({x: 1})
    assert isinstance(specialized, AlgebraicRootFunction)
    assert specialized.base_variables == (z,)
    assert sp.expand(specialized.polynomial - (y**2 - z - 1)) == 0
    assert specialized.certificate is None


def test_certified_implicit_derivative_and_fiber_derivative():
    x, y = sp.symbols("x y", real=True)
    poly = y**2 - x
    cert = _cert(poly, y, 1, sample_value=1, regular=True, base=(x,))
    root = AlgebraicRootFunction(poly, y, 1, (x,), certificate=cert)
    assert root.is_regular
    r = root.as_expr()
    assert sp.simplify(root.fiber_derivative_expr() - 2 * r) == 0
    assert sp.simplify(2 * r * root.derivative_expr(x) - 1) == 0
    assert root.derivative_expr(sp.Symbol("z", real=True)) == 0


def test_derivative_requires_regular_certificate_by_default():
    x, y = sp.symbols("x y", real=True)
    root = AlgebraicRootFunction(y**2 - x, y, 1, (x,))
    assert not root.is_regular
    with pytest.raises(ValueError, match="certified regularity"):
        root.derivative_expr(x)
    r = root.as_expr()
    assert sp.simplify(2 * r * root.derivative_expr(x, require_regular=False) - 1) == 0
