"""Native quotient arithmetic through parameter-dependent semialg towers."""

import pytest
import sympy as sp
from algroots import QuotientAlgebra, groebner_basis
from algroots.native_f4 import compiled_available

from semialg.algebraic_function_fields import (
    MonogenicFunctionField,
    RationalFunctionField,
    to_coefficient_domain,
)

x, y, t, a, b = sp.symbols("x y t a b")


@pytest.mark.parametrize("backend", ["sympy", "native"])
def test_tower_domain(backend):
    if backend == "native" and not compiled_available():
        pytest.skip("native F4 is not installed")
    F = RationalFunctionField((t,))
    K = MonogenicFunctionField(F, a, (-t, 0, 1))
    L = MonogenicFunctionField(K, b, (-K.alpha - 1, K.zero, K.one))
    domain = to_coefficient_domain(L)
    source = [x * x - b, y - x]
    basis = groebner_basis(source, (x, y), domain=domain, backend=backend)
    assert set(basis.polys) == {
        sp.Poly(y * y - b, x, y, domain=domain),
        sp.Poly(x - y, x, y, domain=domain),
    }
    quotient = QuotientAlgebra.from_polynomials(
        source, (x, y), domain=domain, groebner_backend=backend
    )
    assert quotient.dimension == 2
    assert quotient.normal_form(x * x - b) == 0
    assert domain.convert(domain.to_sympy(domain.convert(b * b - a - 1))) == domain.zero


def test_reducible_bridge():
    F = RationalFunctionField((t,))
    K = MonogenicFunctionField(F, a, (-t * t, 0, 1))
    with pytest.raises(ValueError, match="reducible"):
        to_coefficient_domain(K)
