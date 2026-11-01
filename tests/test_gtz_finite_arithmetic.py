"""Prescribed local lengths and independent nilpotence checks over QQ(a)."""

from dataclasses import replace

import pytest
import sympy as sp

from semialg.algebraic.gtz_primary import (
    _compute_fraction_field_zero_dimensional,
    verify_fraction_field_zero_dim_certificate,
)

x, y, a = sp.symbols("x y a")


@pytest.mark.parametrize(
    "polynomial,lengths,rank",
    [((x - a) ** 2 * (x + a) ** 3, [4, 6], 2), ((x * x - a) ** 2, [8], 2)],
)
def test_prescribed_factors(polynomial, lengths, rank):
    source = (polynomial, y * y)
    result = _compute_fraction_field_zero_dimensional(source, (x, y), (a,))
    certificate = result.certificate
    assert result.complete and verify_fraction_field_zero_dim_certificate(certificate)
    assert len(certificate.standard_monomials) == sum(lengths)
    assert certificate.trace_rank == rank
    separator = next(iter(certificate.characteristic_polynomial.free_symbols - {a}))
    assert (
        sorted(
            sp.Poly(component.irreducible_factor, separator, domain=sp.QQ.frac_field(a)).degree()
            * component.multiplicity
            for component in result.components
        )
        == lengths
    )
    basis = sp.groebner(source, x, y, domain=sp.QQ.frac_field(a))
    for value in certificate.nilradical_generators:
        assert basis.reduce(sp.expand(value ** sum(lengths)))[1] == 0
    assert len(certificate.nilradical_generators) == sum(lengths) - rank
    for changes in (
        {"trace_rank": rank + 1},
        {"standard_monomials": certificate.standard_monomials[:-1]},
        {"nilradical_generators": ()},
        {"characteristic_polynomial": certificate.characteristic_polynomial + 1},
        {"components": ()},
    ):
        assert not verify_fraction_field_zero_dim_certificate(replace(certificate, **changes))


def test_zero_algebra():
    result = _compute_fraction_field_zero_dimensional((1,), (x,), (a,))
    certificate = result.certificate
    assert result.components == () and result.complete
    assert verify_fraction_field_zero_dim_certificate(certificate)
    assert not verify_fraction_field_zero_dim_certificate(replace(certificate, trace_rank=1))
    assert not verify_fraction_field_zero_dim_certificate(replace(certificate, primitive_element=1))
