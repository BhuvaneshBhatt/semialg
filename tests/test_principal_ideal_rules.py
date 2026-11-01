"""UFD identities avoid general elimination while preserving exact ideals."""

from dataclasses import replace

import sympy as sp

from semialg.algebraic.ideal_ops import ideal_intersection_qq, principal_ideal_intersection
from semialg.algebraic_decomposition import (
    primary_decomposition,
    verify_primary_decomposition_certificate,
)
from semialg.cache_control import clear_caches


def test_principal_intersection_lcm():
    x, y = sp.symbols("x y")
    first = (x * (x + y)) ** 2
    second = (x + y) ** 3 * y
    result = ideal_intersection_qq((first,), (second,), (x, y))
    expected = sp.expand(x * x * y * (x + y) ** 3)
    assert result == (expected,)
    assert principal_ideal_intersection((), (x,), (x, y)) == ()
    assert ideal_intersection_qq((1,), (x,), (x, y)) == (x,)


def test_algebraic_principal_intersection():
    x, y = sp.symbols("x y")
    factor = x - sp.sqrt(2) * y
    result = principal_ideal_intersection((factor**2,), (factor * y,), (x, y))
    assert (
        sp.Poly(result[0], x, y, extension=True).monic()
        == sp.Poly(factor**2 * y, x, y, extension=True).monic()
    )


def test_cached_primary_proof_is_bound_to_canonical_ideal():
    x, y = sp.symbols("x y")
    clear_caches(collect=False)
    first = primary_decomposition((x * x, x * y), (x, y))
    second = primary_decomposition((-3 * x * y, 2 * x * x), (x, y))
    assert first is second
    assert {tuple(component.radical) for component in first.components} == {(x,), (x, y)}
    assert verify_primary_decomposition_certificate(first.certificate)
    component = first.certificate.components[0]
    forged = replace(
        first.certificate,
        components=(
            replace(component, degree=component.degree + 1),
            *first.certificate.components[1:],
        ),
    )
    assert not verify_primary_decomposition_certificate(forged)
    clear_caches(collect=False)
    third = primary_decomposition((x * x, x * y), (x, y))
    assert third is not first
    assert verify_primary_decomposition_certificate(third.certificate)
