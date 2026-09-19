"""Generated exact invariants for algebraic decomposition and modular algebra."""

from __future__ import annotations

import pytest
import sympy as sp

import semialg
from semialg.algebraic.groebner_utils import compute_groebner_basis
from semialg.algebraic_decomposition import associated_primes, primary_decomposition

pytestmark = pytest.mark.slow

hypothesis = pytest.importorskip("hypothesis")
assume = hypothesis.assume
given = hypothesis.given
settings = hypothesis.settings
st = pytest.importorskip("hypothesis.strategies")

_SMALL = st.integers(min_value=-3, max_value=3)


@settings(max_examples=16, deadline=None)
@given(a=_SMALL, b=_SMALL, c=_SMALL, d=_SMALL, e1=st.integers(1, 3), e2=st.integers(1, 3))
def test_generated_principal_primary_reconstruction(a, b, c, d, e1, e2):
    x, y = sp.symbols("x y")
    p = x + a * y + b
    q = x + c * y + d
    assume(sp.expand(p - q) != 0)
    source = sp.expand(p**e1 * q**e2)

    result = primary_decomposition((source,), (x, y))
    assert result.complete
    assert semialg.replay_certificate(result).verified
    associated = associated_primes((source,), (x, y))
    assert associated.complete
    assert len(associated.primes) == 2


@settings(max_examples=20, deadline=None)
@given(a=_SMALL, b=_SMALL, c=_SMALL, d=_SMALL)
def test_generated_modular_groebner_matches_direct_exact(a, b, c, d):
    x, y = sp.symbols("x y")
    generators = (x**2 + a * x * y + b * y + c, y**2 + d * x + a)
    direct = compute_groebner_basis(
        generators, (x, y), order="grevlex", domain=sp.QQ, modular=False
    )
    modular = compute_groebner_basis(
        generators, (x, y), order="grevlex", domain=sp.QQ, modular=True
    )
    assert tuple(sp.expand(p.as_expr()) for p in modular.polys) == tuple(
        sp.expand(p.as_expr()) for p in direct.polys
    )


@settings(max_examples=16, deadline=None)
@given(
    a=st.integers(-2, 2),
    b=st.integers(-2, 2),
    c=st.integers(-2, 2),
    d=st.integers(-2, 2),
)
def test_primary_decomposition_invariant_under_generator_permutation_and_units(a, b, c, d):
    x, y = sp.symbols("x y")
    generators = (x**2 + a * x * y + b * y, x * y + c * x + d * y)
    assume(any(sp.expand(g) != 0 for g in generators))
    base = primary_decomposition(generators, (x, y))
    transformed = primary_decomposition((-3 * generators[1], 2 * generators[0]), (x, y))
    assert base.complete == transformed.complete
    if base.complete:
        assert semialg.replay_certificate(base).verified
        assert semialg.replay_certificate(transformed).verified
        assert {tuple(map(sp.expand, component.radical)) for component in base.components} == {
            tuple(map(sp.expand, component.radical)) for component in transformed.components
        }


@pytest.mark.parametrize("shift", range(-2, 3))
def test_embedded_primary_decomposition_covariant_under_affine_translation(shift):
    x, y = sp.symbols("x y")
    u = x + y + shift
    source = (u**2, y * u)
    result = primary_decomposition(source, (x, y))
    assert result.complete
    assert result.irredundant
    assert len(result.components) == 2
    assert semialg.replay_certificate(result).verified
    primes = associated_primes(source, (x, y))
    assert primes.complete
    assert len(primes.primes) == 2
