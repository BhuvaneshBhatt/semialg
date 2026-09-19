from dataclasses import replace

import sympy as sp

from semialg.algebraic import (
    modular_groebner_basis_fraction_field,
    verify_modular_fraction_field_groebner_certificate,
)
from semialg.algebraic.groebner_utils import (
    _prefer_modular_fraction_field,
    compute_groebner_basis,
)
from semialg.algebraic.gtz import certified_independent_localization
from semialg.algebraic_decomposition import _splitter_search_cost
from semialg.algebraic_function_fields import (
    MonogenicFunctionField,
    RationalFunctionField,
    compress_primitive_element,
    verify_primitive_element_compression,
)


def test_modular_fraction_field_groebner_reconstructs_denominators():
    u, x, y = sp.symbols("u x y")
    generators = (x / (u + 1) + y, y**2 - u / (u + 2))
    result = modular_groebner_basis_fraction_field(generators, (x, y), (u,), max_primes=8)
    assert result is not None
    assert result.certified
    assert len(result.certificate.primes) >= 2
    assert verify_modular_fraction_field_groebner_certificate(result.certificate)

    domain = sp.QQ.frac_field(u)
    direct = sp.groebner(generators, x, y, order="grevlex", domain=domain)
    direct_basis = tuple(sp.cancel(poly.as_expr()) for poly in direct.polys)
    assert result.basis == direct_basis


def test_fraction_field_dispatcher_uses_certified_candidate_semantics():
    u, x, y = sp.symbols("u x y")
    domain = sp.QQ.frac_field(u)
    generators = ((u + 1) * x + y, x * y - u)
    accelerated = compute_groebner_basis(
        generators, (x, y), order="grevlex", domain=domain, modular=None
    )
    direct = compute_groebner_basis(
        generators, (x, y), order="grevlex", domain=domain, modular=False
    )
    assert tuple(sp.cancel(p.as_expr()) for p in accelerated.polys) == tuple(
        sp.cancel(p.as_expr()) for p in direct.polys
    )


def test_localization_planner_prefers_lower_parameter_degree():
    x, y = sp.symbols("x y")
    result = certified_independent_localization((x * y**3 - 1,), (x, y))
    assert result.independent_variables == (x,)
    assert result.dependent_variables == (y,)
    assert result.quotient_dimension == 3


def test_regular_chain_splitter_cost_prefers_small_sparse_condition():
    x, y, z = sp.symbols("x y z")
    assert _splitter_search_cost(x + 1, (x, y, z)) < _splitter_search_cost(
        x**4 + x * y + y**3 + z, (x, y, z)
    )


def test_reconstructed_primitive_guard_replays_degree_four_tower():
    a, b = sp.symbols("a b")
    base = RationalFunctionField(tuple())
    lower = MonogenicFunctionField(base, a, (-2, 0, 1))
    top = MonogenicFunctionField(lower, b, (lower.convert(-3), lower.zero, lower.one))
    compression = compress_primitive_element(top)
    assert compression.complete
    assert compression.field.degree == 4
    assert len(compression.steps) == 1
    assert compression.steps[0].bad_parameter_guard != 0
    assert verify_primitive_element_compression(compression)


def test_fraction_field_cost_gate_skips_tiny_and_accepts_larger_systems():
    u, x, y, z = sp.symbols("u x y z")
    domain = sp.QQ.frac_field(u)
    assert not _prefer_modular_fraction_field((x - u,), (x,), domain)
    larger = (
        x**3 + y**2 + z + u,
        y**3 + z**2 + x + u**2,
        z**3 + x**2 + y + u**3,
    )
    assert _prefer_modular_fraction_field(larger, (x, y, z), domain)


def test_fraction_field_modular_certificate_rejects_membership_tamper():
    u, x, y = sp.symbols("u x y")
    result = modular_groebner_basis_fraction_field(
        ((u + 1) * x + y, x * y - u), (x, y), (u,), max_primes=8
    )
    assert result is not None
    representations = list(result.certificate.membership_representations)
    first = list(representations[0])
    first[0] = sp.cancel(first[0] + 1)
    representations[0] = tuple(first)
    tampered = replace(result.certificate, membership_representations=tuple(representations))
    assert not verify_modular_fraction_field_groebner_certificate(tampered)
