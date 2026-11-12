"""Production primary decomposition and RUR routing over coefficient towers."""

import pickle
from dataclasses import replace

import pytest
import sympy as sp
from algroots.native_f4 import compiled_available

from semialg.algebraic.gtz_zero_dim import (
    verify_zero_dim_primary_certificate,
    zero_dimensional_primary_decomposition,
)
from semialg.algebraic.rational_univariate import compute_rur
from semialg.algebraic_decomposition import (
    primary_decomposition,
    verify_primary_decomposition_certificate,
)
from semialg.algebraic_function_fields import MonogenicFunctionField, RationalFunctionField

x, y, t, a, b = sp.symbols("x y t a b")


def coefficient_fields():
    F = RationalFunctionField((t,))
    K = MonogenicFunctionField(F, a, (-t, 0, 1))
    L = MonogenicFunctionField(K, b, (-K.alpha - 1, K.zero, K.one))
    return [(F, t), (K, a), (L, b), (sp.QQ.algebraic_field(sp.sqrt(2)), sp.sqrt(2))]


@pytest.mark.parametrize("field,root", coefficient_fields())
@pytest.mark.parametrize("backend", ["sympy", "native"])
def test_primary_field_route(field, root, backend):
    if backend == "native" and not compiled_available():
        pytest.skip("native F4 is not installed")
    result = primary_decomposition(
        [(x - root) ** 2 * (x + root) ** 3, y * y],
        (x, y),
        coefficient_field=field,
        groebner_backend=backend,
    )
    assert result.complete and result.irredundant
    assert sorted(component.degree for component in result.components) == [4, 6]
    assert verify_primary_decomposition_certificate(result.certificate)
    restored = pickle.loads(pickle.dumps(result.certificate))
    assert verify_primary_decomposition_certificate(restored)


@pytest.mark.parametrize("field,root", coefficient_fields()[:3])
def test_rur_field_route(field, root):
    result = compute_rur([x * x - root, y - x], (x, y), domain=field)
    assert result.dimension == 2 and result.solution_count == 2
    coordinates = result.normalized_coordinate_polynomials()
    D = result.defining_polynomial.domain
    assert (
        (coordinates[0] ** 2 - sp.Poly(root, result.parameter, domain=D))
        .rem(result.defining_polynomial)
        .is_zero
    )
    assert coordinates[0] == coordinates[1]
    assert field.as_domain() is field.as_domain()


def test_localized_primary_entry():
    F, K = coefficient_fields()[0][0], coefficient_fields()[1][0]
    for field in (F, K):
        result = zero_dimensional_primary_decomposition(
            [x * x - t, y], (x, y), coefficient_field=field
        )
        assert verify_zero_dim_primary_certificate(result.certificate)
        assert sum(c.degree for c in result.components) == 2


def test_primary_mutations():
    field, root = coefficient_fields()[1]
    result = primary_decomposition(
        [(x - root) ** 2 * (x + root) ** 3, y * y], (x, y), coefficient_field=field
    )
    certificate = result.certificate
    first = certificate.components[0]
    bad = (
        replace(certificate, complete=False),
        replace(certificate, irredundant=False),
        replace(certificate, coefficient_field=sp.QQ),
        replace(certificate, variables=(y, x)),
        replace(certificate, source_equations=(x, y)),
        replace(certificate, components=()),
        replace(
            certificate,
            components=(replace(first, degree=first.degree + 1), *certificate.components[1:]),
        ),
        replace(
            certificate, components=(replace(first, primary=False), *certificate.components[1:])
        ),
        replace(
            certificate, basis_certificate=replace(certificate.basis_certificate, memberships=())
        ),
    )
    assert all(not verify_primary_decomposition_certificate(proof) for proof in bad)


def test_finite_coefficient_field():
    result = primary_decomposition([x * x - x, y * y], (x, y), coefficient_field=sp.GF(2))
    assert sorted(c.degree for c in result.components) == [2, 2]
    assert verify_primary_decomposition_certificate(result.certificate)


def test_field_domain_cache_detaches_inputs():
    parameters = [t]
    base = RationalFunctionField(parameters)
    parameters[0] = a
    modulus = [-t, 0, 1]
    field = MonogenicFunctionField(base, a, modulus)
    domain = field.as_domain()
    modulus[0] = -t - 1
    assert base.parameters == (t,)
    assert field.modulus == (-t, 0, 1)
    assert domain.convert(a * a - t) == domain.zero
    assert field.as_domain() is domain
