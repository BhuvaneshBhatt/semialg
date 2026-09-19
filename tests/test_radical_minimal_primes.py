import sympy as sp

from semialg.algebraic_decomposition import certified_radical_minimal_prime_decomposition


def test_coordinate_minimal_primes_are_certified():
    x, y = sp.symbols("x y")
    result = certified_radical_minimal_prime_decomposition((x * y,), (x, y))
    assert result.radical_complete
    assert result.minimal_primes_complete
    assert len(result.components) == 2
    assert all(c.radical and c.prime for c in result.components)


def test_embedded_monomial_component_is_removed():
    x, y = sp.symbols("x y")
    result = certified_radical_minimal_prime_decomposition((x**2, x * y), (x, y))
    assert result.radical_complete
    assert result.minimal_primes_complete
    assert len(result.components) == 1
    assert result.components[0].prime


def test_irreducible_hypersurface_prime():
    x, y = sp.symbols("x y")
    result = certified_radical_minimal_prime_decomposition((x**2 + y**2 - 1,), (x, y))
    assert result.radical_complete
    assert result.minimal_primes_complete
    assert result.components[0].prime_method == "irreducible_hypersurface"


def test_monomial_curve_triangular_graph_prime():
    x, y, z = sp.symbols("x y z")
    result = certified_radical_minimal_prime_decomposition((x**2 - y, x * y - z), (x, y, z))
    assert result.radical_complete
    assert result.minimal_primes_complete
    assert len(result.components) == 1
    assert result.components[0].prime_method == "triangular_fraction_field_graph_prime"


def test_radical_but_unproved_prime_is_not_promoted():
    x, y = sp.symbols("x y")
    # squarefree reducible hypersurface is split first, so choose a regular-chain
    # branch whose primality is intentionally outside the conservative prover.
    result = certified_radical_minimal_prime_decomposition((y**2 - x**3 - x,), (x, y))
    assert result.radical_complete
    # This curve may be recognized as irreducible; the invariant we care about
    # is that minimal-prime completeness never exceeds exact prime certification.
    assert result.minimal_primes_complete == all(c.prime for c in result.components)


def test_successive_fraction_field_splits_after_linear_quotient():
    x, y, z = sp.symbols("x y z")
    result = certified_radical_minimal_prime_decomposition((y - x**2, z**2 - y), (x, y, z))
    assert result.radical_complete
    assert result.minimal_primes_complete
    assert len(result.components) == 2
    equations = [set(component.equations) for component in result.components]
    assert any(x - z in branch or -x + z in branch for branch in equations)
    assert any(x + z in branch or -x - z in branch for branch in equations)
    assert all(component.prime for component in result.components)


def test_successive_fraction_field_prime_with_nonunit_linear_initial():
    x, y, z = sp.symbols("x y z")
    result = certified_radical_minimal_prime_decomposition((x * y - 1, z**2 - y), (x, y, z))
    assert result.radical_complete
    assert result.minimal_primes_complete
    assert len(result.components) == 1
    assert result.components[0].prime_method == "successive_fraction_field_prime"


def test_second_nonlinear_algebraic_extension_is_certified_by_tower():
    x, y, z = sp.symbols("x y z")
    # The second nonlinear stage is decided in QQ(x)(y), y**2=x.
    result = certified_radical_minimal_prime_decomposition((y**2 - x, z**2 - y), (x, y, z))
    assert result.radical_complete
    assert result.minimal_primes_complete
    assert result.components[0].prime_method == "algebraic_function_field_tower_prime"
    assert len(result.components) == 1
    assert result.components[0].radical
    assert result.components[0].prime


def test_second_nonlinear_tower_detects_reducible_stage():
    x, y, z = sp.symbols("x y z")
    result = certified_radical_minimal_prime_decomposition((y**2 - x, z**2 - x), (x, y, z))
    assert result.radical_complete
    assert result.minimal_primes_complete
    assert len(result.components) == 2
    for component in result.components:
        assert component.prime
        assert component.prime_method == "algebraic_function_field_tower_prime"
    equations = [set(map(sp.expand, component.equations)) for component in result.components]
    assert any(sp.expand(y - z) in eqs or sp.expand(-y + z) in eqs for eqs in equations)
    assert any(sp.expand(y + z) in eqs or sp.expand(-y - z) in eqs for eqs in equations)


def test_general_radical_ideal_api_reconstructs_reduced_generators():
    from semialg.algebraic_decomposition import radical_ideal

    x, y = sp.symbols("x y")
    result = radical_ideal((x**2, x * y), (x, y))
    assert result.complete
    assert result.generators == (x,)
    assert result.decomposition.complete
    assert result.decomposition.reconstruction_complete


def test_triangular_primality_certificate_records_successive_field_stages():
    from semialg.algebraic_decomposition import certify_triangular_primality

    x, y, z = sp.symbols("x y z")
    certificate = certify_triangular_primality((y**2 - x, z**2 - y), (x, y, z))
    assert certificate.complete
    assert certificate.prime
    assert [stage.leader for stage in certificate.stages] == [y, z]
    assert all(stage.reconstruction_verified for stage in certificate.stages)
    assert all(stage.factor_degrees == (2,) for stage in certificate.stages)


def test_triangular_primality_certificate_exposes_reducible_stage():
    from semialg.algebraic_decomposition import certify_triangular_primality

    x, y, z = sp.symbols("x y z")
    certificate = certify_triangular_primality((y**2 - x, z**2 - x), (x, y, z))
    assert certificate.complete
    assert not certificate.prime
    assert len(certificate.splitters) == 2
