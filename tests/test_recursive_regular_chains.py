import sympy as sp

from semialg.algebraic_decomposition import (
    certified_radical_minimal_prime_decomposition,
    recursive_regular_chain_decomposition,
)


def test_regular_chain_split_conserves_degree():
    x, y = sp.symbols("x y")
    result = recursive_regular_chain_decomposition((x * y,), (x, y))
    assert result.complete
    assert result.degree_complete
    assert result.source_dimension == 1
    assert result.source_degree == 2
    assert len(result.components) == 2
    assert {component.degree for component in result.components} == {1}
    assert all(component.squarefree_regular for component in result.components)


def test_recursive_regular_chain_handles_quotient_field_factor_split():
    x, y, z = sp.symbols("x y z")
    result = recursive_regular_chain_decomposition((y**2 - x, z**2 - x), (x, y, z))
    assert result.complete
    assert result.reconstruction_complete
    assert result.degree_complete
    assert len(result.components) == 1
    assert result.components[0].dimension == 1
    assert result.components[0].degree == result.reduced_degree


def test_minimal_prime_result_carries_hilbert_degree_certificate():
    x, y = sp.symbols("x y")
    result = certified_radical_minimal_prime_decomposition((x * y,), (x, y))
    assert result.radical_complete
    assert result.minimal_primes_complete
    assert result.degree_complete
    assert sorted(component.degree for component in result.components) == [1, 1]


def test_regular_chain_degree_uses_reduced_union():
    x, y = sp.symbols("x y")
    result = recursive_regular_chain_decomposition((x**2 * y,), (x, y))
    assert result.complete
    assert result.reconstruction_complete
    assert result.source_degree == 3
    assert result.reduced_degree == 2
    assert result.degree_complete
    assert sorted(component.degree for component in result.components) == [1, 1]


def test_recursive_regular_chain_uses_competing_leader_regular_gcd_split():
    x, y, z = sp.symbols("x y z")
    result = recursive_regular_chain_decomposition((x * y + z, y * z), (x, y, z))
    assert result.complete
    assert result.reconstruction_complete
    assert "regular_competing_leader_gcd_split" in result.methods
    equations = {component.equations for component in result.components}
    assert equations == {(x, z), (y, z)}


def test_recursive_regular_chain_splits_zero_divisor_separant_at_failed_prefix():
    x, y, z = sp.symbols("x y z")
    result = recursive_regular_chain_decomposition((x * y - z, x * z - y), (x, y, z))
    assert result.complete
    assert result.reconstruction_complete
    assert "regular_zero_divisor_separant_split" in result.methods
    assert len(result.components) == 2
    assert all(component.squarefree_regular for component in result.components)
