import pytest
import sympy as sp

from semialg.algebraic_decomposition import (
    EquidimensionalDecomposition,
    equidimensional_decomposition,
)
from semialg.region_analysis import (
    _region_singular_locus_result,
    region_singular_locus_result,
)


def _incomplete_decomposition(equations, variables, **kwargs):
    result = equidimensional_decomposition(equations, variables, **kwargs)
    return EquidimensionalDecomposition(
        variables=result.variables,
        pieces=result.pieces,
        complete=False,
    )


def test_incomplete_equality_decomposition_is_explicit_not_global_rank_fallback():
    x, y, z = sp.symbols("x y z", real=True)
    equations = (x * (x - 1), x * (z**2 - y**3))
    region = sp.And(*(sp.Eq(equation, 0) for equation in equations))
    result = _region_singular_locus_result(
        region, (x, y, z), decomposition_provider=_incomplete_decomposition
    )

    assert not result.complete
    assert result.formula is None
    assert result.known_singular_formula is sp.false
    assert "global-rank fallback is intentionally disabled" in result.diagnostics[0]
    with pytest.raises(NotImplementedError, match="global-rank fallback"):
        result.require_complete()


def test_incomplete_boundary_component_decomposition_propagates_unknown():
    x, y, z = sp.symbols("x y z", real=True)
    region = sp.And(sp.Eq(x * y, 0), z >= 0)
    result = _region_singular_locus_result(
        region, (x, y, z), decomposition_provider=_incomplete_decomposition
    )

    assert not result.complete
    assert result.formula is None
    assert isinstance(result.known_singular_formula, sp.Basic)
    assert any("equidimensional decomposition" in message for message in result.diagnostics)
    with pytest.raises(NotImplementedError):
        result.require_complete()


def test_complete_result_exposes_exact_formula_and_certified_lower_bound():
    x, y = sp.symbols("x y", real=True)
    result = region_singular_locus_result(sp.Eq(x * y, 0), (x, y))

    assert result.complete
    assert result.formula == result.known_singular_formula
    assert sp.simplify(result.formula.subs({x: 0, y: 0})) is sp.true
    assert sp.simplify(result.formula.subs({x: 1, y: 0})) is sp.false
