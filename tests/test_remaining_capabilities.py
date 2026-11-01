"""Exact coefficient fields, bounded linear dispatch and public resource controls."""

from unittest.mock import patch

import pytest
import sympy as sp

from semialg import cad, matrix_pd_on, quantifier_eliminate
from semialg._critical_loci import project_kkt_locus
from semialg.algebraic.groebner_utils import algebraic_coefficient_field
from semialg.cache_utils import cache_disabled
from semialg.errors import ResourceLimitError
from semialg.geometry_queries import semialgebraic_preimage
from semialg.optimization_geometry import polynomial_locus_dimension
from semialg.qe.linear import try_linear_elimination


@pytest.mark.parametrize(
    "alpha,beta",
    [(sp.sqrt(2), sp.sqrt(3)), (-sp.sqrt(2), sp.sqrt(2)), (sp.real_root(2, 3), sp.sqrt(5))],
)
def test_critical_number_fields(alpha, beta):
    x, y, multiplier = sp.symbols("x y multiplier", real=True)
    equations = (multiplier - alpha * x, multiplier - beta * y)
    projected = project_kkt_locus(equations, (multiplier,), (x, y))
    field = algebraic_coefficient_field(equations, (multiplier, x, y))
    assert field.is_AlgebraicField
    assert polynomial_locus_dimension(projected, (x, y)) == 1
    basis = sp.groebner(projected, x, y, domain=field)
    assert basis.reduce(alpha * x - beta * y)[1] == 0
    assert polynomial_locus_dimension((*projected, x - 1), (x, y)) == 0
    assert polynomial_locus_dimension((*projected, x, y - 1), (x, y)) == -1


def test_field_embedding_and_parameters():
    x, y, a = sp.symbols("x y a", real=True)
    positive = sp.CRootOf(x**3 - x - 1, 0)
    assert polynomial_locus_dimension((x - positive * y,), (x, y)) == 1
    assert algebraic_coefficient_field((2 * x - y,), (x, y)) == sp.QQ
    assert polynomial_locus_dimension((a * x - y,), (x, y)) is None
    assert polynomial_locus_dimension((sp.pi * x - y,), (x, y)) is None


@pytest.mark.parametrize(
    "prefix,expected", [(("forall", "exists"), True), (("exists", "forall"), False)]
)
def test_linear_alternation(prefix, expected):
    x, y = sp.symbols("x y", real=True)
    with patch(
        "semialg.solve.reduce.try_quadratic_vs_qe", side_effect=AssertionError("unnecessary VS")
    ):
        result = quantifier_eliminate(
            x < y, tuple(zip(prefix, (x, y), strict=True)), return_result=True
        )
    assert result.formula is (sp.true if expected else sp.false)
    assert result.method == "linear-boolean"


def test_linear_boolean_boundaries():
    a, x = sp.symbols("a x", real=True)
    matrix = sp.Or(sp.And(x > 0, x < a), sp.And(x < -1, x > a))
    result = try_linear_elimination(matrix, (("exists", x),), (a, x))
    for value in (-2, -1, 0, 1):
        assert bool(result.formula.subs(a, value)) == (value < -1 or value > 0)
    assert try_linear_elimination(matrix, (("exists", x),), (a, x), max_terms=1) is None
    assert try_linear_elimination(a * x > 1, (("exists", x),), (a, x)) is None


@pytest.mark.parametrize("strategy", ["collins", "lazard-complete", "auto"])
def test_public_cell_budget(strategy):
    x = sp.Symbol("x", real=True)
    with cache_disabled():
        result = cad(x > 0, (x,), strategy=strategy, max_cells=3, return_result=True)
    assert result.status == "complete"
    assert result.as_function()(1) and not result.as_function()(0)
    # Cached and fresh results obey the same bound.
    for _ in range(2):
        assert (
            cad(x > 0, (x,), strategy=strategy, max_cells=2, return_result=True).status == "unknown"
        )
    with pytest.raises(ResourceLimitError):
        cad(x > 0, (x,), max_cells=2)


def test_budget_survives_preprocessing():
    x = sp.Symbol("x", real=True)
    source = sp.sqrt(x) > 1
    limited = cad(source, (x,), max_cells=3, return_result=True)
    assert limited.status == "unknown"
    assert limited.formula == source
    with pytest.raises(ResourceLimitError):
        limited.as_function()
    sufficient = cad(source, (x,), max_cells=100, return_result=True)
    assert sufficient.status == "complete"
    assert sufficient.as_function()(2) and not sufficient.as_function()(1)


def test_public_deadline():
    x = sp.Symbol("x", real=True)
    with patch("semialg.cad_algorithms.resource_limits.perf_counter", side_effect=[0, 2]):
        result = cad(x > 0, (x,), timeout=1, return_result=True)
    assert result.status == "unknown"
    with patch("semialg.cad_algorithms.resource_limits.perf_counter", side_effect=[0, 2]):
        with pytest.raises(ResourceLimitError):
            quantifier_eliminate(x > 0, (("exists", x),), timeout=1)


@pytest.mark.parametrize(
    "option,value",
    [("max_cells", True), ("max_cells", 1.5), ("timeout", float("nan")), ("timeout", True)],
)
def test_budget_validation(option, value):
    x = sp.Symbol("x", real=True)
    with pytest.raises(ValueError):
        cad(x > 0, (x,), **{option: value})


def test_empty_domain_and_column_maps():
    assert matrix_pd_on([[-1]], domain=sp.false) is True
    x, u, v = sp.symbols("x u v", real=True)
    matrix = sp.Matrix([x, 1 / x])
    target = sp.And(u > 0, v < 2)
    column = semialgebraic_preimage(matrix, target, (x,), target_variables=(u, v))
    vector = semialgebraic_preimage(tuple(matrix), target, (x,), target_variables=(u, v))
    assert column == vector
    assert column.subs(x, 0) is sp.false
    assert column.subs(x, 1) is sp.true


def test_number_field_critical_image():
    from semialg import critical_value_image
    from semialg.optimization_active_sets import constant_gradient_rank, pruned_active_subsets

    x, y = sp.symbols("x y", real=True)
    result = critical_value_image(sp.sqrt(2) * x, sp.Eq((y - sp.sqrt(3)) ** 2, 0), (x, y))
    assert result.formula is sp.true
    assert result.component_images == (sp.true,)
    subsets = pruned_active_subsets((x,), (sp.sqrt(2) * x + y,), (x, y))
    assert subsets == ((), (x + sp.sqrt(2) * y / 2,))
    root = sp.CRootOf(x**3 - x - 1, 0)
    assert constant_gradient_rank(((root**3 - root - 1) * x,), (x,)) == 0


@pytest.mark.parametrize(
    "requested",
    ["positive_definite", "positive_semidefinite", "negative_definite", "negative_semidefinite"],
)
def test_empty_matrix_domains(requested):
    from semialg.matrix_analysis import matrix_definiteness

    assert matrix_definiteness([[1]], domain=sp.false, requested=requested) is True
