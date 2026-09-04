from __future__ import annotations

import sympy as sp

from semialg import SemialgebraicContext, is_satisfiable
from semialg.ec.selection import choose_designated_ec
from semialg.incidence import analyze_incidence, sparse_variable_order
from semialg.parameter_stratification import exceptional_parameter_analysis
from semialg.preprocess.algebraic import normalize_polynomial_atoms
from semialg.presolve import presolve_semialgebraic


def test_incidence_analysis_finds_independent_components_and_sparse_order() -> None:
    x, y, z, w = sp.symbols("x y z w", real=True)
    polys = (x**2 + y, z**2 + w)
    analysis = analyze_incidence(polys, (x, y, z, w))
    assert {frozenset(c.variables) for c in analysis.components} == {
        frozenset((x, y)),
        frozenset((z, w)),
    }
    order = sparse_variable_order(polys, (x, y, z, w))
    assert set(order) == {x, y, z, w}


def test_satisfiability_uses_independent_incidence_components() -> None:
    x, y = sp.symbols("x y", real=True)
    result = is_satisfiable(sp.And(x**2 <= 1, y**2 >= 4), (x, y), return_result=True)
    assert result.satisfiable is True
    assert result.method.startswith("incidence_decomposition[")
    assert result.witness is not None
    assert abs(result.witness[x]) <= 1
    assert result.witness[y] ** 2 >= 4


def test_projection_aware_ec_selection_prefers_low_burden_main_variable() -> None:
    x, y = sp.symbols("x y", real=True)
    dense_high_degree = y**4 + x * y**3 + x**2 * y**2 + x**3 * y + 1
    sparse_linear = y + x**5
    choice = choose_designated_ec(
        (dense_high_degree, sparse_linear), policy="auto", variables=(x, y)
    )
    assert choice.expr == sparse_linear
    assert choice.ranked[0].projection_burden <= choice.ranked[1].projection_burden


def test_exceptional_parameter_analysis_tracks_degree_and_discriminant_strata() -> None:
    x, a, b, c = sp.symbols("x a b c", real=True)
    analysis = exceptional_parameter_analysis(sp.Eq(a * x**2 + b * x + c, 0), (x,), (a, b, c))
    by_source = {cause.source: cause.polynomial for cause in analysis.causes}
    assert sp.expand(by_source["degree_drop"] - a) == 0
    assert sp.expand(by_source["discriminant"] ** 2 - (b**2 - 4 * a * c) ** 2) == 0
    assert a in analysis.exceptional_condition.free_symbols


def test_parameterized_quadratic_degree_drop_is_not_misclassified() -> None:
    x, a, b, c = sp.symbols("x a b c", real=True)
    ctx = SemialgebraicContext(sp.Eq(a * x**2 + b * x + c, 0), (x,))
    analysis = exceptional_parameter_analysis(ctx.formula, (x,), (a, b, c))
    assert any(cause.source == "degree_drop" for cause in analysis.causes)
    # The constant nonzero degree-drop fiber has no real root.
    assert is_satisfiable(sp.Eq(0 * x**2 + 0 * x + 1, 0), (x,)) is False


def test_semialgebraic_context_reuses_structural_and_exact_work() -> None:
    x, y = sp.symbols("x y", real=True)
    ctx = SemialgebraicContext(sp.And(x**2 + y**2 <= 1, x >= 0), (x, y))
    first = ctx.complete_cad()
    second = ctx.complete_cad()
    assert first is second
    assert set(ctx.variable_order) == {x, y}
    assert ctx.incidence.components
    stats = ctx.stats()
    assert stats["structural_cache_entries"] >= 4


def test_factor_normalization_and_affine_substitution_preserve_semantics() -> None:
    x, y = sp.symbols("x y", real=True)
    formula = sp.And(sp.Eq(6 * (x + y) ** 2, 0), -4 * x <= 0)
    normalized = normalize_polynomial_atoms(formula)
    assert sp.Eq(x + y, 0) in normalized.args
    assert (x >= 0) in normalized.args
    result = presolve_semialgebraic(formula, (x, y), eliminate=(x,))
    assert result.substitutions == ((x, -y),)
    assert sp.simplify_logic(sp.Equivalent(result.formula, y <= 0)) is sp.true
