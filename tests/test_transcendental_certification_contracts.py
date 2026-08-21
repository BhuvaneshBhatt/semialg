from __future__ import annotations

import pytest
import sympy as sp

from semialg.errors import DimensionMismatchError
from semialg.solve.transcendental.roots import isolate_univar_roots
from semialg.solve.transcendental.system_roots import (
    SearchBox,
    _certify_point,
    orchestrate_trans_search,
)


def test_numeric_residual_is_not_called_a_certificate():
    x = sp.Symbol("x", real=True)
    point = _certify_point((x - sp.sqrt(2),), (x,), (sp.Float("1.41421356237"),))
    assert not point.certified


def test_exact_symbolic_residual_is_certified():
    x = sp.Symbol("x", real=True)
    point = _certify_point((x**2 - 2,), (x,), (sp.sqrt(2),))
    assert point.certified


def test_polynomial_real_root_fallback_does_not_return_complex_roots(monkeypatch):
    x = sp.Symbol("x", real=True)

    def unavailable(*args, **kwargs):
        raise NotImplementedError

    monkeypatch.setattr(sp, "solveset", unavailable)
    result = isolate_univar_roots(x**2 + 1, x, domain=sp.S.Reals)

    assert result.complete
    assert result.roots == ()


def test_finite_transcendental_nonlinsolve_result_is_not_assumed_complete(monkeypatch):
    x, y = sp.symbols("x y", real=True)

    monkeypatch.setattr(
        sp, "nonlinsolve", lambda eqs, vars_: sp.FiniteSet((sp.Integer(0), sp.Integer(0)))
    )
    result = orchestrate_trans_search(sp.And(sp.Eq(sp.sin(x), 0), sp.Eq(y, 0)), (x, y))

    assert result.points
    assert not result.complete
    assert not result.completeness_certificate.complete


def test_root_isolation_rejects_inequality_boolean_input():
    x = sp.Symbol("x", real=True)
    import pytest

    with pytest.raises(TypeError, match="requires an equality or scalar residual"):
        isolate_univar_roots(sp.sin(x) > 0, x)


def test_periodic_inequality_reconstruction_precedes_algebraization():
    from semialg.solve.transcendental import build_trans_state, reduce_trans_problem

    x = sp.Symbol("x", real=True)
    result = reduce_trans_problem(build_trans_state(sp.sin(x) > 0, (x,)))
    assert result.method == "periodic_interval_reconstruction"
    assert result.result_semantics == "periodic_window_approximation"
    assert result.validity_window is not None
    assert not result.complete
    # Regression: the old bug returned a finite set of fake point roots.
    assert not isinstance(result.formula, sp.Equality)


def test_search_box_and_grid_validation():
    x, y = sp.symbols("x y", real=True)
    with pytest.raises(DimensionMismatchError):
        SearchBox((0.0,), (1.0, 2.0))
    with pytest.raises(ValueError, match="lower bound exceeds"):
        SearchBox((2.0,), (1.0,))
    with pytest.raises(ValueError, match="finite"):
        SearchBox((float("-inf"),), (1.0,))
    with pytest.raises(DimensionMismatchError):
        orchestrate_trans_search(
            sp.And(sp.Eq(x, 0), sp.Eq(y, 0)),
            (x, y),
            search_box=SearchBox((-1.0,), (1.0,)),
        )
    with pytest.raises(ValueError, match="at least 1"):
        orchestrate_trans_search(sp.And(sp.Eq(x, 0), sp.Eq(y, 0)), (x, y), grid_points_per_axis=0)
    with pytest.raises(TypeError, match="integer"):
        orchestrate_trans_search(sp.And(sp.Eq(x, 0), sp.Eq(y, 0)), (x, y), grid_points_per_axis=2.5)


def test_result_semantics_distinguish_global_and_window_coverage():
    from semialg.solve.transcendental import ResultSemantics

    assert ResultSemantics.EXACT.is_exact
    assert ResultSemantics.SUBSET.is_subset
    assert ResultSemantics.SUPERSET.is_superset
    assert ResultSemantics.WINDOW_SUBSET.is_subset
    assert ResultSemantics.WINDOW_SUBSET.is_window_scoped
    assert ResultSemantics.WINDOW_SUPERSET.is_superset
    assert ResultSemantics.WINDOW_SUPERSET.is_window_scoped
    assert ResultSemantics.PERIODIC_WINDOW_APPROX.is_window_scoped
    assert not ResultSemantics.WINDOW_APPROXIMATION.is_subset
    assert not ResultSemantics.WINDOW_APPROXIMATION.is_superset
