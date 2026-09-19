"""Failure categories stay distinct from valid negative mathematical answers."""

from __future__ import annotations

import pytest
import sympy as sp

from semialg import (
    IntervalRegion,
    ResourceLimitError,
    cad,
    classify_real_roots,
    equivalent,
    implies,
    integrate_over_region,
    is_satisfiable,
    is_tautology,
    matrix_pd_on,
    matrix_psd_on,
    sample_points,
    semialgebraic_image,
    semialgebraic_maximize,
    semialgebraic_minimize,
    semialgebraic_preimage,
    semialgebraic_projection,
    translate,
)

X, Y, U = sp.symbols("x y u", real=True)


@pytest.mark.parametrize("structured", [False, True], ids=["direct", "structured"])
@pytest.mark.parametrize("domain", ["integers", "complexes"])
def test_decision_family_rejects_unsupported_domains(domain, structured):
    options = {"domain": domain, "return_result": structured}
    with pytest.raises(NotImplementedError, match="real domain"):
        is_satisfiable(X >= 0, (X,), **options)
    with pytest.raises(NotImplementedError, match="real domain"):
        is_tautology(X >= 0, (X,), **options)
    with pytest.raises(NotImplementedError, match="real domain"):
        implies(X > 1, X > 0, (X,), **options)
    with pytest.raises(NotImplementedError, match="real domain"):
        equivalent(X > 0, X >= 0, (X,), **options)
    assert is_satisfiable(sp.Eq(X**2, -1), (X,)) is False


@pytest.mark.parametrize("structured", [False, True], ids=["direct", "structured"])
@pytest.mark.parametrize(
    "bounds, message",
    [
        pytest.param({X: (1, 0)}, "lower bound exceeds", id="reversed"),
        pytest.param({Y: (0, 1)}, "not in the variable list", id="undeclared"),
        pytest.param({X: (0, 1, 2)}, "pair", id="malformed"),
        pytest.param(((X, 0, 1), (X, 0, 2)), "duplicate bound", id="duplicate"),
    ],
)
def test_integral_rejects_invalid_bounds_without_returning_zero(bounds, message, structured):
    with pytest.raises(ValueError, match=message):
        integrate_over_region(1, X >= 0, (X,), bounds=bounds, return_result=structured)
    assert integrate_over_region(1, sp.false, (X,)) == 0


@pytest.mark.parametrize("structured", [False, True], ids=["direct", "structured"])
def test_integral_declines_unsupported_modes(structured):
    with pytest.raises(NotImplementedError, match="numeric intrinsic"):
        integrate_over_region(
            1, sp.Eq(X, 0), (X,), method="numeric", measure_dimension=0, return_result=structured
        )
    with pytest.raises(NotImplementedError, match="extra bounds"):
        integrate_over_region(
            1, IntervalRegion(0, 1), (X,), bounds={X: (0, 1)}, return_result=structured
        )
    with pytest.raises(ValueError, match="method"):
        integrate_over_region(
            1, (X >= 0) & (X <= 1), (X,), method="invalid", return_result=structured
        )


@pytest.mark.parametrize("structured", [False, True], ids=["direct", "structured"])
@pytest.mark.parametrize(
    "optimize", [semialgebraic_minimize, semialgebraic_maximize], ids=["minimum", "maximum"]
)
def test_optimization_branch_limit_failure_and_recovery(optimize, structured):
    domain = sp.Or(sp.Eq(X, -1), sp.Eq(X, 2))
    with pytest.raises(NotImplementedError, match="max_boolean_branches=1"):
        optimize(X, domain, (X,), max_boolean_branches=1, return_result=structured)
    recovered = optimize(X, domain, (X,), max_boolean_branches=2, return_result=True)
    assert recovered.value == (-1 if optimize is semialgebraic_minimize else 2)
    assert recovered.attained and recovered.certified
    assert recovered.points
    assert all(domain.subs(point) is sp.true for point in recovered.points)


@pytest.mark.parametrize(
    "invalid", [sp.Matrix([[1, 2]]), sp.Matrix([[1, 1], [0, 1]])], ids=["nonsquare", "nonsymmetric"]
)
def test_definiteness_rejects_invalid_matrices(invalid):
    with pytest.raises(ValueError, match="square|symmetric"):
        matrix_pd_on(invalid)
    with pytest.raises(ValueError, match="square|symmetric"):
        matrix_psd_on(invalid)
    assert matrix_pd_on(sp.diag(-1, 1)) is False
    assert matrix_psd_on(sp.diag(-1, 1)) is False


@pytest.mark.parametrize(
    "options, message",
    [
        ({"default_sampling_radius": -1}, "radius"),
        ({"max_random_denominator": 0}, "denominator"),
        ({"numeric_precision": 0}, "precision"),
        ({"random_attempts_min": -1}, "attempt"),
        ({"strategy": "invalid"}, "strategy"),
    ],
)
def test_sampling_rejects_invalid_controls(options, message):
    with pytest.raises(ValueError, match=message):
        sample_points(X >= 0, (X,), **options)
    assert sample_points(sp.false, (X,)) == ()


@pytest.mark.parametrize(
    "polynomial", [sp.sin(X), 1 / X, sp.sqrt(X)], ids=["transcendental", "rational", "radical"]
)
def test_root_classifier_rejects_nonpolynomials(polynomial):
    with pytest.raises(sp.PolynomialError):
        classify_real_roots(polynomial, X)
    assert classify_real_roots(X**2 + 1, X).generic_root_count == 0


def test_geometry_rejects_mismatched_coordinates():
    with pytest.raises(ValueError, match="same dimension"):
        translate(X >= 0, (1, 2), (X,))
    with pytest.raises(ValueError, match="not problem variables"):
        semialgebraic_projection(X >= 0, (Y,), (X,))
    with pytest.raises(ValueError, match="same length"):
        semialgebraic_image((X, X + 1), X >= 0, (X,), image_variables=(U,))
    with pytest.raises(ValueError, match="same length"):
        semialgebraic_preimage((X, X + 1), U >= 0, (X,), target_variables=(U,))


@pytest.mark.parametrize("names", [(U, U), ("u", "u")], ids=["symbols", "strings"])
def test_images_reject_duplicate_target_coordinates(names):
    with pytest.raises(ValueError, match="distinct"):
        semialgebraic_image((X, X + 1), X >= 0, (X,), image_variables=names)
    with pytest.raises(ValueError, match="distinct"):
        semialgebraic_preimage((X, X + 1), U >= 0, (X,), target_variables=names)


@pytest.mark.parametrize(
    "options", [{"max_cells": 1}, {"timeout": 1}], ids=["cell-budget", "time-budget"]
)
@pytest.mark.parametrize("output", ["formula", "cells", "function", "tree"])
def test_unimplemented_cad_budgets_are_never_silently_ignored(options, output):
    formula = X**2 <= 1
    # Warm caches must not bypass validation of a requested budget.
    assert cad(formula, (X,)).as_set() == sp.Interval(-1, 1)
    with pytest.raises(NotImplementedError, match="not implemented"):
        cad(formula, (X,), output=output, **options)
    with pytest.raises(NotImplementedError, match="not implemented"):
        cad(formula, (X,), output=output, strict=True, return_result=True, **options)
    result = cad(formula, (X,), output=output, return_result=True, **options)
    assert result.status == "unknown"
    assert result.formula == formula
    assert "not implemented" in result.diagnostics["reason"]
    with pytest.raises(ResourceLimitError):
        result.as_function()


@pytest.mark.parametrize("output", ["formula", "cells", "function"])
def test_limited_cad_cannot_be_used_as_a_complete_representation(output):
    formula = sp.Abs(X) <= 1
    result = cad(formula, (X,), output=output, max_preprocess_aux_vars=0, return_result=True)
    assert result.status == "unknown"
    assert result.formula == formula
    with pytest.raises(ResourceLimitError):
        result.as_function()
    with pytest.raises(ResourceLimitError):
        cad(formula, (X,), output=output, max_preprocess_aux_vars=0)
    with pytest.raises(ResourceLimitError):
        cad(
            formula, (X,), output=output, max_preprocess_aux_vars=0, strict=True, return_result=True
        )
    recovered = cad(formula, (X,), return_result=True)
    assert recovered.status == "complete"
    assert recovered.formula.as_set() == sp.Interval(-1, 1)


@pytest.mark.parametrize("error", [AssertionError, RuntimeError], ids=["assertion", "runtime"])
def test_unexpected_witness_errors_propagate(error):
    from semialg.decision import _witnesses

    def broken_sampler(*args, **kwargs):
        raise error("injected sampler defect")

    with pytest.raises(error, match="injected sampler defect"):
        _witnesses.find_validated_witness(X >= 0, (X,), sampler=broken_sampler)


def test_expected_witness_failure_returns_no_witness():
    from semialg.decision import _witnesses

    def unsupported_sampler(*args, **kwargs):
        raise NotImplementedError("injected unsupported sampling case")

    assert _witnesses.find_validated_witness(X >= 0, (X,), sampler=unsupported_sampler) is None


@pytest.mark.parametrize("point", [{X: -1}, {}, {X: 1}], ids=["outside", "incomplete", "valid"])
def test_witness_hook_validates_candidates(point):
    from semialg.decision import _witnesses

    def candidate_sampler(*args, **kwargs):
        return (point,)

    result = _witnesses.find_validated_witness(X >= 0, (X,), sampler=candidate_sampler)
    assert result == ({X: 1} if point == {X: 1} else None)
