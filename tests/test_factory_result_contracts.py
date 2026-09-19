"""Result/certificate contracts exercised through the public operation that creates them."""

from __future__ import annotations

import pickle
from dataclasses import FrozenInstanceError

import pytest
import sympy as sp

import semialg
from semialg.convexity import ConvexityCertificate
from semialg.decomposition.cylindrical import CADResult
from semialg.optimization_results import FunctionRangeResult, OptimizationResult
from semialg.parameters import (
    RootCountConditionsResult,
    SolvabilityConditionsResult,
    root_count_conditions,
)
from semialg.region_integral_results import RegionIntegralResult
from semialg.root_classification import RootClassificationResult

pytestmark = pytest.mark.contract


def test_cad_result_is_created_by_public_cad_operation() -> None:
    x = sp.Symbol("x", real=True)

    result = semialg.cad(sp.true, (x,), return_result=True)

    assert isinstance(result, CADResult)
    assert result.status == "complete"
    assert result.formula is sp.true
    assert result.variables == (x,)
    assert len(result.cell_set.cells) == 1


def test_parameter_results_are_created_by_parameter_operations() -> None:
    x, a = sp.symbols("x a", real=True)

    solvability = semialg.solvability_conditions(sp.Eq(x, a), (x,), (a,), return_result=True)
    root_counts = root_count_conditions(x**2 - a, x, (a,), return_result=True)

    assert isinstance(solvability, SolvabilityConditionsResult)
    assert solvability.formula is sp.true
    assert solvability.variables == (x,)
    assert solvability.parameters == (a,)
    assert solvability.is_unconditionally_solvable

    assert isinstance(root_counts, RootCountConditionsResult)
    assert root_counts.parameters == (a,)
    assert sp.simplify(root_counts.condition_for_count(2) ^ (a > 0)) is sp.false
    assert sp.simplify(root_counts.condition_for_count(1) ^ sp.Eq(a, 0)) is sp.false
    assert sp.simplify(root_counts.condition_for_count(0) ^ (a < 0)) is sp.false

    restored = pickle.loads(pickle.dumps(solvability))
    assert restored == solvability


def test_optimization_and_range_results_come_from_public_operations() -> None:
    x = sp.Symbol("x", real=True)
    interval = sp.And(x >= 0, x <= 1)

    optimum = semialg.semialgebraic_minimize(x, interval, (x,), return_result=True)
    image = semialg.function_range(x, interval, (x,), return_result=True)

    assert isinstance(optimum, OptimizationResult)
    assert optimum.value == 0
    assert optimum.attained is True
    assert optimum.certified is True
    assert optimum.point == {x: 0}

    assert isinstance(image, FunctionRangeResult)
    assert image.infimum == 0
    assert image.supremum == 1
    assert image.minimum_attained is True
    assert image.maximum_attained is True
    assert image.is_interval is True


def test_region_integral_result_comes_from_public_integration_operation() -> None:
    x = sp.Symbol("x", real=True)

    result = semialg.integrate_over_region(
        x,
        sp.And(x >= 0, x <= 1),
        (x,),
        return_result=True,
    )

    assert isinstance(result, RegionIntegralResult)
    assert result.value == sp.Rational(1, 2)
    assert result.exact is True
    assert result.evaluated is True
    assert result.variables == (x,)


def test_convexity_certificate_is_created_by_public_certificate_operation() -> None:
    x = sp.Symbol("x", real=True)

    certificate = semialg.convexity_certificate(sp.true, (x,))

    assert isinstance(certificate, ConvexityCertificate)
    assert certificate.outcome is True
    assert certificate.certified is True
    assert certificate.method == "trivial"
    assert certificate.variables == (x,)


def test_root_classification_result_is_created_by_public_operation() -> None:
    x = sp.Symbol("x", real=True)

    result = semialg.classify_real_roots((x - 1) ** 2 * (x + 2), x)

    assert isinstance(result, RootClassificationResult)
    assert result.cells[0].root_count == 2
    assert result.cells[0].multiplicity_pattern == (1, 2)
    assert result.method == "univariate_exact_roots"


def test_factory_results_preserve_frozen_result_contracts() -> None:
    x = sp.Symbol("x", real=True)
    result = semialg.integrate_over_region(x, sp.And(x >= 0, x <= 1), (x,), return_result=True)

    with pytest.raises(FrozenInstanceError):
        result.value = sp.Integer(7)  # type: ignore[misc]
