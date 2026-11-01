"""Unified replay and diagnostic views for exact semialg results."""

from __future__ import annotations

from collections.abc import Callable, Mapping
from dataclasses import dataclass, field

import sympy as sp

from ._zero_testing import certified_equal


@dataclass(frozen=True)
class CertificateReplayResult:
    """Outcome of the checks supported by a result family's replay handler."""

    verified: bool | None
    kind: str
    method: str
    details: Mapping[str, object] = field(default_factory=dict)


@dataclass(frozen=True)
class ResultDiagnostics:
    """Uniform public diagnostic summary across solver result families."""

    kind: str
    method: str
    certified: bool | None
    counters: Mapping[str, int] = field(default_factory=dict)
    cache_stats: Mapping[str, object] = field(default_factory=dict)
    details: Mapping[str, object] = field(default_factory=dict)


ReplayHandler = Callable[[object], CertificateReplayResult]
_REPLAY_REGISTRY: dict[tuple[str, str], ReplayHandler] = {}


def _register_replay(module: str, *class_names: str):
    """Register a lazy result-family replay handler by defining module and class."""

    def decorator(handler: ReplayHandler) -> ReplayHandler:
        for class_name in class_names:
            _REPLAY_REGISTRY[(module, class_name)] = handler
        return handler

    return decorator


def _bool_equal(left: object, right: object) -> bool:
    return left is right or left == right


@_register_replay("semialg.cad_algorithms.decomposition", "CompleteCAD")
def _replay_native_cad(result):
    from .cad_algorithms import verify_cad

    if result.tower.metadata.get("variety_only") or result.tower.metadata.get("ec_lifting"):
        return CertificateReplayResult(
            None,
            "cad",
            result.backend,
            {"reason": "formula-scoped or finite-variety proof requires its statement"},
        )
    replay = verify_cad(result, return_result=True)
    return CertificateReplayResult(
        replay.valid,
        "cad",
        result.backend,
        {
            "scope": "ambient geometry",
            "geometric_replay": replay,
        },
    )


@_register_replay("semialg.decomposition.cylindrical", "CADResult")
def _replay_public_cad(result):
    """Bind native geometry, input truth, selected cells and public representations."""
    from .cache_utils import cache_disabled
    from .cad_algorithms import verify_cad
    from .decomposition.cylindrical import CADOptions, _make_result, build_cad_tree
    from .errors import ExactComputationUnavailable
    from .formula import formula_polynomials, parse_formula

    statement = result.diagnostics.get("input_formula")
    if (
        result.status != "complete"
        or statement is None
        or result.operation not in {None, "components"}
        or result.diagnostics.get("preprocessed")
        or result.cad.tower.metadata.get("variety_only")
    ):
        return CertificateReplayResult(
            None,
            "cad",
            result.cad.backend,
            {
                "reason": "adapter supports complete native polynomial cell selection",
            },
        )
    try:
        formula = parse_formula(statement)
        geometry = verify_cad(result.cad, formula=formula, polynomials=formula_polynomials(formula))
        with cache_disabled():
            expected = _make_result(
                statement,
                formula,
                result.variables,
                CADOptions(
                    operation=result.operation,
                    strategy="collins",
                    output=result.output,
                    formula_form=result.diagnostics.get("formula_form", "nested"),
                    max_formula_terms=result.diagnostics.get("max_formula_terms", 512),
                ),
            )
        verified = (
            geometry
            and result.variables == result.cad.tower.variables
            and result.cell_set == expected.cell_set
            and result.formula == expected.formula
        )
        if result.function is not None:
            function = result.function
            verified = verified and (
                function.variables == result.variables
                and function.formula == result.formula
                and function.operation == result.operation
                and function.cad == result.cad
                and function.cell_set == result.cell_set
                and function.tree == build_cad_tree(result.cad, result.cell_set)
            )
    except (ExactComputationUnavailable, sp.PolynomialError, TypeError, ValueError):
        verified = False
    return CertificateReplayResult(
        bool(verified),
        "cad",
        result.cad.backend,
        {
            "scope": "input statement, ambient geometry, cell truth and canonical output",
        },
    )


@_register_replay("semialg.partial.qe", "LazyQEResult")
def _replay_partial_qe(result):
    """Reconstruct the partial trace, then check truth using full geometric lifting."""
    from .cache_utils import cache_disabled
    from .cad_algorithms import verify_cad
    from .cad_algorithms.decomposition import decomp_from_proj_tower
    from .errors import ExactComputationUnavailable
    from .formula import formula_polynomials, parse_formula
    from .partial.qe import lazy_quantifier_eliminate
    from .qe.complete import qe_from_cad

    try:
        with cache_disabled():
            matrix = parse_formula(result.matrix)
            reconstructed = lazy_quantifier_eliminate(
                result.variables,
                result.quantifiers,
                matrix,
                return_result=True,
            )
            full = decomp_from_proj_tower(result.tower)
            geometry = verify_cad(full, polynomials=formula_polynomials(matrix))
            reference = qe_from_cad(
                full, result.variables, result.quantifiers, matrix, return_result=True
            )
        verified = (
            geometry
            and result.variables == result.tower.variables
            and result.status == "complete"
            and result.free_variables == reference.free_variables
            and result.truth_by_index == reference.cell_truth
            and result.formula == reconstructed.formula
            and result.tower == reconstructed.tower
            and result.cells_by_level == reconstructed.cells_by_level
        )
    except (ExactComputationUnavailable, sp.PolynomialError, ValueError, TypeError, KeyError):
        verified = False
    return CertificateReplayResult(
        bool(verified),
        "quantifier_elimination",
        result.backend,
        {
            "scope": "partial trace, projection geometry and full-lifting truth reference",
        },
    )


@_register_replay("semialg.model", "CADResult")
def _replay_cad(result: object) -> CertificateReplayResult:
    from .invariants import validate_cad_result

    issues = tuple(validate_cad_result(result))
    return CertificateReplayResult(
        not issues,
        "cad",
        result.projection_mode,
        {"invariant_issues": issues, "cell_levels": len(result.cells_by_level)},
    )


@_register_replay("semialg.model", "QEResult")
def _replay_qe(result: object) -> CertificateReplayResult:
    cad_replay = replay_certificate(result.cad)
    return CertificateReplayResult(
        cad_replay.verified,
        "quantifier_elimination",
        result.cad.projection_mode,
        {"cad": cad_replay, "is_sentence": result.is_sentence},
    )


@_register_replay("semialg.optimization_results", "OptimizationResult")
def _replay_optimization(result: object) -> CertificateReplayResult:
    condition = result.diagnostics.get("constraint_formula")
    if condition is None:
        return CertificateReplayResult(
            None,
            "optimization",
            result.method,
            {"reason": "constraint_formula is unavailable for this result"},
        )
    from .optimization import _is_feasible, semialgebraic_maximize, semialgebraic_minimize

    solver = semialgebraic_minimize if result.kind == "min" else semialgebraic_maximize
    replayed = solver(
        result.objective,
        condition,
        result.variables,
        return_result=True,
        certification="complete" if result.certified else "auto",
    )
    value_equal = certified_equal(replayed.value, result.value)
    verified = (
        result.kind in {"min", "max"}
        and (result.attained or not result.points)
        and value_equal is True
        and replayed.attained == result.attained
        and (not result.certified or replayed.certified)
    )
    if result.attained and result.points:
        verified = verified and all(
            set(result.variables).issubset(point)
            and _is_feasible(condition, point)
            and certified_equal(result.objective.subs(point), result.value) is True
            for point in result.points
        )
    return CertificateReplayResult(
        verified,
        "optimization",
        result.method,
        {"replayed_method": replayed.method, "replayed_certified": replayed.certified},
    )


@_register_replay("semialg.optimization_results", "FunctionRangeResult")
def _replay_function_range(result: object) -> CertificateReplayResult:
    value = result.value_symbol
    checks: list[bool] = []
    if result.infimum is not None and result.infimum not in (-sp.oo, sp.oo):
        at_lower = sp.simplify(result.formula.subs(value, result.infimum))
        if result.minimum_attained is True:
            checks.append(_bool_equal(at_lower, sp.true))
    if result.supremum is not None and result.supremum not in (-sp.oo, sp.oo):
        at_upper = sp.simplify(result.formula.subs(value, result.supremum))
        if result.maximum_attained is True:
            checks.append(_bool_equal(at_upper, sp.true))
    return CertificateReplayResult(
        all(checks) if checks else None,
        "function_range",
        result.method,
        {
            "endpoint_checks": tuple(checks),
            "scope": "endpoint consistency",
            "complete_image_proof": False,
        },
    )


@_register_replay(
    "semialg.algebraic_decomposition",
    "RegularChainDecomposition",
    "RadicalIdealResult",
    "PrimaryDecompositionResult",
    "RadicalMinimalPrimeDecomposition",
)
@_register_replay(
    "semialg._algebraic_decomposition_types",
    "RegularChainDecomposition",
    "RadicalIdealResult",
    "PrimaryDecompositionResult",
    "RadicalMinimalPrimeDecomposition",
)
def _replay_algebraic_decomposition(result: object) -> CertificateReplayResult:
    from .algebraic_decomposition import (
        PrimaryDecompositionResult,
        RadicalIdealResult,
        RadicalMinimalPrimeDecomposition,
        RegularChainDecomposition,
        verify_minimal_prime_decomposition_certificate,
        verify_primary_decomposition_certificate,
        verify_radical_ideal_certificate,
        verify_regular_chain_certificate,
    )

    if isinstance(result, RegularChainDecomposition):
        if result.certificate is None:
            return CertificateReplayResult(
                None,
                "regular_chain_decomposition",
                "recursive_regular_chain",
                {"reason": "certificate unavailable"},
            )
        return CertificateReplayResult(
            verify_regular_chain_certificate(result.certificate),
            "regular_chain_decomposition",
            "recursive_regular_chain",
            {"components": len(result.components)},
        )
    if isinstance(result, RadicalIdealResult):
        if result.certificate is None:
            return CertificateReplayResult(
                None, "radical_ideal", result.method, {"reason": "certificate unavailable"}
            )
        return CertificateReplayResult(
            verify_radical_ideal_certificate(result.certificate),
            "radical_ideal",
            result.method,
            {"generators": len(result.generators)},
        )
    if isinstance(result, PrimaryDecompositionResult):
        if result.certificate is None:
            return CertificateReplayResult(
                None,
                "primary_decomposition",
                result.method,
                {"reason": "certificate unavailable"},
            )
        return CertificateReplayResult(
            verify_primary_decomposition_certificate(result.certificate),
            "primary_decomposition",
            result.method,
            {"components": len(result.components)},
        )
    if not isinstance(result, RadicalMinimalPrimeDecomposition):
        return CertificateReplayResult(None, type(result).__name__, "unknown", {})
    if result.certificate is None:
        return CertificateReplayResult(
            None,
            "minimal_prime_decomposition",
            "certified_minimal_primes",
            {"reason": "certificate unavailable"},
        )
    return CertificateReplayResult(
        verify_minimal_prime_decomposition_certificate(result.certificate),
        "minimal_prime_decomposition",
        "certified_minimal_primes",
        {"components": len(result.components)},
    )


@_register_replay(
    "semialg.algebraic.gtz",
    "IndependentLocalizationResult",
    "LocalizationContractionResult",
    "SaturationStabilizationResult",
)
@_register_replay(
    "semialg.algebraic.gtz_zero_dim",
    "ZeroDimensionalPrimaryResult",
)
@_register_replay(
    "semialg.algebraic.gtz_primary",
    "GTZPrimaryDecompositionResult",
)
@_register_replay(
    "semialg.algebraic._gtz_primary_types",
    "GTZPrimaryDecompositionResult",
)
def _replay_gtz(result: object) -> CertificateReplayResult:
    from .algebraic.gtz import (
        IndependentLocalizationResult,
        LocalizationContractionResult,
        SaturationStabilizationResult,
        verify_contraction_certificate,
        verify_localization_certificate,
        verify_saturation_certificate,
    )
    from .algebraic.gtz_primary import (
        GTZPrimaryDecompositionResult,
        verify_gtz_primary_decomposition_certificate,
    )
    from .algebraic.gtz_zero_dim import (
        ZeroDimensionalPrimaryResult,
        verify_zero_dim_primary_certificate,
    )

    if isinstance(result, IndependentLocalizationResult):
        return CertificateReplayResult(
            verify_localization_certificate(result.certificate),
            "gtz_independent_localization",
            result.method,
            {
                "independent_variables": result.independent_variables,
                "quotient_dimension": result.quotient_dimension,
            },
        )
    if isinstance(result, LocalizationContractionResult):
        return CertificateReplayResult(
            verify_contraction_certificate(result.certificate),
            "gtz_localization_contraction",
            result.method,
            {"denominator_product": result.denominator_product},
        )
    if isinstance(result, SaturationStabilizationResult):
        return CertificateReplayResult(
            verify_saturation_certificate(result.certificate),
            "gtz_saturation_stabilization",
            result.method,
            {"exponent": result.exponent},
        )
    if isinstance(result, ZeroDimensionalPrimaryResult):
        if result.certificate is None:
            return CertificateReplayResult(
                None,
                "gtz_zero_dimensional_primary",
                result.method,
                {"reason": "certificate unavailable"},
            )
        return CertificateReplayResult(
            verify_zero_dim_primary_certificate(result.certificate),
            "gtz_zero_dimensional_primary",
            result.method,
            {"components": len(result.components)},
        )
    if not isinstance(result, GTZPrimaryDecompositionResult):
        return CertificateReplayResult(None, type(result).__name__, "unknown", {})
    if result.certificate is None:
        return CertificateReplayResult(
            None,
            "gtz_primary_decomposition",
            result.method,
            {"reason": "certificate unavailable"},
        )
    return CertificateReplayResult(
        verify_gtz_primary_decomposition_certificate(result.certificate),
        "gtz_primary_decomposition",
        result.method,
        {"components": len(result.components)},
    )


@_register_replay(
    "semialg.convexity",
    "ConvexityCertificate",
    "PolynomialConvexityCertificate",
    "QuadraticConvexityCertificate",
)
def _replay_convexity(result: object) -> CertificateReplayResult:
    from .convexity import (
        ConvexityCertificate,
        PolynomialConvexityCertificate,
        QuadraticConvexityCertificate,
        convexity_certificate,
        polynomial_convexity_certificate,
        quadratic_convexity_certificate,
    )

    if isinstance(result, ConvexityCertificate):
        fresh = convexity_certificate(result.formula, result.variables)
        return CertificateReplayResult(
            fresh.outcome == result.outcome,
            "convexity",
            result.method,
            {"replayed_method": fresh.method},
        )
    if isinstance(result, PolynomialConvexityCertificate):
        fresh = polynomial_convexity_certificate(
            result.expression,
            result.variables,
            domain=result.domain,
            sense=result.sense,
        )
        return CertificateReplayResult(
            fresh.outcome == result.outcome,
            "polynomial_convexity",
            result.method,
            {"replayed_method": fresh.method},
        )
    if not isinstance(result, QuadraticConvexityCertificate):
        return CertificateReplayResult(None, type(result).__name__, "unknown", {})
    fresh = quadratic_convexity_certificate(result.formula, result.variables)
    return CertificateReplayResult(
        fresh.outcome == result.outcome,
        "quadratic_convexity",
        result.method,
        {"replayed_method": fresh.method},
    )


def replay_certificate(result: object) -> CertificateReplayResult:
    """Replay a supported exact certificate through its registered verifier.

    Each registered handler has its own documented verification scope.
    ``verified=None`` means an unregistered type or insufficient replay
    payload; it never means that a certificate failed.
    """

    key = (type(result).__module__, type(result).__name__)
    handler = _REPLAY_REGISTRY.get(key)
    if handler is None:
        return CertificateReplayResult(None, type(result).__name__, "unknown", {})
    return handler(result)


def result_diagnostics(result: object) -> ResultDiagnostics:
    """Return one stable diagnostic schema for CAD/QE/optimization certificates."""

    from .model import CADResult, QEResult

    if isinstance(result, CADResult):
        return ResultDiagnostics(
            "cad",
            result.projection_mode,
            bool(result.lifting_validated and result.well_oriented),
            dict(result.diagnostics.counters),
            {
                name: {"hits": item.hits, "misses": item.misses}
                for name, item in result.diagnostics.cache_stats.items()
            },
            {
                "nullification_events": len(result.nullification_events),
                "fallback_projection": result.used_fallback_projection,
                "lifting_certificates": len(result.lifting_certificates),
            },
        )
    if isinstance(result, QEResult):
        cad = result_diagnostics(result.cad)
        return ResultDiagnostics(
            "quantifier_elimination",
            cad.method,
            cad.certified,
            cad.counters,
            cad.cache_stats,
            {
                **dict(cad.details),
                "cells_visited": result.cells_visited,
                "partial_evaluation": result.partial_evaluation,
                "guided_pruning": result.guided_pruning,
            },
        )
    method = str(getattr(result, "method", type(result).__name__))
    certified = getattr(result, "certified", None)
    details = dict(getattr(result, "diagnostics", {}) or {})
    return ResultDiagnostics(type(result).__name__, method, certified, {}, {}, details)


__all__ = [
    "CertificateReplayResult",
    "ResultDiagnostics",
    "replay_certificate",
    "result_diagnostics",
]
