"""Correctness-first benchmark execution with reproducible structural metrics."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field
from time import perf_counter, process_time

import sympy as sp

from .cases import BenchmarkCase, BenchmarkDomain


@dataclass(frozen=True)
class AlgorithmMetrics:
    wall_seconds: float
    cpu_seconds: float
    counters: dict[str, int | float | str] = field(default_factory=dict)


@dataclass(frozen=True)
class ConformanceResult:
    case: BenchmarkCase
    algorithm: str
    correct: bool | None
    value: object
    metrics: AlgorithmMetrics


def polynomial_data(case: BenchmarkCase) -> tuple[tuple[sp.Symbol, ...], tuple[sp.Expr, ...]]:
    symbols = sp.symbols(" ".join(case.variables), real=True)
    if isinstance(symbols, sp.Symbol):
        symbols = (symbols,)
    local = dict(zip(case.variables, symbols, strict=True))
    polys = tuple(sp.sympify(text.replace("^", "**"), locals=local) for text in case.polynomials)
    return tuple(symbols), polys


def run_groebner(case: BenchmarkCase) -> ConformanceResult:
    if case.domain is not BenchmarkDomain.GROEBNER:
        raise ValueError("run_groebner requires a Gröbner benchmark")
    source_variables, polynomials = polynomial_data(case)
    # Published CAD banks list variables highest/eliminate-first to lowest; semialg
    # stores a projection tower from lowest to highest.
    variables = tuple(reversed(source_variables))
    local = dict(zip(case.variables, source_variables, strict=True))
    polynomials = tuple(
        sp.sympify(text.replace("^", "**"), locals=local) for text in case.polynomials
    )
    cpu0, wall0 = process_time(), perf_counter()
    basis = sp.groebner(polynomials, *variables)
    metrics = AlgorithmMetrics(
        perf_counter() - wall0,
        process_time() - cpu0,
        {
            "input_polynomials": len(polynomials),
            "basis_polynomials": len(basis.polys),
            "variables": len(variables),
        },
    )
    # A Gröbner basis must reduce every input generator to zero.
    correct = all(basis.reduce(poly)[1] == 0 for poly in polynomials)
    return ConformanceResult(case, "sympy.groebner", correct, basis, metrics)


def _cad_reference(
    case: BenchmarkCase, configuration: str, observed: int
) -> dict[str, int | float | str]:
    refs = case.expected.get("reference_cell_counts", {})
    if not isinstance(refs, dict) or configuration not in refs:
        return {}
    expected = int(refs[configuration])
    return {
        "reference_cells": expected,
        "cell_delta": observed - expected,
        "cell_ratio": observed / expected if expected else float("nan"),
    }


def run_cad(case: BenchmarkCase) -> ConformanceResult:
    """Build a Collins sign-invariant CAD using the case's prescribed variable order."""
    if case.domain is not BenchmarkDomain.CAD:
        raise ValueError("run_cad requires a CAD benchmark")
    from ..cad_algorithms.decomposition import decomp_collins_complete

    variables, polynomials = polynomial_data(case)
    cpu0, wall0 = process_time(), perf_counter()
    decomposition = decomp_collins_complete(polynomials or (sp.Integer(1),), variables)
    wall, cpu = perf_counter() - wall0, process_time() - cpu0
    check = decomposition.verify_sign_invariance()
    cells = len(decomposition.cells)
    from ..algebraic.cache import algebraic_cache_stats

    algebraic_stats = algebraic_cache_stats()
    if (
        case.expected.get("source_identifier") in {"CMXYExamples:4", "CMXYExamples:6"}
        and algebraic_stats.tower_escapes != 0
    ):
        raise AssertionError(
            "fiber-certified Wilson CAD lifting escaped the algebraic tower "
            f"{algebraic_stats.tower_escapes} time(s)"
        )
    counters: dict[str, int | float | str] = {
        "cells": cells,
        "max_stack_size": decomposition.max_stack_size(),
        "variables": len(variables),
        "input_polynomials": len(polynomials),
        "cell_count_by_level": str(decomposition.cell_count_by_level()),
        "projection_poly_count_by_level": str(decomposition.proj_poly_count_by_level()),
        "root_isolation_calls": algebraic_stats.calls,
        "root_cache_hits": algebraic_stats.cache_hits,
        "comparison_hits": algebraic_stats.comparison_hits,
        "comparison_misses": algebraic_stats.comparison_misses,
        "comparison_refinements": algebraic_stats.comparison_refinements,
        "specialization_hits": algebraic_stats.specialization_hits,
        "specialization_misses": algebraic_stats.specialization_misses,
        "tower_escapes": algebraic_stats.tower_escapes,
        "family_relationship_hits": algebraic_stats.family_relationship_hits,
        "family_relationship_misses": algebraic_stats.family_relationship_misses,
        "number_field_fallbacks": algebraic_stats.number_field_fallbacks,
        "number_field_refinement_fallbacks": algebraic_stats.number_field_refinement_fallbacks,
        "refinement_requests": str(algebraic_stats.refinement_requests),
        "refinement_advances": str(algebraic_stats.refinement_advances),
        "refinement_decisive_advances": str(algebraic_stats.refinement_decisive_advances),
        "sign_refinement_advances": str(algebraic_stats.sign_refinement_advances),
        "sign_refinement_chain_lengths": str(algebraic_stats.sign_refinement_chain_lengths),
        "sign_refinement_coordinates": str(algebraic_stats.sign_refinement_coordinates),
    }
    counters.update(_cad_reference(case, "collins_sign_invariant", cells))
    return ConformanceResult(
        case,
        "collins_sign_invariant",
        check.ok,
        decomposition,
        AlgorithmMetrics(wall, cpu, counters),
    )


def _tticad_formula(case: BenchmarkCase, local: dict[str, sp.Symbol]):
    from ..formula import parse_formula

    branches = []
    for family in case.expected.get("formulae", ()):
        atoms = []
        for text in family["equational_constraints"]:
            atoms.append(sp.Eq(sp.sympify(text.replace("^", "**"), locals=local), 0))
        # The published ProjectionCAD input stores non-EC polynomials rather than
        # relation directions. A fixed >0 atom preserves the projection/sign data
        # and family structure needed by this TTICAD benchmark.
        for text in family["other_constraints"]:
            atoms.append(sp.Gt(sp.sympify(text.replace("^", "**"), locals=local), 0))
        branches.append(sp.And(*atoms) if atoms else sp.true)
    return parse_formula(sp.Or(*branches) if len(branches) > 1 else branches[0])


def run_tticad(case: BenchmarkCase) -> ConformanceResult:
    """Run the family-aware certified TTICAD path on a published TTICAD case."""
    if "tticad" not in case.tags:
        raise ValueError("run_tticad requires a TTICAD benchmark")
    from ..tticad.safe import decompose_tticad_safe

    source_variables = sp.symbols(" ".join(case.variables), real=True)
    if isinstance(source_variables, sp.Symbol):
        source_variables = (source_variables,)
    source_variables = tuple(source_variables)
    local = dict(zip(case.variables, source_variables, strict=True))
    variables = tuple(reversed(source_variables))
    formula = _tticad_formula(case, local)
    cpu0, wall0 = process_time(), perf_counter()
    result = decompose_tticad_safe(formula, variables)
    wall, cpu = perf_counter() - wall0, process_time() - cpu0
    cells = len(result.cad.cells)
    counters: dict[str, int | float | str] = {
        "cells": cells,
        "families": result.family_count,
        "used_fallback": str(result.used_fallback),
        "effective_backend": result.effective_backend,
        "cell_count_by_level": str(result.cad.cell_count_by_level()),
        "projection_poly_count_by_level": str(result.cad.proj_poly_count_by_level()),
    }
    counters.update(_cad_reference(case, "published_tticad", cells))
    correct = bool(result.complete and (result.used_fallback or result.validity.valid))
    return ConformanceResult(
        case, "tticad_certified", correct, result, AlgorithmMetrics(wall, cpu, counters)
    )


def run_with(
    case: BenchmarkCase,
    algorithm: str,
    operation: Callable[[BenchmarkCase], object],
    checker: Callable[[BenchmarkCase, object], bool] | None = None,
    counters: Callable[[object], dict[str, int | float | str]] | None = None,
) -> ConformanceResult:
    """Run a semialg algorithm while keeping correctness separate from timing."""
    cpu0, wall0 = process_time(), perf_counter()
    value = operation(case)
    wall = perf_counter() - wall0
    cpu = process_time() - cpu0
    correct = checker(case, value) if checker is not None else None
    return ConformanceResult(
        case,
        algorithm,
        correct,
        value,
        AlgorithmMetrics(wall, cpu, counters(value) if counters else {}),
    )


__all__ = [
    "AlgorithmMetrics",
    "ConformanceResult",
    "polynomial_data",
    "run_cad",
    "run_groebner",
    "run_tticad",
    "run_with",
]
