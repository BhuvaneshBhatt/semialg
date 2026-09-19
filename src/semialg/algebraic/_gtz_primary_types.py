"""Immutable data model and replay certificates for certified GTZ decomposition."""

from __future__ import annotations

from dataclasses import dataclass

import sympy as sp

from .gtz import (
    IndependentLocalizationCertificate,
    LocalizationContractionCertificate,
    SaturationStabilizationCertificate,
)
from .gtz_zero_dim import ZeroDimensionalPrimaryCertificate


@dataclass(frozen=True)
class LocalPrimaryFactor:
    generators: tuple[sp.Expr, ...]
    radical: tuple[sp.Expr, ...]
    irreducible_factor: sp.Expr
    multiplicity: int
    separator: sp.Expr


@dataclass(frozen=True)
class FractionFieldZeroDimensionalCertificate:
    parameters: tuple[sp.Symbol, ...]
    variables: tuple[sp.Symbol, ...]
    source_generators: tuple[sp.Expr, ...]
    standard_monomials: tuple[tuple[int, ...], ...]
    trace_rank: int
    nilradical_generators: tuple[sp.Expr, ...]
    primitive_coefficients: tuple[sp.Expr, ...]
    primitive_element: sp.Expr
    characteristic_polynomial: sp.Expr
    factors: tuple[tuple[sp.Expr, int], ...]
    components: tuple[LocalPrimaryFactor, ...]


@dataclass(frozen=True)
class FractionFieldZeroDimensionalResult:
    parameters: tuple[sp.Symbol, ...]
    variables: tuple[sp.Symbol, ...]
    components: tuple[LocalPrimaryFactor, ...]
    complete: bool
    certificate: FractionFieldZeroDimensionalCertificate
    method: str = "gtz_trace_primitive_artinian"


@dataclass(frozen=True)
class GTZContractedComponentCertificate:
    equations: tuple[sp.Expr, ...]
    radical: tuple[sp.Expr, ...]
    dimension: int
    degree: int
    primary_contraction: LocalizationContractionCertificate
    radical_contraction: LocalizationContractionCertificate


@dataclass(frozen=True)
class GTZNodeCertificate:
    variables: tuple[sp.Symbol, ...]
    source_generators: tuple[sp.Expr, ...]
    dimension: int
    kind: str
    zero_dimensional_certificate: ZeroDimensionalPrimaryCertificate | None = None
    localization_certificate: IndependentLocalizationCertificate | None = None
    local_primary_certificate: FractionFieldZeroDimensionalCertificate | None = None
    contracted_components: tuple[GTZContractedComponentCertificate, ...] = tuple()
    hull_generators: tuple[sp.Expr, ...] = tuple()
    split_certificate: SaturationStabilizationCertificate | None = None
    residual: GTZNodeCertificate | None = None


@dataclass(frozen=True)
class GTZPrimaryDecompositionCertificate:
    variables: tuple[sp.Symbol, ...]
    source_generators: tuple[sp.Expr, ...]
    root: GTZNodeCertificate
    final_components: tuple[object, ...]
    irredundant: bool


@dataclass(frozen=True)
class GTZExecutionPlan:
    variables: tuple[sp.Symbol, ...]
    generator_count: int
    dimension: int
    degree: int | None
    maximum_total_degree: int
    local_degree_estimate: int | None
    use_modular_groebner: bool
    memoize_subproblems: bool
    primitive_search: str
    reason: str


@dataclass(frozen=True)
class GTZPrimaryDecompositionResult:
    variables: tuple[sp.Symbol, ...]
    components: tuple[object, ...]
    complete: bool
    irredundant: bool
    method: str
    certificate: GTZPrimaryDecompositionCertificate | None
    plan: GTZExecutionPlan | None = None
