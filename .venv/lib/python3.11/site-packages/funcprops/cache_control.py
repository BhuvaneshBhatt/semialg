"""Unified cache management for long-lived and benchmarked sessions."""

from __future__ import annotations


def _cache_registry():
    from .assumptions import _entails_cached
    from .branchcuts.computation import _mpmath_numeric_callable
    from .conditions import parse_catalogue_expression
    from .dlmf import _manifest
    from .proof_cache import (
        _derivative_cached,
        _domain_cached,
        _endpoint_limit_cached,
        _factored_derivative_cached,
        _proof_context_cached,
        _sign_cached,
        _singularity_cached,
        _zero_geometry_cached,
    )
    from .real_analysis.monotonic_certificates import _monotonic_partition_on_set
    from .registry import registry
    from .singularities import (
        _cached_complexify_equality_locus,
        _canonical_multivar_records,
    )
    from .singularities.records import _denominator_singularity_records
    from .verification import _certificate_verifiers

    return {
        "assumption_entailment": _entails_cached,
        "branch_numeric_callable": _mpmath_numeric_callable,
        "catalogue_parse": parse_catalogue_expression,
        "dlmf_manifest": _manifest,
        "derivative": _derivative_cached,
        "domain": _domain_cached,
        "endpoint_limit": _endpoint_limit_cached,
        "factored_derivative": _factored_derivative_cached,
        "proof_context": _proof_context_cached,
        "sign": _sign_cached,
        "singularity": _singularity_cached,
        "zero_geometry": _zero_geometry_cached,
        "monotonic_partition": _monotonic_partition_on_set,
        "registry": registry,
        "multivariate_registry": _canonical_multivar_records,
        "complexified_locus": _cached_complexify_equality_locus,
        "denominator_singularities": _denominator_singularity_records,
        "certificate_verifiers": _certificate_verifiers,
    }


def clear_caches():
    """Clear all bounded package caches registered for explicit lifecycle control."""
    for cached in _cache_registry().values():
        cached.cache_clear()


def cache_info():
    """Return ``functools`` cache statistics for every registered package cache."""
    return {name: cached.cache_info() for name, cached in _cache_registry().items()}
