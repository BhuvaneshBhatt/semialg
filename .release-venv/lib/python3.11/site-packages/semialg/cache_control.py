"""Unified lifecycle controls for process-local semialg performance caches.

The package uses bounded caches in several independent exact
subsystems.  :func:`clear_caches` provides one supported teardown hook for
long-running processes and test suites; :func:`cache_report` exposes compact
cache occupancy/statistics without making cache contents public.
"""

from __future__ import annotations

import gc
import importlib
from collections.abc import Callable
from dataclasses import dataclass
from typing import Any

BOUNDED_CACHE_REGISTRY = frozenset(
    {
        "algebraic.boundary-root-ranks",
        "algebraic.comparisons",
        "algebraic.compiled-polys",
        "algebraic.descartes-nodes",
        "algebraic.descartes-variations",
        "algebraic.family-relationships",
        "algebraic.fiber-common-gcds",
        "algebraic.fiber-sign-certificates",
        "algebraic.fiber-sturm-data",
        "algebraic.fiber-subresultants",
        "algebraic.fiber-zero-certificates",
        "algebraic.interval-poly-terms",
        "algebraic.parent-signs",
        "algebraic.root-index-intervals",
        "algebraic.root-interval-splits",
        "algebraic.root-refinement-nodes",
        "algebraic.root-sets",
        "algebraic.root-sturm-sequences",
        "algebraic.roots",
        "algebraic.rur",
        "algebraic.sector-separators",
        "algebraic.fiber-subresultant-section-scans",
        "algebraic.ranked-boundary-certificates",
        "algebraic.ranked-partition-nodes",
        "algebraic.root-interval-certificates",
        "algebraic.root-rank-certificates",
        "algebraic.shared-fiber-partitions",
        "algebraic.sign-certificates",
        "algebraic.sign-interval-certificates",
        "algebraic.sign-tables",
        "algebraic.signs",
        "algebraic.specializations",
        "algebraic.sturm-endpoint-certificates",
        "algebraic.sturm-variations",
        "algebraic.tower-reductions",
        "cad.complete_cads",
        "cad.projection_steps",
        "cad.projection_towers",
        "cad.squarefree_bases",
        "decision.metadata",
    }
)


@dataclass(frozen=True)
class FunctoolsCacheSpec:
    """Declarative registration for one process-local ``functools`` cache."""

    namespace: str
    module: str
    attribute: str

    @property
    def source(self) -> tuple[str, str]:
        module_path = self.module.replace(".", "/") + ".py"
        return (module_path, self.attribute)


FUNCTOOLS_CACHE_REGISTRY = (
    FunctoolsCacheSpec(
        "optimization.range", "_optimization_certification", "_cached_complete_range_certification"
    ),
    FunctoolsCacheSpec(
        "integration.standard_polyhedron",
        "standard_region_integrate",
        "_polyhedron_measure_tetrahedra",
    ),
    FunctoolsCacheSpec("simplify.implication", "simplify.implication", "_is_unsatisfiable_cached"),
    FunctoolsCacheSpec(
        "planner.variable_orders", "planner.heuristics", "_candidate_variable_orders_cached"
    ),
    FunctoolsCacheSpec("cad.polynomial_keys", "cad_algorithms.polynomial_utils", "polynomial_key"),
    FunctoolsCacheSpec(
        "cad.normalized_polynomials", "cad_algorithms.polynomial_utils", "normalize_poly"
    ),
    FunctoolsCacheSpec(
        "cad.discriminants", "cad_algorithms.polynomial_utils", "_discriminant_expression_cached"
    ),
    FunctoolsCacheSpec(
        "cad.resultants", "cad_algorithms.polynomial_utils", "_resultant_expression_cached"
    ),
    FunctoolsCacheSpec("convexity.hessian", "convexity", "_hessian_data"),
    FunctoolsCacheSpec(
        "convexity.polynomial_certificate", "convexity", "_cached_polynomial_convexity_certificate"
    ),
    FunctoolsCacheSpec("convexity.definition", "convexity", "_is_convex_by_definition"),
    FunctoolsCacheSpec(
        "integer.groebner", "solve.integer.groebner_recursion", "_groebner_basis_cached"
    ),
    FunctoolsCacheSpec("gtz.qq_basis", "algebraic.gtz_primary", "_canonical_qq_basis_cached"),
    FunctoolsCacheSpec("gtz.ff_basis", "algebraic.gtz_primary", "_canonical_ff_basis_cached"),
    FunctoolsCacheSpec("sampling.ellipsoid", "region_sampling", "_ellipsoid_data"),
    FunctoolsCacheSpec("sampling.polytope", "region_sampling", "_polytope_sampling_data"),
    FunctoolsCacheSpec("sampling.weighted_pieces", "region_sampling", "_weighted_region_pieces"),
    FunctoolsCacheSpec(
        "topology.closed_assignments", "topology.incidence", "_closed_assignment_items"
    ),
    FunctoolsCacheSpec(
        "topology.fully_closed_assignments", "topology.incidence", "_fully_closed_assignment_items"
    ),
)

REGISTERED_LRU_CACHES = frozenset(spec.source for spec in FUNCTOOLS_CACHE_REGISTRY)


def _cache_info_dict(info: object) -> dict[str, int]:
    fields = ("hits", "misses", "maxsize", "currsize")
    return {name: int(getattr(info, name)) for name in fields if hasattr(info, name)}


def _functools_caches() -> tuple[tuple[str, Callable[[], None], Callable[[], object]], ...]:
    """Resolve the declarative cache registry lazily.

    This function is called only by explicit lifecycle/reporting APIs, so the
    dynamic imports do not affect ordinary mathematical execution.
    """

    resolved = []
    for spec in FUNCTOOLS_CACHE_REGISTRY:
        module = importlib.import_module(f"{__package__}.{spec.module}")
        cache = getattr(module, spec.attribute)
        resolved.append((spec.namespace, cache.cache_clear, cache.cache_info))
    return tuple(resolved)


def clear_caches(*, include_sympy: bool = False, collect: bool = True) -> None:
    """Clear all process-local semialg caches.

    Parameters
    ----------
    include_sympy:
        Also clear SymPy's global expression cache.  This is useful at explicit
        lifecycle boundaries in long-running workers, but is opt-in because it
        affects other SymPy-using code in the same process.
    collect:
        Run a full Python garbage collection after releasing cache references.

    Clearing caches changes performance only; it never changes mathematical
    semantics or certificate validity.
    """

    from .algebraic.cache import clear_algebraic_caches
    from .algebraic.gtz_primary import clear_gtz_caches
    from .cad_algorithms.performance_cache import clear_cad_caches
    from .decision._metadata import clear_solution_metadata_cache
    from .simplify.implication import clear_implication_minimization_stats

    clear_algebraic_caches()
    clear_cad_caches()
    clear_gtz_caches()
    clear_solution_metadata_cache()
    clear_implication_minimization_stats()

    # Some local clearers above overlap these LRUs.  Calling cache_clear twice
    # is harmless and keeps this registry authoritative if a local helper
    # changes independently.
    for _name, clearer, _info in _functools_caches():
        clearer()

    if include_sympy:
        from sympy.core.cache import clear_cache as clear_sympy_cache

        clear_sympy_cache()
    if collect:
        gc.collect()


def cache_report() -> dict[str, Any]:
    """Return a compact snapshot of cache occupancy and hit/miss statistics."""

    from .algebraic.cache import CACHE, algebraic_cache_stats
    from .algebraic.gtz_primary import gtz_cache_info
    from .cad_algorithms.performance_cache import (
        COMPLETE_CADS,
        PROJECTION_STEPS,
        PROJECTION_TOWERS,
        SQUAREFREE_BASES,
        cad_cache_stats,
    )

    algebraic_stats = algebraic_cache_stats()
    cad_stats = cad_cache_stats()
    from .decision import _metadata

    gtz = gtz_cache_info()
    report: dict[str, Any] = {
        "algebraic": {
            "sizes": {
                "roots": len(CACHE.roots),
                "signs": len(CACHE.signs),
                "comparisons": len(CACHE.comparisons),
                "specializations": len(CACHE.specializations),
                "rur": len(CACHE.rur),
            },
            "limits": {
                "roots": CACHE.roots.maxsize,
                "signs": CACHE.signs.maxsize,
                "comparisons": CACHE.comparisons.maxsize,
                "specializations": CACHE.specializations.maxsize,
                "rur": CACHE.rur.maxsize,
            },
            "stats": dict(algebraic_stats.__dict__),
        },
        "cad": {
            "sizes": {
                "projection_towers": len(PROJECTION_TOWERS),
                "squarefree_bases": len(SQUAREFREE_BASES),
                "projection_steps": len(PROJECTION_STEPS),
                "complete_cads": len(COMPLETE_CADS),
            },
            "limits": {
                "projection_towers": PROJECTION_TOWERS.maxsize,
                "squarefree_bases": SQUAREFREE_BASES.maxsize,
                "projection_steps": PROJECTION_STEPS.maxsize,
                "complete_cads": COMPLETE_CADS.maxsize,
            },
            "stats": dict(cad_stats.__dict__),
        },
        "metadata": {
            "size": len(_metadata._METADATA_CACHE),
            "limit": _metadata._METADATA_CACHE.maxsize,
        },
        "gtz": {name: _cache_info_dict(info) for name, info in gtz.items()},
        "functools": {},
    }
    report["functools"] = {
        name: _cache_info_dict(info()) for name, _clear, info in _functools_caches()
    }
    return report


__all__ = ["cache_report", "clear_caches"]
