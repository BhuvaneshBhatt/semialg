from __future__ import annotations

from dataclasses import dataclass, field

import sympy as sp

CACHE_VERSION = 4


@dataclass
class RootIsolationStats:
    calls: int = 0
    cache_hits: int = 0
    cache_misses: int = 0
    refinements: int = 0
    comparison_refinements: int = 0


@dataclass
class RootIsolationCache:
    """Small in-memory cache for root isolation.

    The cache is deliberately process-local so imports and test runs do not
    have hidden shutdown side effects.
    """

    roots: dict[tuple[str, str], tuple[object, ...]] = field(default_factory=dict)
    signs: dict[tuple[str, tuple[str, ...]], int] = field(default_factory=dict)
    stats: RootIsolationStats = field(default_factory=RootIsolationStats)

    def clear(self) -> None:
        self.roots.clear()
        self.signs.clear()
        self.stats = RootIsolationStats()


CACHE = RootIsolationCache()


def expr_key(expr: sp.Expr) -> str:
    return sp.srepr(sp.factor(sp.expand(expr)))


def poly_key(poly: sp.Poly) -> str:
    return f"{tuple(str(g) for g in poly.gens)}::{expr_key(poly.as_expr())}"


def root_isolation_costs() -> RootIsolationStats:
    return CACHE.stats
