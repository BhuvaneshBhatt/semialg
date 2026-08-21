from __future__ import annotations

from dataclasses import dataclass

from ..cache_utils import BoundedLRU


@dataclass
class CADCacheStats:
    projection_tower_hits: int = 0
    projection_tower_misses: int = 0
    squarefree_hits: int = 0
    squarefree_misses: int = 0
    projection_step_hits: int = 0
    projection_step_misses: int = 0


STATS = CADCacheStats()
PROJECTION_TOWERS: BoundedLRU[object] = BoundedLRU(64, "cad.projection_towers")
SQUAREFREE_BASES: BoundedLRU[object] = BoundedLRU(256, "cad.squarefree_bases")
PROJECTION_STEPS: BoundedLRU[object] = BoundedLRU(256, "cad.projection_steps")


def clear_cad_caches() -> None:
    PROJECTION_TOWERS.clear()
    SQUAREFREE_BASES.clear()
    PROJECTION_STEPS.clear()
    fresh = CADCacheStats()
    STATS.__dict__.update(fresh.__dict__)


def cad_cache_stats() -> CADCacheStats:
    return CADCacheStats(**STATS.__dict__)


__all__ = ["CADCacheStats", "cad_cache_stats", "clear_cad_caches"]
