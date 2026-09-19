from __future__ import annotations

import gc
import weakref

from semialg.cache_utils import BoundedLRU, CachePolicy
from semialg.context import computation_context


class _Artifact:
    pass


def test_operation_local_cache_does_not_cross_contexts() -> None:
    cache = BoundedLRU.operation_local(4, "test.operation")
    with computation_context():
        value = _Artifact()
        cache.put("key", value)
        assert cache.get("key") is value
    with computation_context():
        assert cache.get("key") is None


def test_heavyweight_cache_does_not_own_artifact_lifetime() -> None:
    cache = BoundedLRU.heavyweight(4, "test.heavy")
    value = _Artifact()
    reference = weakref.ref(value)
    cache.put("key", value)
    del value
    gc.collect()

    assert reference() is None
    assert cache.get("key") is None
    assert len(cache) == 0


def test_clear_generation_cannot_resurrect_previous_value() -> None:
    cache = BoundedLRU.immutable(4, "test.generation")
    old = _Artifact()
    new = _Artifact()
    with computation_context():
        cache.put("key", old)
        assert cache.get("key") is old
        cache.clear()
        assert cache.get("key") is None
        cache.put("key", new)
        assert cache.get("key") is new
    assert cache.info().policy is CachePolicy.IMMUTABLE
