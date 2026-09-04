from __future__ import annotations

from collections.abc import Mapping
from types import MappingProxyType
from typing import TypeVar

K = TypeVar("K")
V = TypeVar("V")


def freeze_value(value: object) -> object:
    """Recursively freeze built-in mutable containers used in cached results."""

    if isinstance(value, Mapping):
        return MappingProxyType({key: freeze_value(item) for key, item in value.items()})
    if isinstance(value, list):
        return tuple(freeze_value(item) for item in value)
    if isinstance(value, tuple):
        return tuple(freeze_value(item) for item in value)
    if isinstance(value, set):
        return frozenset(freeze_value(item) for item in value)
    if isinstance(value, frozenset):
        return frozenset(freeze_value(item) for item in value)
    return value


def freeze_mapping(mapping: Mapping[K, V]) -> Mapping[K, V]:
    """Return a recursively immutable read-only mapping."""

    return freeze_value(mapping)  # type: ignore[return-value]


__all__ = ["freeze_mapping", "freeze_value"]
