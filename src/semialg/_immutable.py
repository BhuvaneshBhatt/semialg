from __future__ import annotations

from collections.abc import Mapping
from dataclasses import fields
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


def _pickle_value(value):
    if isinstance(value, MappingProxyType):
        return _FrozenMappingPayload({key: _pickle_value(item) for key, item in value.items()})
    if isinstance(value, Mapping):
        return {key: _pickle_value(item) for key, item in value.items()}
    if isinstance(value, tuple):
        return tuple(_pickle_value(item) for item in value)
    if isinstance(value, frozenset):
        return frozenset(_pickle_value(item) for item in value)
    return value


class _FrozenMappingPayload:
    def __init__(self, mapping):
        self.mapping = mapping

    def __reduce__(self):
        return freeze_mapping, (self.mapping,)


class FrozenMappingPickle:
    """Serialize dataclass constructor fields and rebuild their immutable mappings.

    Cached locators and other derived attributes are omitted. Reconstruction
    calls each dataclass's normal validation and freezing path.
    """

    def __reduce__(self):
        return type(self), tuple(
            _pickle_value(getattr(self, item.name)) for item in fields(self) if item.init
        )


__all__ = ["freeze_mapping", "freeze_value"]
