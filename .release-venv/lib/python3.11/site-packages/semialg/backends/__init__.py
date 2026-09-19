"""Backend registry and lightweight backend accessors."""

from __future__ import annotations

from importlib import import_module

from .registry import BACKENDS, BackendSpec, get_backend

__all__ = [
    "BACKENDS",
    "BackendSpec",
    "get_backend",
    "CollinsCompleteBackend",
    "FallbackDecision",
    "collins_safe_fallback",
    "maybe_fallback",
]


def __getattr__(name: str):
    if name == "CollinsCompleteBackend":
        return import_module(f"{__name__}.collins_complete").CollinsCompleteBackend
    if name in {"FallbackDecision", "collins_safe_fallback", "maybe_fallback"}:
        return getattr(import_module(f"{__name__}.fallback"), name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
