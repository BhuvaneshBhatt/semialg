"""Semialgebraic reasoning, CAD, and real quantifier elimination for Python."""

from __future__ import annotations

from importlib import import_module

from ._public_api import PUBLIC_EXPORTS

__version__ = "1.2.0"

# Derive ``__all__`` from the same registry used by ``__getattr__`` so public
# star-import/documentation surfaces cannot silently drift out of sync.
__all__ = ["__version__", *PUBLIC_EXPORTS]


def __getattr__(name: str):
    module_name = PUBLIC_EXPORTS.get(name)
    if module_name is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    module = import_module(module_name, __name__)
    value = getattr(module, name)
    globals()[name] = value
    return value


def __dir__() -> list[str]:
    """Return the documented root API and standard module metadata."""
    metadata = {
        "__doc__",
        "__file__",
        "__loader__",
        "__name__",
        "__package__",
        "__path__",
        "__spec__",
    }
    return sorted(metadata | set(__all__))
