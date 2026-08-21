"""Shared structural validation helpers.

These helpers turn low-level Python shape/length errors into package-specific
exceptions with enough context to diagnose malformed points, coordinates, and
parallel symbolic data.
"""

from __future__ import annotations

from collections.abc import Iterable, Iterator, Sized
from typing import TypeVar

from .errors import DimensionMismatchError

_T = TypeVar("_T")


def require_same_length(
    *values: Sized,
    context: str = "parallel values",
    names: tuple[str, ...] | None = None,
) -> int:
    """Require all sized inputs to have the same length and return that length.

    Parameters
    ----------
    values:
        Sized objects whose lengths must agree.
    context:
        Human-readable operation name included in any error.
    names:
        Optional labels corresponding to ``values``.
    """

    if not values:
        return 0
    lengths = tuple(len(value) for value in values)
    if len(set(lengths)) == 1:
        return lengths[0]
    if names is not None and len(names) == len(values):
        details = ", ".join(f"{name}={length}" for name, length in zip(names, lengths, strict=True))
    else:
        details = ", ".join(str(length) for length in lengths)
    raise DimensionMismatchError(f"{context} dimension mismatch ({details})")


def require_point_dimension(
    point: Sized,
    variables: Sized,
    *,
    context: str = "point",
) -> None:
    """Require a point/coordinate sequence to match its variable sequence."""

    require_same_length(point, variables, context=context, names=("point", "variables"))


def zip_equal(
    *iterables: Iterable[_T], context: str = "parallel values"
) -> Iterator[tuple[_T, ...]]:
    """Zip iterables strictly, translating length mismatch to semialg's error type.

    Unlike checking lengths up front, this also works with generators and other
    unsized iterables.  The exception is raised when iteration reaches the first
    mismatched endpoint.
    """

    try:
        yield from zip(*iterables, strict=True)
    except ValueError as exc:
        raise DimensionMismatchError(f"{context} dimension mismatch") from exc


__all__ = ["require_point_dimension", "require_same_length", "zip_equal"]
