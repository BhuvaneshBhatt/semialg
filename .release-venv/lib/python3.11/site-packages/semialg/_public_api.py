"""Registry for the primary top-level :mod:`semialg` API."""

from __future__ import annotations

from ._api_policy import PRIMARY_EXPORTS

PUBLIC_EXPORTS = PRIMARY_EXPORTS
PUBLIC_NAMES = frozenset(PUBLIC_EXPORTS)
