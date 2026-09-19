"""Per-call memoization for structural oracle facts."""

from __future__ import annotations

import sympy as sp

from ._domains import NumberKind, number_kind
from ._features import ExprFeatures, expression_features


class OracleMemo:
    """Lazy memo tables shared by recursive work in one oracle call.

    The tables are intentionally separate because ``nonzero`` and
    ``defined_nonzero`` prove different propositions.  The object is created
    only for expressions above the call-local memo threshold and is discarded
    with the oracle invocation.
    """

    __slots__ = ("defined_cache", "feature_cache", "kind_cache", "nonzero_cache")

    def __init__(self) -> None:
        self.nonzero_cache: dict[tuple[sp.Expr, object], bool | None] | None = None
        self.defined_cache: dict[tuple[sp.Expr, object], bool | None] | None = None
        self.kind_cache: dict[sp.Expr, NumberKind] | None = None
        self.feature_cache: dict[sp.Expr, ExprFeatures] | None = None

    def features(self, term: sp.Expr) -> ExprFeatures:
        """Return structural features, reusing facts within this oracle call."""
        term = sp.sympify(term)
        if self.feature_cache is None:
            self.feature_cache = {}
        cached = self.feature_cache.get(term)
        if cached is None:
            cached = expression_features(term)
            self.feature_cache[term] = cached
        return cached

    def nonzero(self, term: sp.Expr, assumptions=True) -> bool | None:
        """Return a cached cheap nonzero proof for ``term``."""
        from ._nonzero import quick_nonzero

        if self.nonzero_cache is None:
            self.nonzero_cache = {}
        return quick_nonzero(term, assumptions, context=self)

    def defined_nonzero(self, term: sp.Expr, assumptions=True) -> bool | None:
        """Return a cached proof that ``term`` is finite and nonzero."""
        from ._defined import quick_defined_nonzero

        if self.defined_cache is None:
            self.defined_cache = {}
        return quick_defined_nonzero(term, assumptions, context=self)

    def kind(self, term: sp.Expr) -> NumberKind:
        """Return the cached theorem-level number kind for ``term``."""
        term = sp.sympify(term)
        if self.kind_cache is None:
            self.kind_cache = {}
        cached = self.kind_cache.get(term)
        if cached is None:
            cached = number_kind(term)
            self.kind_cache[term] = cached
        return cached
