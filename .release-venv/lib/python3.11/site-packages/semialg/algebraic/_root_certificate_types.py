"""Replay-certificate models for exact real-root isolation."""

from __future__ import annotations

from dataclasses import dataclass

import sympy as sp

from .intervals import RationalInterval
from .samples import AlgebraicRoot


@dataclass(frozen=True)
class PolynomialRootIntervalCertificate:
    """Replayable exact certificate for real polynomial roots in a rational interval.

    ``root_count`` counts distinct real roots of the square-free part.  Endpoint
    roots are recorded separately so open/closed interval semantics can be
    replayed without numerical approximation.
    """

    polynomial: sp.Poly
    left: sp.Rational
    right: sp.Rational
    root_count: int
    left_is_root: bool
    right_is_root: bool
    method: str
    include_left: bool = True
    include_right: bool = True
    descartes_variations: int | None = None

    @property
    def unique(self) -> bool:
        return self.root_count == 1

    @property
    def root_free(self) -> bool:
        return self.root_count == 0


@dataclass(frozen=True)
class SignStableRootNeighborhood:
    """A rational isolating neighborhood with certified companion signs."""

    root: AlgebraicRoot
    interval: RationalInterval
    signs: tuple[int, ...]
    companion_polynomials: tuple[sp.Poly, ...]
