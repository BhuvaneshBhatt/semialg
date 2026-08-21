"""Canonical polynomial helpers shared by CAD projection and lifting."""

from __future__ import annotations

import sympy as sp


def polynomial_key(poly: sp.Poly) -> str:
    """Return the canonical expression key used by CAD metadata and certificates."""

    return sp.sstr(sp.expand(poly.as_expr()))
