"""Collision-free symbols for internal algebraic constructions."""

from __future__ import annotations

import sympy as sp


def fresh_dummy(prefix: str, **assumptions: object) -> sp.Dummy:
    """Return a scope-local symbol that cannot collide with user symbols."""

    return sp.Dummy(prefix, **assumptions)


def fresh_real_dummy(prefix: str) -> sp.Dummy:
    """Return a collision-free real internal symbol."""

    return fresh_dummy(prefix, real=True)


__all__ = ["fresh_dummy", "fresh_real_dummy"]
