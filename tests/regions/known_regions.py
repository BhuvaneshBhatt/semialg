"""Curated semialgebraic regions with mathematically known invariants.

The corpus contains both cheap default cases and deliberately expensive
higher-dimensional examples.  Consumers can mark ``slow`` entries so normal
local runs retain breadth without paying for every CAD stress case.
"""

from __future__ import annotations

from dataclasses import dataclass

import sympy as sp

from semialg import SemialgebraicRegion


@dataclass(frozen=True)
class KnownRegion:
    name: str
    region: SemialgebraicRegion
    dimension: int
    components: int
    bounded: bool
    closed: bool
    euler_compact_support: int | None = None
    slow: bool = False


def corpus() -> tuple[KnownRegion, ...]:
    x, y, z = sp.symbols("x y z", real=True)
    r2 = x**2 + y**2
    rho2 = x**2 + y**2 + z**2
    return (
        KnownRegion("point", SemialgebraicRegion(sp.Eq(x, 0), (x,)), 0, 1, True, True, 1),
        KnownRegion(
            "closed_interval",
            SemialgebraicRegion(sp.And(x >= -1, x <= 1), (x,)),
            1,
            1,
            True,
            True,
            1,
        ),
        KnownRegion(
            "open_interval",
            SemialgebraicRegion(sp.And(x > -1, x < 1), (x,)),
            1,
            1,
            True,
            False,
            -1,
        ),
        KnownRegion(
            "ray",
            SemialgebraicRegion(x >= 0, (x,)),
            1,
            1,
            False,
            True,
            0,
        ),
        KnownRegion(
            "two_points",
            SemialgebraicRegion(sp.Or(sp.Eq(x, -1), sp.Eq(x, 1)), (x,)),
            0,
            2,
            True,
            True,
            2,
        ),
        KnownRegion("disk", SemialgebraicRegion(r2 <= 1, (x, y)), 2, 1, True, True, 1),
        KnownRegion("circle", SemialgebraicRegion(sp.Eq(r2, 1), (x, y)), 1, 1, True, True, 0),
        KnownRegion(
            "annulus",
            SemialgebraicRegion(sp.And(r2 >= 1, r2 <= 4), (x, y)),
            2,
            1,
            True,
            True,
            0,
        ),
        KnownRegion(
            "two_disks",
            SemialgebraicRegion(
                sp.Or((x + 2) ** 2 + y**2 <= 1, (x - 2) ** 2 + y**2 <= 1),
                (x, y),
            ),
            2,
            2,
            True,
            True,
            2,
            slow=True,
        ),
        KnownRegion(
            "filled_triangle",
            SemialgebraicRegion(sp.And(x >= 0, y >= 0, x + y <= 1), (x, y)),
            2,
            1,
            True,
            True,
            1,
        ),
        KnownRegion(
            "unit_cube",
            SemialgebraicRegion(
                sp.And(x >= 0, x <= 1, y >= 0, y <= 1, z >= 0, z <= 1),
                (x, y, z),
            ),
            3,
            1,
            True,
            True,
            1,
            slow=True,
        ),
        KnownRegion(
            "sphere",
            SemialgebraicRegion(sp.Eq(rho2, 1), (x, y, z)),
            2,
            1,
            True,
            True,
            2,
            slow=True,
        ),
        KnownRegion(
            "ball",
            SemialgebraicRegion(rho2 <= 1, (x, y, z)),
            3,
            1,
            True,
            True,
            1,
            slow=True,
        ),
        KnownRegion(
            "spherical_shell",
            SemialgebraicRegion(sp.And(rho2 >= 1, rho2 <= 4), (x, y, z)),
            3,
            1,
            True,
            True,
            2,
            slow=True,
        ),
        KnownRegion(
            "solid_torus",
            SemialgebraicRegion((rho2 + 3) ** 2 <= 16 * (x**2 + y**2), (x, y, z)),
            3,
            1,
            True,
            True,
            0,
            slow=True,
        ),
    )
