from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import sympy as sp

from ...algebraic.signs import sign_at_sample
from .stack import CADCell


@dataclass(frozen=True)
class SignInvarianceCheck:
    ok: bool
    failures: tuple[str, ...] = ()
    checked_cells: int = 0
    checked_polynomials: int = 0


def _poly_key(poly: sp.Poly) -> str:
    return sp.sstr(sp.expand(poly.as_expr()))


def verify_recorded_signs(
    cells: Sequence[CADCell], polys: Sequence[sp.Poly]
) -> SignInvarianceCheck:
    """Recompute all recorded signs and report missing or inconsistent entries."""

    failures: list[str] = []
    expected = {_poly_key(poly): poly for poly in polys}
    checked = 0
    for cell in cells:
        missing = set(expected).difference(cell.signs)
        if missing:
            failures.append(f"cell {cell.index} is missing signs for {sorted(missing)}")
        for key, poly in expected.items():
            if key not in cell.signs:
                continue
            actual = sign_at_sample(poly, cell.sample)
            checked += 1
            if cell.signs[key] != actual:
                failures.append(
                    f"cell {cell.index} recorded sign {cell.signs[key]} for {key}, but recomputation gave {actual}"
                )
    return SignInvarianceCheck(
        ok=not failures,
        failures=tuple(failures),
        checked_cells=len(cells),
        checked_polynomials=checked,
    )


def verify_cad_sign_inv(
    cells_by_level: dict[int, tuple[CADCell, ...]], tower
) -> SignInvarianceCheck:
    """Verify sign-table completeness and consistency at every CAD level."""

    failures: list[str] = []
    checked_cells = 0
    checked_polys = 0
    for level, cells in sorted(cells_by_level.items()):
        level_polys = tower.level(level).polynomials
        result = verify_recorded_signs(cells, level_polys)
        checked_cells += result.checked_cells
        checked_polys += result.checked_polynomials
        failures.extend(f"level {level}: {failure}" for failure in result.failures)
    return SignInvarianceCheck(
        ok=not failures,
        failures=tuple(failures),
        checked_cells=checked_cells,
        checked_polynomials=checked_polys,
    )
