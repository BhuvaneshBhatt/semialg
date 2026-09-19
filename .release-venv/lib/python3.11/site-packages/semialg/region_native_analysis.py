"""Reusable native CAD analysis for repeated region invariants."""

from __future__ import annotations

from dataclasses import dataclass

import sympy as sp

from .formula import parse_formula
from .normalization import conjuncts


@dataclass
class NativeRegionAnalysis:
    formula: sp.Expr
    variables: tuple[sp.Symbol, ...]
    cad: object
    selected_cells: tuple[object, ...]

    @classmethod
    def build(cls, formula: sp.Expr, variables: tuple[sp.Symbol, ...]):
        from .decomposition.cylindrical import _build_cad_for_formula, _select_formula_cells

        parsed = parse_formula(formula)
        cad = _build_cad_for_formula(parsed, variables)
        selected = _select_formula_cells(cad, parsed, variables)
        return cls(formula, variables, cad, selected)

    def cell_dimension(self, cell) -> int:
        return sum(
            next(c for c in self.cad.cells_by_level[level] if c.index == cell.index[:level]).kind
            == "sector"
            for level in range(1, len(self.variables) + 1)
        )

    def dimension(self) -> int:
        if not self.selected_cells:
            return -1
        return max(self.cell_dimension(cell) for cell in self.selected_cells)

    def euler_characteristic(self) -> int:
        return int(sum((-1) ** self.cell_dimension(cell) for cell in self.selected_cells))

    def bounded_certificate(self) -> bool | None:
        """Return a cheap exact sufficient boundedness certificate.

        A sublevel constraint of a coercive polynomial bounds the entire
        conjunction.  We certify the common exact case where the leading
        homogeneous form is a positive multiple of a power of ||x||^2.
        """
        norm = sp.Add(*(v**2 for v in self.variables))
        for atom in conjuncts(self.formula):
            if isinstance(atom, (sp.LessThan, sp.StrictLessThan)):
                expr = sp.expand(atom.lhs - atom.rhs)
            elif isinstance(atom, (sp.GreaterThan, sp.StrictGreaterThan)):
                expr = sp.expand(atom.rhs - atom.lhs)
            else:
                continue
            try:
                poly = sp.Poly(expr, *self.variables)
            except (sp.PolynomialError, TypeError, ValueError):
                continue
            degree = poly.total_degree()
            if degree <= 0 or degree % 2:
                continue
            leading = sp.Add(
                *(
                    coeff
                    * sp.prod(v**power for v, power in zip(self.variables, monom, strict=True))
                    for monom, coeff in poly.terms()
                    if sum(monom) == degree
                )
            )
            k = degree // 2
            ratio = sp.simplify(leading / (norm**k))
            if not ratio.free_symbols and ratio.is_positive is True:
                return True
        return None
