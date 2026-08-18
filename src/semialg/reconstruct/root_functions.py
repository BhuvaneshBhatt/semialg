from __future__ import annotations

from dataclasses import dataclass

import sympy as sp

root_of = sp.Function("root_of")


@dataclass(frozen=True)
class RootFunction:
    """A delineable real root of a fiber polynomial over a CAD base cell.

    ``root_index`` is zero-based and follows the sorted real-root order used by
    the lifting stack. The object is intentionally symbolic: it represents a
    variable-dependent algebraic function, not an algebraic number obtained by
    substituting a sample point into the base variables.
    """

    polynomial: sp.Expr
    fiber_var: sp.Symbol
    root_index: int
    base_vars: tuple[sp.Symbol, ...] = ()
    base_index: tuple[int, ...] | None = None

    def as_expr(self) -> sp.Expr:
        return root_function_expr(self.polynomial, self.fiber_var, self.root_index)


def root_function_expr(poly: sp.Expr, fiber_var: sp.Symbol, root_index: int) -> sp.Expr:
    """Return an exact expression for a fiber root.

    SymPy's ``RootOf`` is used for ordinary univariate polynomials. CAD also
    needs to describe roots of fiber polynomials whose coefficients depend on
    base variables. Those are algebraic functions over a base cell rather than
    algebraic numbers, so they use the package-level ``root_of`` placeholder.
    """

    expanded = sp.expand(poly)
    if expanded.free_symbols <= {fiber_var}:
        try:
            return sp.RootOf(expanded, sp.Integer(root_index - 1 if root_index > 0 else root_index))
        except Exception:
            pass
    return root_of(expanded, fiber_var, sp.Integer(root_index))


__all__ = ["RootFunction", "root_function_expr", "root_of"]
