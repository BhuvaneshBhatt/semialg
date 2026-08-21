from __future__ import annotations

import sympy as sp

from ..cad.bounds import AlgebraicRootFunction


class root_of(sp.Function):
    """Opaque ordered real-root expression used in reconstructed formulas.

    The fiber variable is binder-like: substitutions for that symbol must not
    rewrite the polynomial inside the root selector.  Once all base parameters
    are specialized and the polynomial is genuinely univariate in the fiber,
    the node evaluates to the requested exact real root.
    """

    nargs = 3

    @classmethod
    def eval(cls, polynomial, fiber_var, root_index):
        polynomial = sp.sympify(polynomial)
        fiber_var = sp.sympify(fiber_var)
        root_index = sp.sympify(root_index)
        if not isinstance(fiber_var, sp.Symbol) or root_index.is_Integer is not True:
            return None
        if not polynomial.free_symbols <= {fiber_var}:
            return None
        try:
            roots = sp.real_roots(polynomial, fiber_var)
        except (NotImplementedError, sp.PolynomialError, ValueError):
            return None
        index = int(root_index)
        if 0 <= index < len(roots):
            return roots[index]
        return None

    def _eval_is_real(self):
        # By definition this selector denotes an ordered real root.
        return True

    def _eval_subs(self, old, new):
        # Treat the fiber variable as bound within this root selector.
        if old == self.args[1]:
            return self
        return super()._eval_subs(old, new)


# ``RootFunction`` denotes the typed CAD root-function object; ``root_of`` is
# only its exact symbolic presentation when no simpler expression is certified.
RootFunction = AlgebraicRootFunction


def root_function_expr(poly: sp.Expr, fiber_var: sp.Symbol, root_index: int) -> sp.Expr:
    return AlgebraicRootFunction(sp.expand(poly), fiber_var, int(root_index)).as_expr()


__all__ = ["AlgebraicRootFunction", "RootFunction", "root_function_expr", "root_of"]
