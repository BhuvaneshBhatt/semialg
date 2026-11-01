from __future__ import annotations

import sympy as sp

from ..cad_algorithms.bounds import AlgebraicRootFunction
from ..internal_symbols import fresh_dummy


def primitive_fiber_polynomial(poly: sp.Expr, fiber_var: sp.Symbol) -> sp.Expr:
    """Remove coefficient content that does not affect the fibre roots."""
    poly = sp.expand(sp.sympify(poly))
    try:
        primitive = sp.Poly(poly, fiber_var).primitive()[1].as_expr()
    except (sp.PolynomialError, TypeError, ValueError):
        return poly
    return sp.expand(primitive)


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
        polynomial = primitive_fiber_polynomial(polynomial, fiber_var)
        fiber_var = sp.sympify(fiber_var)
        root_index = sp.sympify(root_index)
        if not isinstance(fiber_var, sp.Symbol) or root_index.is_Integer is not True:
            return None
        if not polynomial.free_symbols <= {fiber_var}:
            return None
        try:
            roots = sp.Poly(polynomial, fiber_var).real_roots(multiple=False)
        except (NotImplementedError, sp.PolynomialError, ValueError):
            return None
        index = int(root_index)
        if 0 <= index < len(roots):
            return roots[index][0]
        # At a collision of all real sections, their common endpoint is
        # unambiguous even though there are fewer distinct roots in that fiber.
        if len(roots) == 1 and 0 <= index < roots[0][1]:
            return roots[0][0]
        return None

    def _eval_is_real(self):
        # By definition this selector denotes an ordered real root.
        return True

    @property
    def free_symbols(self):
        return self.args[0].free_symbols - {self.args[1]}

    def _eval_subs(self, old, new):
        polynomial, fiber_var, root_index = self.args
        if old == fiber_var or fiber_var in old.free_symbols:
            return self
        if fiber_var in new.free_symbols:
            renamed = fresh_dummy(fiber_var.name, **fiber_var.assumptions0)
            polynomial = polynomial.xreplace({fiber_var: renamed})
        else:
            renamed = fiber_var
        substituted = polynomial._subs(old, new)
        if substituted == polynomial:
            return self
        return self.func(substituted, renamed, root_index)


def root_function_expr(poly: sp.Expr, fiber_var: sp.Symbol, root_index: int) -> sp.Expr:
    poly = primitive_fiber_polynomial(poly, fiber_var)
    return AlgebraicRootFunction(poly, fiber_var, int(root_index)).as_expr()


__all__ = [
    "AlgebraicRootFunction",
    "primitive_fiber_polynomial",
    "root_function_expr",
    "root_of",
]
