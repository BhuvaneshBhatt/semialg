"""Statistical special-function heads missing from the supported SymPy release.

The classes intentionally implement only identities that are exact and useful to
funcprops.  Property theorems live in :mod:`funcprops.real_analysis.special_functions`.
"""

from __future__ import annotations

import sympy as sp
from sympy import S


class OwenT(sp.Function):
    """Owen's T function ``T(h, a)``."""

    nargs = 2

    @classmethod
    def eval(cls, h, a):
        h, a = map(sp.sympify, (h, a))
        if a.is_zero:
            return S.Zero
        if h.is_zero:
            return sp.atan(a) / (2 * sp.pi)
        if a.could_extract_minus_sign():
            return -cls(h, -a)
        if h.could_extract_minus_sign():
            return cls(-h, a)

    def fdiff(self, argindex=1):
        h, a = self.args
        if argindex == 1:
            return (
                -sp.exp(-(h**2) / 2)
                * sp.erf(a * h / sp.sqrt(2))
                / (2 * sp.sqrt(2 * sp.pi))
            )
        if argindex == 2:
            return sp.exp(-(h**2) * (1 + a**2) / 2) / (2 * sp.pi * (1 + a**2))
        raise ArgumentIndexError(self, argindex)


class MittagLefflerE(sp.Function):
    """Mittag-Leffler ``E_alpha(z)`` or generalized ``E_{alpha,beta}(z)``."""

    nargs = (2, 3)

    @classmethod
    def eval(cls, *args):
        if len(args) == 2:
            alpha, z = map(sp.sympify, args)
            beta = S.One
        else:
            alpha, beta, z = map(sp.sympify, args)
        if z.is_zero:
            return 1 / sp.gamma(beta)
        if alpha == 1 and beta == 1:
            return sp.exp(z)
        if alpha == 1 and beta == 2:
            return (sp.exp(z) - 1) / z

    def fdiff(self, argindex=1):
        if len(self.args) == 2:
            alpha, z = self.args
            if argindex == 2:
                # E_{a,1}'(z) = E_{a,a}(z)/a.
                return MittagLefflerE(alpha, alpha, z) / alpha
        else:
            alpha, beta, z = self.args
            if argindex == 3:
                return (
                    MittagLefflerE(alpha, beta - 1, z)
                    - (beta - 1) * MittagLefflerE(alpha, beta, z)
                ) / (alpha * z)
        raise ArgumentIndexError(self, argindex)


class RegularizedGammaP(sp.Function):
    """Regularized lower incomplete gamma ``P(a, x)=gamma(a,x)/Gamma(a)``."""

    nargs = 2

    @classmethod
    def eval(cls, a, x):
        a, x = map(sp.sympify, (a, x))
        if x is S.Infinity and a.is_positive is True:
            return S.One
        if x.is_zero and a.is_positive is True:
            return S.Zero

    def fdiff(self, argindex=1):
        a, x = self.args
        if argindex == 2:
            return x ** (a - 1) * sp.exp(-x) / sp.gamma(a)
        raise ArgumentIndexError(self, argindex)

    def _eval_rewrite_as_lowergamma(self, a, x, **kwargs):
        return sp.lowergamma(a, x) / sp.gamma(a)


class RegularizedGammaQ(sp.Function):
    """Regularized upper incomplete gamma ``Q(a, x)=Gamma(a,x)/Gamma(a)``."""

    nargs = 2

    @classmethod
    def eval(cls, a, x):
        a, x = map(sp.sympify, (a, x))
        if x is S.Infinity and a.is_positive is True:
            return S.Zero
        if x.is_zero and a.is_positive is True:
            return S.One

    def fdiff(self, argindex=1):
        a, x = self.args
        if argindex == 2:
            return -(x ** (a - 1)) * sp.exp(-x) / sp.gamma(a)
        raise ArgumentIndexError(self, argindex)

    def _eval_rewrite_as_uppergamma(self, a, x, **kwargs):
        return sp.uppergamma(a, x) / sp.gamma(a)


# Local import avoids exporting a SymPy-internal name as part of our public API.
from sympy.core.function import ArgumentIndexError  # noqa: E402

__all__ = [
    "MittagLefflerE",
    "OwenT",
    "RegularizedGammaP",
    "RegularizedGammaQ",
]
