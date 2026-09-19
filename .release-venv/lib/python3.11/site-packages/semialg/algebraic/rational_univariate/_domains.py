from __future__ import annotations

import sympy as sp

from .representation import RationalUnivariateError


def require_exact_rur_domain(domain):
    """Return an exact field accepted by the RUR backend or raise."""

    if domain not in (sp.ZZ, sp.QQ) and not getattr(domain, "is_AlgebraicField", False):
        raise RationalUnivariateError(
            "RUR supports rational or exact algebraic-number coefficients, "
            f"not coefficient domain {domain}"
        )
    return domain if domain.is_Field else domain.get_field()


__all__ = ["require_exact_rur_domain"]
