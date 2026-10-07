"""Compatibility wrapper for exact quotient coefficient-domain validation.

The equality-only domain policy is owned by :mod:`algroots.quotient`.
"""

from algroots.errors import QuotientAlgebraError
from algroots.quotient import _require_exact_domain

from .representation import RationalUnivariateError


def require_exact_rur_domain(domain):
    try:
        return _require_exact_domain(domain)
    except QuotientAlgebraError as exc:
        raise RationalUnivariateError(str(exc)) from exc


__all__ = ["require_exact_rur_domain"]
