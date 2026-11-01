"""Validate exact coefficient domains for real rational univariate solving.

The underlying quotient-domain check is implemented by ``algroots``.
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
