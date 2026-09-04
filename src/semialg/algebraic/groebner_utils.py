"""Small shared construction policy for exact Groebner bases."""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from typing import Any

import sympy as sp


def compute_groebner_basis(
    generators: Iterable[sp.Expr],
    variables: Sequence[sp.Symbol],
    *,
    order: str = "grevlex",
    domain: Any | None = None,
    extension: bool | None = None,
    modulus: int | None = None,
) -> sp.polys.polytools.GroebnerBasis:
    """Construct a Groebner basis with one consistent zero-variable guard.

    The helper deliberately does not hide SymPy failures or change coefficient
    semantics.  It centralizes option assembly so exact callers do not each
    reinvent domain/extension handling and accidentally request an invalid
    generator-free polynomial ring.
    """

    vars_tuple = tuple(variables)
    if not vars_tuple:
        raise sp.PolynomialError("Groebner basis construction requires variables")
    options: dict[str, object] = {"order": order}
    if domain is not None:
        options["domain"] = domain
    if extension is not None:
        options["extension"] = extension
    if modulus is not None:
        options["modulus"] = modulus
    return sp.groebner(tuple(generators), *vars_tuple, **options)


def lex_basis_generators(
    generators: Iterable[sp.Expr],
    variables: Sequence[sp.Symbol],
    *,
    domain: Any | None = None,
    extension: bool | None = None,
) -> tuple[sp.Expr, ...]:
    """Return a lexicographic basis using semialg's triangular variable order.

    Triangular algorithms in semialg pass the ambient variables to the
    polynomial engine in reverse order, making the last listed ambient
    variable largest.  Keeping that convention here avoids duplicate basis
    construction policy in decomposition and real-algebraic solving.
    """
    basis = compute_groebner_basis(
        generators,
        tuple(reversed(tuple(variables))),
        order="lex",
        domain=domain,
        extension=extension,
    )
    return tuple(sp.expand(poly.as_expr()) for poly in basis.polys)


__all__ = ["compute_groebner_basis", "lex_basis_generators"]
