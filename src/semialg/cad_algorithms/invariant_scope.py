"""Track the formulas and polynomial signs proved invariant by a CAD."""

from dataclasses import replace

import sympy as sp

from ..formula import formula_polynomials, parse_formula
from .polynomial_utils import polynomial_key
from .projection.collins import squarefree_basis


def _missing_boundaries(tower, polynomials, *, original_only=False):
    """Find required factors outside the CAD's certified polynomial family."""

    variables = tower.variables
    required = squarefree_basis(
        poly if isinstance(poly, sp.Poly) else sp.Poly(poly, *variables) for poly in polynomials
    )
    if original_only:
        available_polynomials = squarefree_basis(tower.original_polynomials)
    else:
        available_polynomials = (poly for level in tower.levels for poly in level.polynomials)
    available = {polynomial_key(poly) for poly in available_polynomials}
    return tuple(poly.as_expr() for poly in required if polynomial_key(poly) not in available)


def _with_invariant(cad, invariant, *, formulas=()):
    """Attach the scope of a successful reduced-CAD certificate immutably."""

    if cad.tower.metadata.get("complete", False):
        return cad
    metadata = {**dict(cad.tower.metadata), "certified_invariant": invariant}
    if invariant == "truth":
        metadata["truth_formulas"] = tuple(dict.fromkeys(formulas))
    return replace(cad, complete=True, tower=replace(cad.tower, metadata=metadata))


def _supports_formula(cad, expression):
    """Whether sample evaluation establishes truth on every cell for this formula.

    A truth certificate applies to its recorded formulas and their complements.
    A sign certificate allows Boolean combinations of the proved polynomial
    signs, including nonzero scalar multiples and repeated factors. Full Collins
    CAD also certifies its lower-level projection factors.
    """

    if not cad.complete or expression.free_symbols - set(cad.tower.variables):
        return False
    if expression in (sp.true, sp.false):
        return True
    metadata = cad.tower.metadata
    invariant = metadata.get("certified_invariant")
    if invariant == "truth":
        return any(
            expression == formula or expression == sp.Not(formula)
            for formula in metadata["truth_formulas"]
        )
    if not metadata.get("complete", False) and invariant != "sign":
        return False
    try:
        polynomials = formula_polynomials(parse_formula(expression))
        return not _missing_boundaries(cad.tower, polynomials, original_only=invariant == "sign")
    except (sp.PolynomialError, TypeError, ValueError):
        return False
