"""Reconstruct projected truth from an already certified ambient sign CAD."""

import sympy as sp

from ..cad_algorithms.polynomial_utils import polynomial_key
from .formula import _canonical_cleanup
from .implication import _minimize_disj_with_shared_cad, _SharedCADTruthContext


def formula_from_cad_signs(cad, variables, truth, *, max_polynomials=8, max_vectors=64):
    """Compress cell truth without rebuilding CADs for reconstructed root functions.

    The ambient CAD covers the full projection space and proves each recorded
    polynomial sign invariant. A Boolean formula in those signs exists only
    when every cell with the same full sign vector has the same requested
    truth. This condition is checked across all cells before selecting cubes.
    Exact bitset minimization then uses this same complete cell context.
    Finite-variety and truth-only reduced CADs cannot supply this argument.
    """
    level = len(variables)
    if (
        level < 2
        or not cad.tower.metadata.get("complete", False)
        or cad.tower.metadata.get("variety_only", False)
    ):
        return None
    # Inactive factors can occur at several levels. Their sign is one
    # predicate, irrespective of how many projection rows carry it.
    family = []
    seen = set()
    for k in range(1, level + 1):
        for polynomial in cad.tower.level(k).polynomials:
            key = polynomial_key(polynomial)
            if key not in seen:
                family.append((k, polynomial, key))
                seen.add(key)
    if len(family) > max_polynomials:
        return None
    cells = cad.cells_by_level[level]
    ancestors = {cell.index: cell for k in range(1, level + 1) for cell in cad.cells_by_level[k]}
    vectors = []
    outcomes = {}
    for cell in cells:
        vector = tuple(ancestors[cell.index[:k]].signs[key] for k, _poly, key in family)
        outcome = bool(truth.get(cell.index, False))
        if vector in outcomes and outcomes[vector] != outcome:
            return None
        outcomes[vector] = outcome
        vectors.append(vector)
    if sum(outcomes.values()) > max_vectors:
        return None
    context = _SharedCADTruthContext(tuple(variables), tuple(cells), {})
    atoms = {}
    for position, (_k, poly, _key) in enumerate(family):
        expression = poly.as_expr()
        for sign, atom in ((-1, expression < 0), (0, sp.Eq(expression, 0)), (1, expression > 0)):
            atoms[position, sign] = atom
            context.truth_cache[atom] = sum(
                1 << i for i, vector in enumerate(vectors) if vector[position] == sign
            )
    branches = [
        tuple(atoms[i, sign] for i, sign in enumerate(vector))
        for vector, outcome in outcomes.items()
        if outcome
    ]
    return _canonical_cleanup(_minimize_disj_with_shared_cad(branches, context))
