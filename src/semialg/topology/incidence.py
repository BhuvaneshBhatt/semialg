from __future__ import annotations

from collections.abc import Mapping, Sequence

import sympy as sp

from ..algebraic.roots import isolate_real_roots
from ..algebraic.samples import sample_to_expr
from ..cad.decomposition import CompleteCAD
from ..cad.lifting.stack import CADCell
from ..exact_arithmetic import exact_truth
from ..reconstruct.cylindrical import path_condition
from ..reconstruct.root_functions import root_of


def cell_dimension(cell: CADCell, cells_by_level: Mapping[int, Sequence[CADCell]]) -> int:
    """Return the Euclidean dimension of a final CAD cell.

    A sector contributes one dimension and a section contributes zero. The
    dimension is read from the complete chain of ancestors, not merely from the
    final fiber cell, because a two-dimensional section over a one-dimensional
    base is a curve, while a section over a point is a point.
    """

    dim = 0
    for level in range(1, cell.level + 1):
        ancestor = cell_at_index(cells_by_level, cell.index[:level])
        if ancestor.kind == "sector":
            dim += 1
    return dim


def cell_at_index(
    cells_by_level: Mapping[int, Sequence[CADCell]], index: tuple[int, ...]
) -> CADCell:
    level = len(index)
    for cell in cells_by_level[level]:
        if cell.index == index:
            return cell
    raise KeyError(index)


def final_cells(cad: CompleteCAD) -> tuple[CADCell, ...]:
    return cad.cells_by_level.get(len(cad.tower.variables), tuple())


def cell_sample_subs(cell: CADCell, variables: Sequence[sp.Symbol]) -> dict[sp.Symbol, sp.Expr]:
    return {var: sample_to_expr(sample) for var, sample in zip(variables, cell.sample, strict=True)}


def is_cell_in_closure(
    target: CADCell,
    source: CADCell,
    variables: Sequence[sp.Symbol],
    cells_by_level: Mapping[int, Sequence[CADCell]],
) -> bool:
    """Return whether ``target`` is contained in the closure of ``source``.

    CAD closure/incidence is recursive and can become subtle when fiber root
    functions degenerate over boundary base cells. The implementation uses the
    cylindrical closed path formula for ``source`` and tests it on the sample of
    ``target``. Because the CAD is sign-invariant for all projection/root
    boundary polynomials, this sample test is exact for the cell-level closure
    operations used by this package. A dimension guard prevents sectors from
    being treated as contained in the closure of lower-dimensional sections.
    """

    if target.level != source.level:
        return False
    if cell_dimension(target, cells_by_level) > cell_dimension(source, cells_by_level):
        return False
    if target.index == source.index:
        return True
    condition = path_condition(source, variables, cells_by_level, closed=True)
    return _truth_at_cell_sample(condition, target, variables)


def _specialize_root_functions(expr: sp.Expr, assignments: Mapping[sp.Symbol, sp.Expr]) -> sp.Expr:
    """Specialize opaque ``root_of`` nodes without substituting their fiber variable.

    A plain SymPy ``subs`` would replace the fiber symbol *inside* ``root_of``
    itself, turning ``root_of(p(x, y), y, k)`` into an uninterpretable object
    such as ``root_of(c, y0, k)``.  Specialize base variables first, isolate the
    requested root exactly, then substitute the ordinary coordinates.
    """

    replacements: dict[sp.Expr, sp.Expr] = {}
    for node in sp.preorder_traversal(expr):
        if getattr(node, "func", None) != root_of or len(node.args) != 3:
            continue
        polynomial, fiber, index = node.args
        if not isinstance(fiber, sp.Symbol) or not index.is_Integer:
            raise ValueError("malformed root_of expression in CAD topology formula")
        base_subs = {var: value for var, value in assignments.items() if var != fiber}
        specialized = sp.expand(polynomial.subs(base_subs))
        roots = isolate_real_roots(sp.Poly(specialized, fiber, domain="EX"))
        root_index = int(index)
        if root_index < 0 or root_index >= len(roots):
            raise ValueError("root_of index is invalid after exact specialization")
        replacements[node] = roots[root_index].as_expr()
    return expr.xreplace(replacements) if replacements else expr


def _truth_condition_at_assignments(
    condition: sp.Expr, assignments: Mapping[sp.Symbol, sp.Expr]
) -> bool:
    """Evaluate a CAD path condition exactly with Boolean short-circuiting."""

    if isinstance(condition, sp.And):
        return all(_truth_condition_at_assignments(arg, assignments) for arg in condition.args)
    if isinstance(condition, sp.Or):
        return any(_truth_condition_at_assignments(arg, assignments) for arg in condition.args)
    if isinstance(condition, sp.Not):
        return not _truth_condition_at_assignments(condition.args[0], assignments)
    specialized = _specialize_root_functions(condition, assignments)
    value = sp.simplify(specialized.subs(assignments))
    return exact_truth(value)


def _truth_at_cell_sample(
    condition: sp.Expr, cell: CADCell, variables: Sequence[sp.Symbol]
) -> bool:
    assignments = cell_sample_subs(cell, variables)
    try:
        return _truth_condition_at_assignments(condition, assignments)
    except (ValueError, NotImplementedError, sp.PolynomialError):
        # Incidence is certification-sensitive.  Failure to establish exact
        # truth is not evidence of incidence.
        return False


def closures_intersect(
    left: CADCell,
    right: CADCell,
    variables: Sequence[sp.Symbol],
    cells_by_level: Mapping[int, Sequence[CADCell]],
    candidates: Sequence[CADCell] | None = None,
) -> bool:
    """Return whether two CAD cell closures intersect in the final CAD level."""

    if left.index == right.index:
        return True
    ambient = candidates or tuple(cells_by_level.get(left.level, ()))
    max_dim = min(cell_dimension(left, cells_by_level), cell_dimension(right, cells_by_level))
    for cell in ambient:
        if cell_dimension(cell, cells_by_level) > max_dim:
            continue
        if is_cell_in_closure(cell, left, variables, cells_by_level) and is_cell_in_closure(
            cell, right, variables, cells_by_level
        ):
            return True
    return False


__all__ = [
    "cell_at_index",
    "cell_dimension",
    "cell_sample_subs",
    "closures_intersect",
    "final_cells",
    "is_cell_in_closure",
]
