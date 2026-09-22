from __future__ import annotations

from collections.abc import Mapping, Sequence
from functools import lru_cache

import sympy as sp

from ..algebraic.roots import isolate_real_roots
from ..algebraic.samples import sample_to_expr
from ..cad_algorithms.decomposition import CompleteCAD
from ..cad_algorithms.lifting.stack import CADCell
from ..cad_algorithms.polynomial_utils import exact_univariate_poly
from ..exact_arithmetic import exact_truth
from ..reconstruct.cylindrical import path_condition
from ..reconstruct.root_functions import root_of
from .sample_context import CylindricalSampleContext


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


@lru_cache(maxsize=4096)
def _closed_assignment_items(
    items: tuple[tuple[sp.Symbol, sp.Expr], ...], excluded: sp.Symbol
) -> tuple[tuple[sp.Symbol, sp.Expr], ...]:
    """Return the cached closed cylindrical prefix before ``excluded``."""
    context = CylindricalSampleContext.from_mapping(dict(items))
    return tuple(context.assignment_before(excluded).items())


@lru_cache(maxsize=4096)
def _fully_closed_assignment_items(
    items: tuple[tuple[sp.Symbol, sp.Expr], ...],
) -> tuple[tuple[sp.Symbol, sp.Expr], ...]:
    """Return the cached fully closed cylindrical sample assignment."""
    context = CylindricalSampleContext.from_mapping(dict(items))
    return tuple(context.full_assignment().items())


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
        # Close the cylindrical sample tower before substituting it into this
        # root polynomial.  Higher-level sample coordinates may themselves be
        # root functions over lower coordinates; substituting the raw mapping
        # in one pass leaves those dependencies symbolic.
        base_subs = dict(_closed_assignment_items(tuple(assignments.items()), fiber))
        # Root functions are determined by the primitive part in the fibre
        # variable.  Removing parameter content *before* specialization avoids
        # false degeneration at a boundary such as root_of(x*y, y, 0) at x=0,
        # whose continuous root branch is the root of y, not an undefined root
        # of the identically-zero specialized polynomial.
        try:
            primitive = sp.Poly(polynomial, fiber).primitive()[1].as_expr()
        except (sp.PolynomialError, TypeError, ValueError):
            primitive = polynomial
        specialized = sp.expand(primitive.subs(base_subs))
        unresolved_parameters = specialized.free_symbols - {fiber}
        if unresolved_parameters:
            names = ", ".join(sorted(symbol.name for symbol in unresolved_parameters))
            raise ValueError(
                "root_of remains parameterized after cell-sample specialization: " + names
            )
        root_index = int(index)
        try:
            # ``real_roots`` repeats roots according to multiplicity.  That is
            # essential on closure boundaries where distinct delineable root
            # branches may coalesce (e.g. +/-sqrt(1-x**2) at x=+/-1).
            roots = tuple(sp.real_roots(sp.Poly(specialized, fiber).as_expr(), fiber))
        except (sp.PolynomialError, ValueError, TypeError, NotImplementedError):
            isolated = isolate_real_roots(
                exact_univariate_poly(specialized, fiber, algebraic_extension=False)
            )
            roots = tuple(
                root.as_expr()
                for root in isolated
                for _ in range(max(1, int(getattr(root, "multiplicity", 1))))
            )
        if root_index < 0 or root_index >= len(roots):
            raise ValueError("root_of index is invalid after exact specialization")
        replacements[node] = sp.sympify(roots[root_index])
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
    needed = specialized.free_symbols.intersection(assignments)
    if needed:
        closed = dict(_fully_closed_assignment_items(tuple(assignments.items())))
        substitutions = {symbol: closed[symbol] for symbol in needed}
    else:
        substitutions = {}
    value = sp.simplify(specialized.subs(substitutions))
    return exact_truth(value)


def _truth_at_cell_sample(
    condition: sp.Expr, cell: CADCell, variables: Sequence[sp.Symbol]
) -> bool:
    assignments = cell_sample_subs(cell, variables)
    try:
        return _truth_condition_at_assignments(condition, assignments)
    except (NotImplementedError, ValueError, sp.PolynomialError):
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
