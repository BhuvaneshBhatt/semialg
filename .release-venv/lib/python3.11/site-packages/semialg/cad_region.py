"""Reusable higher-level operations on CAD-defined semialgebraic regions."""

from __future__ import annotations

import hashlib
from collections.abc import Callable, Sequence
from dataclasses import dataclass

import sympy as sp

from .decomposition.cylindrical import CADFunction, CADResult, CellSet, cad
from .normalization import normalize_formula
from .symbolic_regions import SemialgebraicRegion, as_semialgebraic_region


@dataclass(frozen=True)
class CADSignature:
    """Stable structural summary of a CAD-defined region."""

    ambient_dimension: int
    region_dimension: int
    selected_cell_count: int
    total_leaf_count: int
    digest: str


@dataclass(frozen=True)
class CADRegion:
    """A reusable region represented by a public :class:`CADResult`."""

    result: CADResult

    @property
    def variables(self) -> tuple[sp.Symbol, ...]:
        return tuple(self.result.variables)

    @property
    def formula(self) -> sp.Expr:
        return self.result.formula

    @property
    def cells(self):
        return self.result.cell_set.cells

    @property
    def cad(self):
        return self.result.cad

    def signature(self) -> CADSignature:
        return cad_signature(self)

    def as_region(self) -> SemialgebraicRegion:
        return SemialgebraicRegion(self.formula, self.variables, cad=self.result)

    def combine(self, *others: object, op: str | Callable[..., object] = "and") -> CADRegion:
        return cad_combine(self, *others, op=op)

    def extend(self, condition: object, *, strategy: str = "auto") -> CADRegion:
        return cad_extend(self, condition, strategy=strategy)

    def locate_point(self, point):
        """Locate a point and return its exact CAD sign vector and cell metadata."""
        from .cad_algorithms.point_location import locate_cad_point

        return locate_cad_point(self, point)

    def sign_vector(self, point) -> tuple[int, ...]:
        return self.locate_point(point).sign_vector

    def cell_complex(self):
        from .cad_algorithms.cell_complex import build_cad_cell_complex

        return build_cad_cell_complex(self)

    def euler_characteristic(self) -> int:
        return self.cell_complex().euler_characteristic()

    def integrate(self, integrand: object = 1, *, evaluate: bool = False):
        """Construct exact integrals directly from this reusable CAD decomposition."""
        from .cad_algorithms.cells import extract_cylindrical_solution
        from .cad_algorithms.integration import full_dimensional_solution_integrals

        solution = extract_cylindrical_solution(self.result, selected_only=True)
        pieces = full_dimensional_solution_integrals(solution, integrand, evaluate=evaluate)
        values = tuple(piece.integral for piece in pieces)
        return sp.Add(*values) if evaluate else pieces

    def measure(self, *, evaluate: bool = True):
        """Ambient measure from the already-built CAD, without rebuilding the region."""
        pieces = self.integrate(sp.Integer(1), evaluate=False)
        total = sp.Add(*(piece.integral for piece in pieces))
        return sp.simplify(total.doit()) if evaluate else total


def as_cad_region(
    region: object, variables: Sequence[sp.Symbol | str] | None = None, *, strategy: str = "auto"
) -> CADRegion:
    """Coerce a region/formula/CAD result to a reusable :class:`CADRegion`."""
    if isinstance(region, CADRegion):
        return region
    if isinstance(region, CADResult):
        return CADRegion(region)
    reg = as_semialgebraic_region(region, variables)
    return CADRegion(reg.ensure_cad(strategy=strategy))


def _leaf_dimension(cell, result: CADResult) -> int:
    dim = 0
    for level in range(1, len(result.variables) + 1):
        prefix = cell.index[:level]
        ancestor = next(c for c in result.cad.cells_by_level[level] if c.index == prefix)
        dim += int(ancestor.kind == "sector")
    return dim


def cad_signature(
    region: object, variables: Sequence[sp.Symbol | str] | None = None
) -> CADSignature:
    """Return a stable structural signature for a CAD-defined region."""
    reg = as_cad_region(region, variables)
    leaves = tuple(reg.result.cad.cells_by_level.get(len(reg.variables), ()))
    dims = tuple(_leaf_dimension(c, reg.result) for c in reg.cells)
    payload = (
        tuple(sp.srepr(v) for v in reg.variables),
        tuple(sorted(c.index for c in reg.cells)),
        tuple(
            sorted(
                (lvl, tuple(sp.srepr(p.as_expr()) for p in ps))
                for lvl, ps in reg.result.cad.tower.by_level.items()
            )
        ),
        tuple(
            (
                cell.index,
                cell.kind,
                cell.stack_position,
                cell.root_index,
                sp.srepr(cell.section_polynomial) if cell.section_polynomial is not None else None,
            )
            for level in sorted(reg.cad.cells_by_level)
            for cell in reg.cad.cells_by_level[level]
        ),
    )
    digest = hashlib.sha256(repr(payload).encode("utf-8")).hexdigest()[:24]
    return CADSignature(
        len(reg.variables), max(dims, default=-1), len(reg.cells), len(leaves), digest
    )


def _op_expr(op: str | Callable[..., object], formulas: Sequence[sp.Expr]) -> sp.Expr:
    if callable(op):
        return normalize_formula(op(*formulas))
    key = op.lower().replace("_", "-")
    if key in {"and", "intersection"}:
        return sp.And(*formulas)
    if key in {"or", "union"}:
        return sp.Or(*formulas)
    if key in {"xor", "symmetric-difference"}:
        return sp.Xor(*formulas)
    if key in {"difference", "minus"}:
        if len(formulas) != 2:
            raise ValueError("difference requires exactly two regions")
        return sp.And(formulas[0], sp.Not(formulas[1]))
    raise ValueError(f"unsupported CAD Boolean operation: {op!r}")


def _shared_cad_combine(
    regions: Sequence[CADRegion],
    formula: sp.Expr,
    op: str | Callable[..., object],
) -> CADRegion | None:
    """Combine selected leaves when every operand shares one CAD decomposition.

    Returning ``None`` signals that a common refinement is required. The fast
    path performs Boolean set operations on leaf indices and therefore does not
    rebuild projection or lifting data.
    """
    base = regions[0].result
    if any(
        r.variables != regions[0].variables or r.result.cad is not base.cad for r in regions[1:]
    ):
        return None
    selected_sets = [set(r.result.cell_set.indices) for r in regions]
    all_leaves = tuple(base.cad.cells_by_level.get(len(base.variables), ()))
    selected: list[object] = []
    for leaf in all_leaves:
        vals = [leaf.index in selected_set for selected_set in selected_sets]
        if callable(op):
            truth_expr = normalize_formula(op(*(sp.true if value else sp.false for value in vals)))
            truth = truth_expr == sp.true
        else:
            key = op.lower().replace("_", "-")
            if key in {"and", "intersection"}:
                truth = all(vals)
            elif key in {"or", "union"}:
                truth = any(vals)
            elif key in {"xor", "symmetric-difference"}:
                truth = sum(vals) % 2 == 1
            elif key in {"difference", "minus"}:
                truth = vals[0] and not vals[1]
            else:
                raise ValueError(f"unsupported CAD Boolean operation: {op!r}")
        if truth:
            selected.append(leaf)
    cell_set = CellSet(base.variables, tuple(selected), base.cad.cells_by_level, formula)
    result = CADResult(
        formula,
        base.variables,
        base.cad,
        cell_set,
        base.output,
        "combine",
        "complete",
        {**dict(base.diagnostics), "cad_reused": True, "combine_operands": len(regions)},
        CADFunction(base.variables, formula, base.cad, cell_set, "combine"),
    )
    return CADRegion(result)


def cad_combine(
    *regions: object, op: str | Callable[..., object] = "and", strategy: str = "auto"
) -> CADRegion:
    """Boolean-combine CAD regions, reusing a shared CAD when possible."""
    if not regions:
        raise ValueError("at least one region is required")
    regs = tuple(as_cad_region(r, strategy=strategy) for r in regions)
    if any(r.variables != regs[0].variables for r in regs[1:]):
        raise ValueError("all CAD regions must have the same variables")
    formula = _op_expr(op, tuple(r.formula for r in regs))
    shared = _shared_cad_combine(regs, formula, op)
    if shared is not None:
        return shared
    result = cad(
        formula, regs[0].variables, output="formula", strategy=strategy, return_result=True
    )
    diag = {**dict(result.diagnostics), "cad_reused": False, "combine_operands": len(regs)}
    return CADRegion(
        CADResult(
            result.formula,
            result.variables,
            result.cad,
            result.cell_set,
            result.output,
            "combine",
            result.status,
            diag,
            result.function,
        )
    )


def _condition_supported_by_cad(condition: sp.Expr, reg: CADRegion) -> bool:
    """Whether every relational boundary is already sign-invariant in the CAD."""
    from .cad_algorithms.polynomial_utils import polynomial_key

    keys = {polynomial_key(poly) for polys in reg.cad.tower.by_level.values() for poly in polys}
    atoms = condition.atoms(sp.Relational) if hasattr(sp, "Relational") else set()
    if not atoms:
        from sympy.core.relational import Relational

        atoms = condition.atoms(Relational)
    for atom in atoms:
        residual = sp.expand(atom.lhs - atom.rhs)
        try:
            poly = sp.Poly(residual, *reg.variables)
        except (sp.PolynomialError, TypeError, ValueError):
            return False
        if poly.total_degree() > 0 and polynomial_key(poly) not in keys:
            return False
    return True


def _truth_on_leaf(condition: sp.Expr, leaf, variables: Sequence[sp.Symbol]) -> bool | None:
    from .algebraic.samples import sample_to_expr

    subs = {var: sample_to_expr(value) for var, value in zip(variables, leaf.sample, strict=True)}
    try:
        value = sp.simplify(condition.subs(subs))
    except (TypeError, ValueError, NotImplementedError):
        return None
    if value in (sp.true, True):
        return True
    if value in (sp.false, False):
        return False
    return None


def cad_extend(region: object, condition: object, *, strategy: str = "auto") -> CADRegion:
    """Refine a CAD region by an additional condition.

    If the condition is already truth-invariant on all existing full-dimensional
    leaves, the existing CAD is reused. Otherwise one common refinement is built.
    """
    reg = as_cad_region(region, strategy=strategy)
    cond = normalize_formula(condition)
    all_leaves = tuple(reg.cad.cells_by_level.get(len(reg.variables), ()))
    truth = {leaf.index: _truth_on_leaf(cond, leaf, reg.variables) for leaf in all_leaves}
    formula = sp.And(reg.formula, cond)
    if _condition_supported_by_cad(cond, reg) and all(v is not None for v in truth.values()):
        old_selected = set(reg.result.cell_set.indices)
        selected = tuple(
            leaf for leaf in all_leaves if leaf.index in old_selected and truth[leaf.index]
        )
        cell_set = CellSet(reg.variables, selected, reg.cad.cells_by_level, formula)
        result = CADResult(
            formula,
            reg.variables,
            reg.cad,
            cell_set,
            reg.result.output,
            "extend",
            "complete",
            {**dict(reg.result.diagnostics), "cad_reused": True},
            CADFunction(reg.variables, formula, reg.cad, cell_set, "extend"),
        )
        return CADRegion(result)
    result = cad(formula, reg.variables, output="formula", strategy=strategy, return_result=True)
    return CADRegion(
        CADResult(
            result.formula,
            result.variables,
            result.cad,
            result.cell_set,
            result.output,
            "extend",
            result.status,
            {**dict(result.diagnostics), "cad_reused": False},
            result.function,
        )
    )


# Numerical geometry lives in focused modules; these imports define the public
# CAD-region namespace without adding wrapper calls.
from .cad_algorithms.meshing import (  # noqa: E402
    CADMesh,
    _mesh_is_conforming,  # noqa: F401
    triangulate_cad_cell,
    triangulate_cad_cells,
    triangulate_cad_region,
)
from .cad_algorithms.numerical_boundaries import (  # noqa: E402
    evaluate_delineable_curve,
    evaluate_delineable_surfaces,
)

__all__ = [
    "CADMesh",
    "CADRegion",
    "CADSignature",
    "as_cad_region",
    "cad_signature",
    "cad_combine",
    "cad_extend",
    "triangulate_cad_cell",
    "triangulate_cad_cells",
    "triangulate_cad_region",
    "evaluate_delineable_curve",
    "evaluate_delineable_surfaces",
]
