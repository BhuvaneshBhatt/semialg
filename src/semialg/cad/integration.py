from __future__ import annotations

from dataclasses import dataclass

import sympy as sp

from .bounds import AlgebraicRootFunction, bound_expr, verify_cad_cell_bounds


@dataclass(frozen=True)
class CADCellIntegral:
    cell_index: tuple[int, ...]
    dimension: int
    ambient_dimension: int
    integrand: sp.Expr
    limits: tuple[tuple[sp.Symbol, sp.Expr, sp.Expr], ...]
    integral: sp.Expr
    intrinsic: bool = False
    metric_factor: sp.Expr = sp.Integer(1)
    certified_bounds: bool = False

    def doit(self) -> sp.Expr:
        return (
            self.integral.doit() if hasattr(self.integral, "doit") else sp.simplify(self.integral)
        )


@dataclass(frozen=True)
class IntrinsicCellStratum:
    """Regular/singular classification of one cylindrical solution cell."""

    cell: object
    cell_index: tuple[int, ...]
    dimension: int
    regular: bool
    reasons: tuple[str, ...] = ()

    @property
    def singular(self) -> bool:
        return not self.regular


@dataclass(frozen=True)
class IntrinsicStratification:
    """Explicit regular/singular stratification used by intrinsic integration."""

    strata: tuple[IntrinsicCellStratum, ...]
    target_dimension: int | None = None

    @property
    def regular_strata(self) -> tuple[IntrinsicCellStratum, ...]:
        return tuple(s for s in self.strata if s.regular)

    @property
    def singular_strata(self) -> tuple[IntrinsicCellStratum, ...]:
        return tuple(s for s in self.strata if s.singular)

    @property
    def regular_cells(self) -> tuple[object, ...]:
        return tuple(s.cell for s in self.regular_strata)

    @property
    def singular_cells(self) -> tuple[object, ...]:
        return tuple(s.cell for s in self.singular_strata)


def _typed_limit_exprs(cell: object) -> tuple[tuple[sp.Symbol, sp.Expr, sp.Expr], ...]:
    return tuple((var, bound_expr(lo), bound_expr(hi)) for var, lo, hi in cell.cylindrical_bounds())


def full_dimensional_cell_integral(
    cell: object,
    integrand: object = 1,
    *,
    evaluate: bool = False,
    require_verified: bool = True,
) -> CADCellIntegral:
    """Create an arbitrary-dimensional ambient integral for one full CAD cell."""

    if not getattr(cell, "is_full_dimensional", False):
        raise ValueError(
            "full-dimensional CAD integration requires a sector at every coordinate level"
        )
    cert = verify_cad_cell_bounds(cell)
    if require_verified and not cert.verify():
        raise ValueError("CAD cell bounds could not be verified")
    limits = tuple(reversed(_typed_limit_exprs(cell)))
    expr = sp.sympify(integrand)
    integral = sp.Integral(expr, *limits)
    if evaluate:
        integral = integral.doit()
    return CADCellIntegral(
        cell_index=tuple(cell.index),
        dimension=cell.dimension,
        ambient_dimension=len(cell.variables),
        integrand=expr,
        limits=limits,
        integral=integral,
        intrinsic=False,
        certified_bounds=cert.verify(),
    )


def full_dimensional_solution_integrals(
    solution: object,
    integrand: object = 1,
    *,
    evaluate: bool = False,
    require_verified: bool = True,
) -> tuple[CADCellIntegral, ...]:
    decomposition_cert = getattr(solution, "decomposition_cert", None)
    if require_verified and decomposition_cert is not None and not decomposition_cert.verify():
        raise ValueError("cylindrical solution decomposition could not be verified")
    return tuple(
        full_dimensional_cell_integral(
            cell, integrand, evaluate=evaluate, require_verified=require_verified
        )
        for cell in getattr(solution, "full_dimensional_cells", ())
    )


def _cell_intrinsic_regularity(cell: object) -> tuple[bool, tuple[str, ...]]:
    reasons: list[str] = []
    for level in getattr(cell, "levels", ()):
        if not getattr(level, "is_section", False):
            continue
        cert = getattr(level, "delineability", None)
        bound = level.typed_lower
        if isinstance(bound, AlgebraicRootFunction):
            if cert is None:
                cert = bound.certificate
            if cert is None or not cert.verify_regularity():
                reasons.append(
                    f"level {level.level} ({level.variable}) has no cell-wide regularity certificate"
                )
                continue
        elif cert is not None and not cert.verify_regularity():
            reasons.append(
                f"level {level.level} ({level.variable}) has no cell-wide regularity certificate"
            )
    return (not reasons, tuple(reasons))


def _verify_intrinsic_regularity(cell: object) -> bool:
    return _cell_intrinsic_regularity(cell)[0]


def stratify_intrinsic_solution(
    solution: object,
    *,
    dimension: int | None = None,
    require_verified: bool = True,
) -> IntrinsicStratification:
    """Classify CAD cells into certified regular and singular strata.

    Singular lower-dimensional cells are retained explicitly rather than being
    silently discarded.  ``dimension`` restricts the returned strata when a
    particular Hausdorff dimension is being integrated.
    """

    decomposition_cert = getattr(solution, "decomposition_cert", None)
    if require_verified and decomposition_cert is not None and not decomposition_cert.verify():
        raise ValueError("cylindrical solution decomposition could not be verified")
    cells = tuple(getattr(solution, "cells", ()))
    if dimension is not None:
        cells = tuple(cell for cell in cells if cell.dimension == dimension)
    strata: list[IntrinsicCellStratum] = []
    for cell in cells:
        bounds_cert = verify_cad_cell_bounds(cell)
        if require_verified and not bounds_cert.verify():
            raise ValueError(f"CAD cell bounds could not be verified for cell {tuple(cell.index)}")
        regular, reasons = _cell_intrinsic_regularity(cell)
        strata.append(
            IntrinsicCellStratum(cell, tuple(cell.index), cell.dimension, regular, reasons)
        )
    return IntrinsicStratification(tuple(strata), dimension)


def _section_mapping_and_jacobian(
    cell: object,
) -> tuple[
    dict[sp.Symbol, sp.Expr],
    tuple[sp.Symbol, ...],
    sp.Matrix,
    tuple[tuple[sp.Symbol, sp.Expr, sp.Expr], ...],
]:
    variables = tuple(cell.variables)
    free_vars = tuple(level.variable for level in cell.levels if level.is_sector)
    free_pos = {var: i for i, var in enumerate(free_vars)}
    mapping: dict[sp.Symbol, sp.Expr] = {}
    derivative_rows: dict[sp.Symbol, list[sp.Expr]] = {}
    limits: list[tuple[sp.Symbol, sp.Expr, sp.Expr]] = []
    for level in cell.levels:
        var = level.variable
        if level.is_sector:
            mapping[var] = var
            row = [sp.Integer(0)] * len(free_vars)
            row[free_pos[var]] = sp.Integer(1)
            derivative_rows[var] = row
            lo = bound_expr(level.typed_lower).subs(mapping)
            hi = bound_expr(level.typed_upper).subs(mapping)
            limits.append((var, sp.simplify(lo), sp.simplify(hi)))
            continue
        bound = level.typed_lower
        root_expr = bound_expr(bound)
        root_expr = sp.simplify(root_expr.subs(mapping))
        mapping[var] = root_expr
        row: list[sp.Expr] = []
        if isinstance(bound, AlgebraicRootFunction):
            p = sp.sympify(bound.polynomial)
            subs_all = {k: mapping[k] for k in mapping}
            denominator = sp.simplify(sp.diff(p, var).subs(subs_all))
            if denominator == 0:
                raise ValueError(
                    "singular algebraic section cannot be integrated as a regular graph"
                )
            for q in free_vars:
                numerator = sp.Integer(0)
                for base in variables[: level.level - 1]:
                    if base in derivative_rows:
                        numerator += (
                            sp.diff(p, base).subs(subs_all) * derivative_rows[base][free_pos[q]]
                        )
                row.append(sp.simplify(-numerator / denominator))
        else:
            row = [sp.simplify(sp.diff(root_expr, q)) for q in free_vars]
        derivative_rows[var] = row
    jacobian = (
        sp.Matrix([derivative_rows[var] for var in variables])
        if free_vars
        else sp.zeros(len(variables), 0)
    )
    return mapping, free_vars, jacobian, tuple(reversed(limits))


def intrinsic_cell_integral(
    cell: object,
    integrand: object = 1,
    *,
    evaluate: bool = False,
    require_verified: bool = True,
) -> CADCellIntegral:
    """Integrate with intrinsic Hausdorff measure on a cylindrical CAD cell.

    Sections are treated as triangular graph coordinates over the sector
    coordinates.  The induced metric factor is ``sqrt(det(J.T*J))``.
    """

    cert = verify_cad_cell_bounds(cell)
    if require_verified and not cert.verify():
        raise ValueError("CAD cell bounds could not be verified")
    if require_verified and not _verify_intrinsic_regularity(cell):
        raise ValueError(
            "intrinsic graph integration requires certified regular algebraic sections"
        )
    mapping, free_vars, jacobian, limits = _section_mapping_and_jacobian(cell)
    expr = sp.sympify(integrand).subs(mapping)
    if not free_vars:
        integral = sp.simplify(expr)
        metric = sp.Integer(1)
    else:
        gram = sp.simplify(jacobian.T * jacobian)
        metric = sp.simplify(sp.sqrt(gram.det()))
        integral = sp.Integral(sp.simplify(expr * metric), *limits)
        if evaluate:
            integral = integral.doit()
    return CADCellIntegral(
        cell_index=tuple(cell.index),
        dimension=cell.dimension,
        ambient_dimension=len(cell.variables),
        integrand=sp.sympify(integrand),
        limits=limits,
        integral=integral,
        intrinsic=True,
        metric_factor=metric,
        certified_bounds=cert.verify(),
    )


def intrinsic_solution_integrals(
    solution: object,
    integrand: object = 1,
    *,
    dimension: int | None = None,
    evaluate: bool = False,
    require_verified: bool = True,
) -> tuple[CADCellIntegral, ...]:
    cells = tuple(getattr(solution, "cells", ()))
    if dimension is None and cells:
        dimension = max(cell.dimension for cell in cells)
    stratification = stratify_intrinsic_solution(
        solution, dimension=dimension, require_verified=require_verified
    )
    same_dim_singular = tuple(
        s for s in stratification.singular_strata if dimension is None or s.dimension == dimension
    )
    if require_verified and same_dim_singular:
        indices = ", ".join(str(s.cell_index) for s in same_dim_singular)
        raise ValueError(
            f"intrinsic integration encountered singular strata: {indices}; use stratify_intrinsic_solution() to inspect them"
        )
    return tuple(
        intrinsic_cell_integral(
            s.cell, integrand, evaluate=evaluate, require_verified=require_verified
        )
        for s in stratification.regular_strata
    )


__all__ = [
    "CADCellIntegral",
    "IntrinsicCellStratum",
    "IntrinsicStratification",
    "stratify_intrinsic_solution",
    "full_dimensional_cell_integral",
    "full_dimensional_solution_integrals",
    "intrinsic_cell_integral",
    "intrinsic_solution_integrals",
]
