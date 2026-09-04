"""Section-only lifting over a zero-dimensional Groebner variety projection."""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass
from time import perf_counter

import sympy as sp

from ...algebraic.roots import isolate_real_roots
from ...algebraic.samples import AlgebraicRoot, Sample, sample_to_expr
from ...algebraic.signs import sign_at_sample
from ...errors import VarietyCADUnsupported
from ..polynomial_utils import polynomial_key
from ..projection.groebner import GroebnerVarietyProjection
from .stack import CADCell, sign_table


@dataclass(frozen=True)
class GroebnerLiftOps:
    """Injectable exact operations used by finite-variety lifting tests and backends."""

    sign_evaluator: Callable[[sp.Poly, Sequence[Sample]], int] = sign_at_sample
    fiber_builder: Callable[[sp.Expr, sp.Symbol], sp.Poly] = staticmethod(
        lambda expr, variable: sp.Poly(expr, variable, extension=True)
    )
    sign_table_builder: Callable[[Sequence[sp.Poly], Sequence[Sample]], dict[str, int]] = sign_table


DEFAULT_LIFT_OPS = GroebnerLiftOps()


def _exact_bool_at_prefix(
    expression: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    sample: tuple[Sample, ...],
) -> bool | None:
    substitutions = {
        variables[index]: sample_to_expr(sample[index]) for index in range(len(sample))
    }
    try:
        value = sp.simplify(expression.subs(substitutions))
    except (TypeError, ValueError, sp.PolynomialError):
        return None
    if value is sp.true or value == sp.true:
        return True
    if value is sp.false or value == sp.false:
        return False
    return None


def _prefix_rejected(
    residual_formula: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    sample: tuple[Sample, ...],
) -> bool:
    """Return whether already-bound residual constraints prove this branch false."""

    bound = set(variables[: len(sample)])
    factors = residual_formula.args if isinstance(residual_formula, sp.And) else (residual_formula,)
    for factor in factors:
        if not factor.free_symbols <= bound:
            continue
        truth = _exact_bool_at_prefix(factor, variables, sample)
        if truth is False:
            return True
    if residual_formula.free_symbols <= bound:
        return _exact_bool_at_prefix(residual_formula, variables, sample) is False
    return False


def _active_fiber_polynomials(
    polynomials: Sequence[sp.Poly],
    variables: tuple[sp.Symbol, ...],
    prefix: tuple[Sample, ...],
    variable: sp.Symbol,
    *,
    ops: GroebnerLiftOps,
) -> tuple[sp.Poly, ...]:
    substitutions = {
        variables[index]: sample_to_expr(prefix[index]) for index in range(len(prefix))
    }
    active: list[sp.Poly] = []
    for poly in polynomials:
        expr = sp.expand(poly.as_expr().subs(substitutions))
        if expr == 0:
            continue
        try:
            fiber = ops.fiber_builder(expr, variable)
        except (sp.PolynomialError, TypeError, ValueError) as exc:
            raise VarietyCADUnsupported(
                f"cannot construct exact fiber polynomial for {variable!s}"
            ) from exc
        if fiber.degree() > 0:
            active.append(fiber)
        elif fiber.as_expr() != 0:
            return tuple()
    return tuple(active)


def _common_real_roots(
    polynomials: Sequence[sp.Poly],
    variables: tuple[sp.Symbol, ...],
    prefix: tuple[Sample, ...],
    variable: sp.Symbol,
    *,
    ops: GroebnerLiftOps,
) -> tuple[AlgebraicRoot, ...]:
    active = _active_fiber_polynomials(polynomials, variables, prefix, variable, ops=ops)
    if not active:
        return tuple()
    # Start with the lowest-degree active equation and filter its exact roots by
    # every other equation at this lifting level.  This avoids the union-of-roots
    # behavior used by complete sign-invariant lifting.
    seed = min(active, key=lambda poly: (poly.degree(), len(poly.terms())))
    candidates = isolate_real_roots(seed)
    accepted: list[AlgebraicRoot] = []
    for root in candidates:
        sample = (*prefix, root)
        ok = True
        for original in polynomials:
            if variable not in original.as_expr().free_symbols:
                continue
            try:
                sign = ops.sign_evaluator(original, sample)
            except (ValueError, TypeError, sp.PolynomialError) as exc:
                raise VarietyCADUnsupported(
                    f"cannot certify Groebner fiber relation at {variable!s}"
                ) from exc
            if sign != 0:
                ok = False
                break
        if ok:
            accepted.append(root)
    return tuple(accepted)


def _section_cell(
    *,
    parent: CADCell | None,
    root: AlgebraicRoot,
    level: int,
    projection: GroebnerVarietyProjection,
    position: int,
    ops: GroebnerLiftOps,
) -> CADCell:
    tower = projection.tower
    prefix = parent.sample if parent is not None else tuple()
    index_prefix = parent.index if parent is not None else tuple()
    sample = (*prefix, root)
    level_polys = tower.level(level).polynomials
    defining = None
    for poly in level_polys:
        if tower.variables[level - 1] not in poly.as_expr().free_symbols:
            continue
        try:
            sign = ops.sign_evaluator(poly, sample)
        except (ValueError, TypeError, sp.PolynomialError) as exc:
            raise VarietyCADUnsupported(
                f"cannot certify section polynomial at level {level}"
            ) from exc
        if sign == 0:
            defining = poly
            break
    expression = defining.as_expr() if defining is not None else root.polynomial.as_expr()
    try:
        signs = ops.sign_table_builder(level_polys, sample)
    except (ValueError, TypeError, sp.PolynomialError) as exc:
        raise VarietyCADUnsupported(f"cannot certify sign table at level {level}") from exc
    return CADCell(
        level=level,
        index=(*index_prefix, position),
        sample=sample,
        interval=(root, root),
        section_polynomial=expression,
        signs=signs,
        kind="section",
        parent_index=parent.index if parent is not None else None,
        stack_position=position,
        root_index=root.root_index,
        defining_polynomial_key=(polynomial_key(defining) if defining is not None else None),
    )


def lift_groebner_variety(
    projection: GroebnerVarietyProjection,
    *,
    ops: GroebnerLiftOps | None = None,
) -> tuple[dict[int, tuple[CADCell, ...]], dict[str, float]]:
    """Lift only compatible algebraic sections of a finite equality variety."""

    start = perf_counter()
    ops = DEFAULT_LIFT_OPS if ops is None else ops
    tower = projection.tower
    variables = tower.variables
    cells_by_level: dict[int, tuple[CADCell, ...]] = {}
    parents: tuple[CADCell | None, ...] = (None,)
    for level, variable in enumerate(variables, start=1):
        cells: list[CADCell] = []
        for parent in parents:
            prefix = parent.sample if parent is not None else tuple()
            roots = _common_real_roots(
                tower.level(level).polynomials,
                variables,
                prefix,
                variable,
                ops=ops,
            )
            for position, root in enumerate(roots):
                cell = _section_cell(
                    parent=parent,
                    root=root,
                    level=level,
                    projection=projection,
                    position=position,
                    ops=ops,
                )
                if _prefix_rejected(projection.residual_formula, variables, cell.sample):
                    continue
                cells.append(cell)
        cells_by_level[level] = tuple(cells)
        parents = tuple(cells)
        if not parents:
            for remaining in range(level + 1, len(variables) + 1):
                cells_by_level[remaining] = tuple()
            break
    return cells_by_level, {"groebner_lifting": perf_counter() - start}


__all__ = ["GroebnerLiftOps", "lift_groebner_variety"]
