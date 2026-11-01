"""Parameter cells and exact cylindrical fiber data for real region families."""

from dataclasses import dataclass

import sympy as sp

from .cad_algorithms.cells import extract_cylindrical_solution
from .symbolic_regions import as_semialgebraic_region
from .topology.semialgebraic import HardtFiberPiece, HardtStratum, HardtTrivialization


@dataclass(frozen=True)
class RegionFamilyStratum:
    """A parameter cell with fiber CAD cells or a complete structural cover."""

    condition: sp.Expr
    sample: dict
    fiber_dimension: int
    fiber_cells: tuple
    component_count: int | None
    betti_numbers: tuple | None = None
    fiber_cover: object = None


@dataclass(frozen=True)
class RegionFamily:
    """Exact parameter partition including empty and dimension-changing fibers.

    One-dimensional fibers additionally carry Hardt graph/band topology.
    Radial families also have structural component counts and Betti numbers.
    Other higher-dimensional fibers expose cylindrical cells and dimension;
    their cell counts do not stand for connected components or Betti numbers.
    """

    formula: sp.Expr
    variables: tuple
    parameters: tuple
    strata: tuple
    empty_condition: sp.Expr
    hardt: object = None


def region_family(region, variables, parameters, *, max_cells=None, timeout=None):
    """Partition a polynomial semialgebraic family by CAD fiber structure.

    Parameters precede fiber coordinates in CAD order. Empty fibers and all
    selected lower-dimensional cells are retained. Unknown symbols are rejected
    so a coefficient parameter cannot accidentally become a decision variable.
    """
    coordinates, params = tuple(variables), tuple(parameters)
    if (
        not coordinates
        or not params
        or any(not isinstance(v, sp.Symbol) or v.is_real is False for v in (*params, *coordinates))
        or len(set((*params, *coordinates))) != len(params) + len(coordinates)
    ):
        raise ValueError("parameters and fiber variables must be distinct nonempty real symbols")
    if any(v.is_integer is True and v.is_zero is not True for v in (*params, *coordinates)):
        raise NotImplementedError("family topology requires continuous real coordinates")
    from .region_covers import _inexact_data

    if _inexact_data(region):
        raise NotImplementedError("exact family topology requires exact coefficients")
    symbolic = as_semialgebraic_region(region, coordinates)
    if symbolic.formula.has(sp.Float):
        raise NotImplementedError("exact family topology requires exact coefficients")
    formula = symbolic.quantifier_free_formula()
    if formula.free_symbols - set((*params, *coordinates)):
        raise ValueError("declare every coefficient parameter")
    from .standard_regions import Ball, Sphere

    if isinstance(region, (Ball, Sphere)) and all(
        expr.is_polynomial(*params) for expr in (region.radius, *region.center)
    ):
        return _radial_family(region, formula, coordinates, params, max_cells, timeout)
    decomposition = _family_cad(formula, (*params, *coordinates), max_cells, timeout)
    solution = extract_cylindrical_solution(
        decomposition, (*params, *coordinates), selected_only=False
    )
    if solution.decomposition_cert is None or not solution.decomposition_cert.certified:
        raise NotImplementedError("family analysis requires a complete certified CAD")
    groups = {}
    for cell in solution.cells:
        groups.setdefault(cell.index[: len(params)], []).append(cell)
    strata = []
    for cells in groups.values():
        first = cells[0]
        selected = tuple(cell for cell in cells if cell.selected)
        condition = sp.And(*(level.as_formula() for level in first.levels[: len(params)]))
        if condition is sp.false:
            continue
        count = None
        if len(coordinates) == 1:
            positions = sorted(cell.index[-1] for cell in selected)
            count = int(bool(positions)) + sum(
                b - a > 1 for a, b in zip(positions, positions[1:], strict=False)
            )
        strata.append(
            RegionFamilyStratum(
                condition,
                {v: first.sample[v] for v in params},
                max(
                    (
                        sum(level.dimension for level in cell.levels[len(params) :])
                        for cell in selected
                    ),
                    default=-1,
                ),
                selected,
                count,
            )
        )
    empty = sp.Or(*(stratum.condition for stratum in strata if stratum.fiber_dimension < 0))
    hardt = (
        _hardt_from_strata(formula, params, coordinates[0], strata)
        if len(coordinates) == 1
        else None
    )
    return RegionFamily(formula, coordinates, params, tuple(strata), empty, hardt)


def _radial_family(region, formula, coordinates, parameters, max_cells, timeout):
    """Classify ball/sphere fibers using a parameter CAD and radial topology."""
    from .region_covers import optimization_cover
    from .standard_regions import Sphere

    parameter_cad = _family_cad(
        sp.Eq(region.radius, 0, evaluate=False), parameters, max_cells, timeout
    )
    cells = extract_cylindrical_solution(parameter_cad, parameters, selected_only=False)
    if cells.decomposition_cert is None or not cells.decomposition_cert.certified:
        raise NotImplementedError("radial family analysis requires a complete parameter CAD")
    cover = optimization_cover(region, coordinates)
    strata = []
    for cell in cells.cells:
        if cell.as_formula() is sp.false:
            continue
        radius = sp.simplify(region.radius.subs(cell.sample, simultaneous=True))
        if radius.is_negative is True:
            dimension, components, betti = -1, 0, ()
        elif radius.is_zero is True:
            dimension, components, betti = 0, 1, (1,)
        elif radius.is_positive is True:
            if isinstance(region, Sphere):
                dimension = len(coordinates) - 1
                betti = (2,) if dimension == 0 else (1, *(0 for _ in range(dimension - 1)), 1)
            else:
                dimension = len(coordinates)
                betti = (1, *(0 for _ in range(dimension)))
            components = betti[0]
        else:
            raise NotImplementedError("parameter cell did not determine radius sign")
        strata.append(
            RegionFamilyStratum(
                cell.as_formula(),
                cell.sample_point(),
                dimension,
                (),
                components,
                betti,
                cover if dimension >= 0 else None,
            )
        )
    empty = sp.Or(*(s.condition for s in strata if s.fiber_dimension < 0))
    return RegionFamily(formula, coordinates, parameters, tuple(strata), empty)


def _hardt_from_strata(formula, parameters, variable, strata):
    """Reuse an existing complete scalar-fiber CAD for graph/band topology."""
    pieces = []
    for stratum in strata:
        if not stratum.fiber_cells:
            continue
        first = stratum.fiber_cells[0]
        fiber = tuple(
            HardtFiberPiece(
                cell.levels[-1].kind,
                cell.levels[-1].lower,
                cell.levels[-1].upper,
                cell.levels[-1].sample,
                cell.index[-1],
            )
            for cell in stratum.fiber_cells
        )
        pieces.append(
            HardtStratum(
                first.index[:-1],
                parameters,
                variable,
                stratum.condition,
                stratum.sample,
                sum(level.dimension for level in first.levels[:-1]),
                fiber,
            )
        )
    image = sp.Or(*(stratum.condition for stratum in strata if stratum.fiber_dimension >= 0))
    empty = sp.Or(*(stratum.condition for stratum in strata if stratum.fiber_dimension < 0))
    return HardtTrivialization(formula, parameters, variable, tuple(pieces), image, empty)


def _family_cad(formula, variables, max_cells, timeout):
    """Build an ambient sign-invariant CAD, including all false formula cells."""
    from .cad_algorithms.decomposition import decomp_collins_complete
    from .decomposition.cylindrical import CADResult, CellSet
    from .formula import formula_polynomials, parse_formula
    from .qe.complete import evaluate_formula_on_cell

    guards = []
    for variable in variables:
        for predicate, relation in (
            ("positive", sp.Gt),
            ("negative", sp.Lt),
            ("zero", sp.Eq),
            ("nonnegative", sp.Ge),
            ("nonpositive", sp.Le),
            ("nonzero", sp.Ne),
        ):
            if getattr(variable, "is_" + predicate) is True:
                guards.append(relation(variable, 0, evaluate=False))
                break
    parsed = parse_formula(sp.And(formula, *guards, evaluate=False))
    cad = decomp_collins_complete(
        formula_polynomials(parsed) or [sp.Integer(1)],
        variables,
        max_cells=max_cells,
        timeout=timeout,
    )
    leaves = cad.cells_by_level[len(variables)]
    selected = tuple(cell for cell in leaves if evaluate_formula_on_cell(parsed, cell, variables))
    cells = CellSet(tuple(variables), selected, cad.cells_by_level, formula)
    return CADResult(formula, tuple(variables), cad, cells, "cells", None, "complete")
