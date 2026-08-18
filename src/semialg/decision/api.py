from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence

import sympy as sp
from sympy.logic.boolalg import Boolean
from sympy.polys.polyerrors import PolynomialError

from ..algebraic.rational_univariate import solve_formula_with_rur
from ..decision_diagnostics import solution_capability_diagnostics
from ..formula import ParsedPrenexFormula, parse_formula
from ..instances.real_fallbacks import satisfies_formula
from ..qe import qe_by_complete_cad
from ..sampling import sample_points
from ..solve import reduce_formula
from .sampling_helpers import (
    _collect_structural_samples,
    _normalize_sample_request,
)
from .solution import (
    EquivalenceResult,
    ImplicationResult,
    IntervalComponent,
    SatisfiabilityResult,
    SemialgebraicSolution,
    TautologyResult,
)

FormulaLike = sp.Expr | Boolean | bool


def _safe_simplify_expr(expr: sp.Expr) -> sp.Expr:
    """Simplify algebraic expressions without sending Boolean formulas to radsimp."""

    if isinstance(expr, Boolean):
        try:
            return sp.simplify_logic(expr, form="dnf")
        except Exception:
            return expr
    return sp.simplify(expr)


def _merge_variables(
    first: Sequence[sp.Symbol],
    second: Sequence[sp.Symbol],
    *formulas: sp.Expr,
) -> tuple[sp.Symbol, ...]:
    """Merge explicit variables and free symbols while preserving order."""

    out: list[sp.Symbol] = []
    seen: set[sp.Symbol] = set()
    for source in (first, second):
        for sym in source:
            if sym not in seen:
                out.append(sym)
                seen.add(sym)
    for formula in formulas:
        for sym in sorted(getattr(formula, "free_symbols", set()), key=lambda item: item.name):
            if sym not in seen:
                out.append(sym)
                seen.add(sym)
    return tuple(out)


def _as_real_symbol(var: sp.Symbol | str) -> sp.Symbol:
    return sp.Symbol(var, real=True) if isinstance(var, str) else var


def _normalize_formula(formula: FormulaLike | Iterable[FormulaLike]) -> sp.Expr:
    if isinstance(formula, (list, tuple, set, frozenset)):
        pieces = [sp.sympify(piece) for piece in formula]
        return sp.And(*pieces) if pieces else sp.true
    if formula is True:
        return sp.true
    if formula is False:
        return sp.false
    if not isinstance(formula, (sp.Basic, Boolean)):
        return sp.sympify(formula)
    return formula


def _normalize_variables(
    variables: Sequence[sp.Symbol | str] | None,
    formula: sp.Expr,
) -> tuple[sp.Symbol, ...]:
    out: list[sp.Symbol] = []
    seen: set[sp.Symbol] = set()
    if variables is not None:
        for var in variables:
            sym = _as_real_symbol(var)
            if sym not in seen:
                out.append(sym)
                seen.add(sym)
    for sym in sorted(formula.free_symbols, key=lambda item: item.name):
        if sym not in seen:
            out.append(sym)
            seen.add(sym)
    return tuple(out)


def _normalize_symbols(symbols: Sequence[sp.Symbol | str] | None) -> tuple[sp.Symbol, ...]:
    out: list[sp.Symbol] = []
    seen: set[sp.Symbol] = set()
    for item in symbols or ():
        sym = _as_real_symbol(item)
        if sym not in seen:
            out.append(sym)
            seen.add(sym)
    return tuple(out)


def _normalize_solve_variables(
    variables: Sequence[sp.Symbol | str] | None,
    formula: sp.Expr,
    parameters: Sequence[sp.Symbol] = (),
) -> tuple[sp.Symbol, ...]:
    params = set(parameters)
    if variables is not None:
        return tuple(sym for sym in _normalize_symbols(variables) if sym not in params)
    return tuple(sorted(formula.free_symbols - params, key=lambda item: item.name))


def _conjuncts(expr: sp.Expr) -> tuple[sp.Expr, ...]:
    if expr is sp.true or expr == sp.true:
        return ()
    if isinstance(expr, sp.And):
        out: list[sp.Expr] = []
        for arg in expr.args:
            out.extend(_conjuncts(arg))
        return tuple(out)
    return (expr,)


def _interval_components_from_set(
    set_expr: sp.Set, var: sp.Symbol
) -> tuple[IntervalComponent, ...]:
    """Convert a one-dimensional SymPy set to exact interval components."""

    if set_expr is sp.S.EmptySet or set_expr == sp.S.EmptySet:
        return ()
    pieces = set_expr.args if isinstance(set_expr, sp.Union) else (set_expr,)
    components: list[IntervalComponent] = []
    for piece in pieces:
        if isinstance(piece, sp.Interval):
            components.append(
                IntervalComponent(
                    var,
                    piece.start,
                    piece.end,
                    not bool(piece.left_open),
                    not bool(piece.right_open),
                )
            )
        elif isinstance(piece, sp.FiniteSet):
            for point in sorted(piece, key=sp.default_sort_key):
                components.append(IntervalComponent(var, point, point, True, True))
        else:
            raise NotImplementedError(f"unsupported one-dimensional solution set piece: {piece!r}")
    components.sort(
        key=lambda comp: (sp.default_sort_key(comp.lower), sp.default_sort_key(comp.upper))
    )
    return tuple(components)


def _one_dim_components(expr: sp.Expr, var: sp.Symbol) -> tuple[IntervalComponent, ...] | None:
    """Return exact connected components for supported one-dimensional formulas."""

    if expr is sp.false or expr == sp.false:
        return ()
    if expr is sp.true or expr == sp.true:
        return (IntervalComponent(var, -sp.oo, sp.oo, False, False),)
    try:
        set_expr = expr.as_set()
        if isinstance(set_expr, sp.ConditionSet):
            return None
        return _interval_components_from_set(set_expr, var)
    except Exception:
        return None


def _components_formula(components: Sequence[IntervalComponent]) -> sp.Expr:
    if not components:
        return sp.false
    formulas = [component.as_formula() for component in components]
    return sp.Or(*formulas) if len(formulas) > 1 else formulas[0]


def _safe_call(default, func, *args, **kwargs):
    try:
        return func(*args, **kwargs)
    except Exception:
        return default


def _has_strict_atom(expr: sp.Expr) -> bool:
    if isinstance(expr, (sp.StrictLessThan, sp.StrictGreaterThan)):
        return True
    if isinstance(expr, (sp.And, sp.Or)):
        return any(_has_strict_atom(arg) for arg in expr.args)
    return False


def _one_dim_bounds(expr: sp.Expr, var: sp.Symbol) -> tuple[sp.Expr, sp.Expr] | None:
    """Extract conservative 1D interval bounds from a conjunction."""

    lo: sp.Expr = -sp.oo
    hi: sp.Expr = sp.oo
    try:
        reduced = sp.reduce_inequalities(list(_conjuncts(expr)), var)
    except Exception:
        reduced = expr
    atoms = _conjuncts(reduced)
    for atom in atoms:
        if not isinstance(atom, sp.core.relational.Relational):
            continue
        lhs, rhs = atom.lhs, atom.rhs
        if isinstance(atom, (sp.LessThan, sp.StrictLessThan)):
            if lhs == var and not rhs.has(var):
                hi = sp.Min(hi, rhs) if hi != sp.oo else rhs
            elif rhs == var and not lhs.has(var):
                lo = sp.Max(lo, lhs) if lo != -sp.oo else lhs
        elif isinstance(atom, (sp.GreaterThan, sp.StrictGreaterThan)):
            if lhs == var and not rhs.has(var):
                lo = sp.Max(lo, rhs) if lo != -sp.oo else rhs
            elif rhs == var and not lhs.has(var):
                hi = sp.Min(hi, lhs) if hi != sp.oo else lhs
        elif isinstance(atom, sp.Equality):
            if lhs == var and not rhs.has(var):
                lo = hi = rhs
            elif rhs == var and not lhs.has(var):
                lo = hi = lhs
    return (lo, hi)


def _fast_parameter_conditions(
    expr: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    parameters: tuple[sp.Symbol, ...],
) -> sp.Expr | None:
    """Fast parameter conditions for common one-variable polynomial atoms."""

    if len(variables) != 1 or not isinstance(expr, sp.core.relational.Relational):
        return None
    var = variables[0]
    try:
        poly = sp.Poly(sp.expand(expr.lhs - expr.rhs), var)
    except Exception:
        return None
    degree = poly.degree()
    if isinstance(expr, sp.Equality):
        if degree == 0:
            return sp.Eq(poly.as_expr(), 0)
        if degree == 1:
            return sp.Ne(poly.LC(), 0)
        if degree == 2:
            coeffs = poly.all_coeffs()
            a2, a1, a0 = coeffs
            return _safe_simplify_expr(a1**2 - 4 * a2 * a0 >= 0)
    return None


def _fast_solution_formula(
    expr: sp.Expr, variables: tuple[sp.Symbol, ...]
) -> tuple[sp.Expr, str, bool | None]:
    """Fast conservative formula/satisfiability path for common solves."""

    if expr is sp.true or expr == sp.true:
        return sp.true, "trivial", True
    if expr is sp.false or expr == sp.false:
        return sp.false, "trivial", False
    if len(variables) == 1:
        components = _one_dim_components(expr, variables[0])
        if components is not None:
            reduced = _components_formula(components)
            return reduced, "one_dimensional_components", bool(components)
        try:
            reduced = sp.reduce_inequalities(list(_conjuncts(expr)), variables[0])
            if reduced is not None:
                return (
                    reduced,
                    "sympy_reduce_inequalities",
                    reduced is not sp.false and reduced != sp.false,
                )
        except Exception:
            pass
    if len(variables) == 2:
        try:
            from ..implicit_utils import decompose_cylindrical_formula_to_vertical_bounds_2d

            cells = tuple(decompose_cylindrical_formula_to_vertical_bounds_2d(expr, variables))
            if cells:
                return expr, "vertical_bounds_2d", True
        except Exception:
            pass
    return expr, "cad", None


def _point_formula(point: Mapping[sp.Symbol, sp.Expr], variables: Sequence[sp.Symbol]) -> sp.Expr:
    """Return an exact conjunction describing one finite solution point."""

    atoms = [sp.Eq(var, sp.sympify(point[var])) for var in variables]
    return sp.And(*atoms) if atoms else sp.true


def _finite_points_formula(
    points: Sequence[Mapping[sp.Symbol, sp.Expr]], variables: Sequence[sp.Symbol]
) -> sp.Expr:
    """Return an exact finite-set formula from point assignments."""

    if not points:
        return sp.false
    formulas = [_point_formula(point, variables) for point in points]
    return sp.Or(*formulas) if len(formulas) > 1 else formulas[0]


def _try_rur_formula(
    formula: sp.Expr,
    variables: Sequence[sp.Symbol],
    *,
    max_solutions: int | None = None,
):
    """Try the exact RUR finite-system backend for a Boolean formula.

    Returns ``None`` when the formula is outside the supported finite equality
    fragment. A returned object with ``status == 'unknown'`` means RUR saw at
    least one unsupported branch, so callers must not use it as an UNSAT proof.
    """

    if not variables:
        return None
    try:
        result = solve_formula_with_rur(
            formula, tuple(variables), real=True, max_solutions=max_solutions
        )
    except Exception:
        return None
    if result is None:
        return None
    if str(result.status).lower() == "unknown":
        return None
    if result.partial and not result.assignments:
        return None
    return result


def _collect_solution_metadata(
    formula: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    *,
    strategy: str | None = None,
) -> dict[str, object]:
    """Best-effort lightweight region metadata for ``solve_semialgebraic``."""

    metadata: dict[str, object] = {
        "dimension": None,
        "bounded": None,
        "closed": None,
        "compact": None,
        "components": (),
        "cells": (),
        "cylindrical_solution": None,
        "connectivity": None,
    }
    if formula is sp.false or formula == sp.false:
        metadata.update({"dimension": None, "bounded": True, "closed": True, "compact": True})
        return metadata
    if formula is sp.true or formula == sp.true:
        metadata.update(
            {
                "dimension": len(variables),
                "bounded": len(variables) == 0,
                "closed": True,
                "compact": len(variables) == 0,
            }
        )
        return metadata

    # Exact interval/component metadata for one-dimensional systems.
    if len(variables) == 1:
        components = _one_dim_components(formula, variables[0])
        if components is not None:
            if not components:
                metadata.update(
                    {
                        "dimension": None,
                        "bounded": True,
                        "closed": True,
                        "compact": True,
                        "components": (),
                    }
                )
            else:
                dimension = max(component.dimension for component in components)
                bounded = all(component.bounded for component in components)
                closed = all(component.closed for component in components)
                metadata.update(
                    {
                        "dimension": dimension,
                        "bounded": bounded,
                        "closed": closed,
                        "compact": bounded and closed,
                        "components": components,
                    }
                )
            return metadata
        bds = _one_dim_bounds(formula, variables[0])
        if bds is not None:
            lo, hi = bds
            finite = lo != -sp.oo and hi != sp.oo
            closed = not _has_strict_atom(formula)
            metadata.update(
                {"dimension": 1, "bounded": finite, "closed": closed, "compact": finite and closed}
            )

    # Fast box/interval metadata, including ordinary 1D inequalities. Keep
    # box extraction and 2D cell extraction independent: disjoint unions are not
    # boxes, but they may still decompose into supported vertical cells.
    try:
        from ..implicit_utils import extract_symbolic_box_bounds

        box = extract_symbolic_box_bounds(formula, variables)
        if box is not None:
            finite = all(lo != -sp.oo and hi != sp.oo for _, lo, hi in box.limits)
            closed = not _has_strict_atom(formula)
            metadata.update(
                {
                    "dimension": len(variables),
                    "bounded": finite,
                    "closed": closed,
                    "compact": finite and closed,
                }
            )
    except Exception:
        pass

    if len(variables) == 2:
        try:
            from ..implicit_utils import decompose_cylindrical_formula_to_vertical_bounds_2d

            cells_all = tuple(
                decompose_cylindrical_formula_to_vertical_bounds_2d(formula, variables)
            )
            cells = (
                tuple(cell for cell in cells_all if getattr(cell, "dimension", 2) == 2) or cells_all
            )
            if cells:
                metadata["cells"] = cells
                metadata["dimension"] = 2
                finite_cells = all(
                    cell.x_interval[0] != -sp.oo
                    and cell.x_interval[1] != sp.oo
                    and all(lower != -sp.oo and upper != sp.oo for lower, upper in cell.y_bounds)
                    for cell in cells
                )
                metadata["bounded"] = finite_cells
                metadata["closed"] = not _has_strict_atom(formula)
                metadata["compact"] = finite_cells and bool(metadata["closed"])
        except Exception:
            try:
                from ..cad.cells import extract_vertical_bounds_from_cad_2d

                cells_all = tuple(
                    extract_vertical_bounds_from_cad_2d(
                        formula, variables, full_dimensional_only=False
                    )
                )
                cells = (
                    tuple(cell for cell in cells_all if getattr(cell, "dimension", 2) == 2)
                    or cells_all
                )
                if cells:
                    metadata["cells"] = cells
                    metadata["dimension"] = max(getattr(cell, "dimension", 2) for cell in cells)
                    finite_cells = all(
                        cell.x_interval[0] != -sp.oo
                        and cell.x_interval[1] != sp.oo
                        and all(
                            lower != -sp.oo and upper != sp.oo for lower, upper in cell.y_bounds
                        )
                        for cell in cells
                    )
                    metadata["bounded"] = finite_cells
                    metadata["closed"] = not _has_strict_atom(formula)
                    metadata["compact"] = finite_cells and bool(metadata["closed"])
            except Exception:
                pass

    # Full cylindrical CAD solution representation. This is deliberately a
    # metadata layer rather than a replacement for the simplified formula: it
    # preserves nested algebraic bounds in variable order for arbitrary CAD cells
    # when the complete CAD engine can extract them.

    has_nonlinear_vertical_cell = False
    if len(variables) >= 2 and metadata.get("cells"):
        y_var = variables[1]
        for atom in getattr(formula, "args", (formula,)):
            if getattr(atom, "is_Relational", False) and y_var in getattr(
                atom, "free_symbols", set()
            ):
                try:
                    if sp.Poly(atom.lhs - atom.rhs, y_var).degree() > 1:
                        has_nonlinear_vertical_cell = True
                        break
                except Exception:
                    pass
    if (
        len(variables) >= 2
        and metadata.get("cylindrical_solution") is None
        and not has_nonlinear_vertical_cell
    ):
        try:
            from ..cad.cells import (
                extract_cylindrical_solution,
                extract_explicit_cylindrical_solution,
            )

            cyl = extract_explicit_cylindrical_solution(formula, variables)
            if cyl is None:
                cyl = extract_cylindrical_solution(formula, variables, selected_only=True)
            metadata["cylindrical_solution"] = cyl
            if cyl.cells:
                metadata["dimension"] = cyl.dimension
                metadata["bounded"] = cyl.bounded
                if not metadata.get("cells"):
                    metadata["cells"] = cyl.cells
                try:
                    from ..connectivity import build_cad_adjacency_graph

                    connectivity = build_cad_adjacency_graph(cyl, formula=formula)
                    metadata["connectivity"] = connectivity
                    if connectivity.components:
                        metadata["components"] = connectivity.components
                except Exception:
                    pass
        except Exception:
            pass

    if metadata["dimension"] is None:
        equalities = [atom for atom in _conjuncts(formula) if isinstance(atom, sp.Equality)]
        if equalities:
            metadata["dimension"] = max(0, len(variables) - len(equalities))
        else:
            metadata["dimension"] = len(variables)
    return metadata


def _make_quantified_sentence(
    formula: sp.Expr, variables: Sequence[sp.Symbol]
) -> tuple[tuple[str, sp.Symbol], ...]:
    return tuple(("exists", var) for var in variables)


def _truth_from_qe_result(result) -> bool:
    if result.is_sentence:
        return bool(result.truth_value)
    simplified = _safe_simplify_expr(result.formula)
    if simplified is sp.true or simplified == sp.true:
        return True
    if simplified is sp.false or simplified == sp.false:
        return False
    # A non-sentence result means parameters escaped the requested variable set.
    # Treat satisfiability existentially over remaining free symbols.
    remaining = tuple(sorted(simplified.free_symbols, key=lambda item: item.name))
    if not remaining:
        return bool(simplified)
    return bool(
        qe_by_complete_cad(
            remaining, _make_quantified_sentence(simplified, remaining), parse_formula(simplified)
        ).truth_value
    )


def _validate_witness(
    formula: sp.Expr, point: Mapping[sp.Symbol, sp.Expr] | None, variables: Sequence[sp.Symbol]
) -> Mapping[sp.Symbol, sp.Expr] | None:
    """Return a normalized satisfying witness, or ``None`` if validation fails."""

    if point is None:
        return None
    normalized = {var: sp.sympify(point[var]) for var in variables if var in point}
    if len(normalized) != len(tuple(variables)):
        return None
    try:
        if satisfies_formula(formula, normalized):
            return normalized
    except Exception:
        pass
    return None


def _find_validated_witness(
    formula: sp.Expr, variables: Sequence[sp.Symbol], *, strategy: str | None = None
) -> Mapping[sp.Symbol, sp.Expr] | None:
    """Find one certified witness using the public sampling path."""

    try:
        points = sample_points(formula, variables, count=1, strategy=strategy or "auto")
    except Exception:
        points = ()
    for point in points:
        validated = _validate_witness(formula, point, variables)
        if validated is not None:
            return validated
    return None


def is_satisfiable(
    formula: FormulaLike | Iterable[FormulaLike],
    variables: Sequence[sp.Symbol | str] | None = None,
    *,
    domain: str = "reals",
    strategy: str | None = None,
    return_result: bool = False,
) -> bool | SatisfiabilityResult:
    """Return whether a real semialgebraic formula has a satisfying point.

    By default this preserves the historical boolean API. With
    ``return_result=True`` it returns a ``SatisfiabilityResult`` containing the
    normalized formula, variable order, backend method, and a validated witness
    when one is cheaply available.
    """

    if domain.lower() not in {"real", "reals", "r", "rr"}:
        raise NotImplementedError("is_satisfiable currently supports only the real domain")
    expr = _normalize_formula(formula)
    vars_ = _normalize_variables(variables, expr)
    method = "trivial"
    witness: Mapping[sp.Symbol, sp.Expr] | None = None
    if expr is sp.true or expr == sp.true:
        sat = True
        witness = {var: sp.Integer(0) for var in vars_}
    elif expr is sp.false or expr == sp.false:
        sat = False
    else:
        # First try exact finite-system dispatch. This proves SAT/UNSAT
        # for supported zero-dimensional equality branches without constructing
        # a full CAD of the ambient space.
        rur_result = _try_rur_formula(expr, vars_, max_solutions=1)
        if rur_result is not None and not rur_result.partial:
            sat = bool(rur_result.assignments)
            method = "rational_univariate"
            witness = dict(rur_result.assignments[0]) if rur_result.assignments else None
        else:
            # Try a validated witness before paying for full QE. This never proves
            # unsatisfiability, but it gives cheap structured results for common
            # full-dimensional feasible regions.
            witness = _find_validated_witness(expr, vars_, strategy=strategy)
            if witness is not None:
                sat = True
                method = "validated_sample"
            else:
                if len(vars_) == 1:
                    try:
                        reduced = sp.reduce_inequalities(list(_conjuncts(expr)), vars_[0])
                        comps = _one_dim_components(reduced, vars_[0])
                        if comps is not None:
                            sat = bool(comps)
                            method = "sympy_reduce_inequalities"
                            witness = {vars_[0]: comps[0].sample_point()} if comps else None
                        else:
                            raise ValueError("univariate reduction did not yield components")
                    except Exception:
                        result = qe_by_complete_cad(
                            vars_, _make_quantified_sentence(expr, vars_), parse_formula(expr)
                        )
                        sat = _truth_from_qe_result(result)
                        method = getattr(result, "method", "complete_cad_qe")
                else:
                    result = qe_by_complete_cad(
                        vars_, _make_quantified_sentence(expr, vars_), parse_formula(expr)
                    )
                    sat = _truth_from_qe_result(result)
                    method = getattr(result, "method", "complete_cad_qe")
                if sat and witness is None:
                    witness = _find_validated_witness(expr, vars_, strategy=strategy)
    if return_result:
        return SatisfiabilityResult(
            bool(sat),
            expr,
            vars_,
            witness=witness,
            method=method,
            diagnostics={"domain": domain, "strategy": strategy},
        )
    return bool(sat)


def is_tautology(
    formula: FormulaLike | Iterable[FormulaLike],
    variables: Sequence[sp.Symbol | str] | None = None,
    *,
    domain: str = "reals",
    strategy: str | None = None,
    return_result: bool = False,
) -> bool | TautologyResult:
    """Return whether a real semialgebraic formula is true for all variables.

    With ``return_result=True``, a false result includes a validated
    counterexample whenever the sampling layer can provide one.
    """

    expr = _normalize_formula(formula)
    vars_ = _normalize_variables(variables, expr)
    negated = sp.Not(expr)
    sat = is_satisfiable(negated, vars_, domain=domain, strategy=strategy, return_result=True)
    taut = not bool(sat)
    if return_result:
        return TautologyResult(
            taut,
            expr,
            vars_,
            counterexample=sat.witness if not taut else None,
            method=sat.method,
            diagnostics={"satisfiability": sat.diagnostics},
        )
    return taut


def implies(
    assumptions: FormulaLike | Iterable[FormulaLike],
    conclusion: FormulaLike,
    variables: Sequence[sp.Symbol | str] | None = None,
    *,
    domain: str = "reals",
    strategy: str | None = None,
    return_result: bool = False,
) -> bool | ImplicationResult:
    """Return whether ``assumptions`` imply ``conclusion`` over the reals.

    With ``return_result=True``, invalid implications include a validated
    counterexample satisfying the premise and falsifying the conclusion when
    available.
    """

    premise = _normalize_formula(assumptions)
    consequent = _normalize_formula(conclusion)
    universe = _normalize_variables(variables, sp.And(premise, consequent))
    counterexample_formula = sp.And(premise, sp.Not(consequent))
    sat = is_satisfiable(
        counterexample_formula, universe, domain=domain, strategy=strategy, return_result=True
    )
    valid = not bool(sat)
    if return_result:
        return ImplicationResult(
            valid,
            premise,
            consequent,
            universe,
            counterexample=sat.witness if not valid else None,
            method=sat.method,
            diagnostics={
                "counterexample_formula": sp.sstr(counterexample_formula),
                "satisfiability": sat.diagnostics,
            },
        )
    return valid


def equivalent(
    lhs: FormulaLike,
    rhs: FormulaLike,
    variables: Sequence[sp.Symbol | str] | None = None,
    *,
    domain: str = "reals",
    strategy: str | None = None,
    return_result: bool = False,
) -> bool | EquivalenceResult:
    """Return whether two semialgebraic formulas define the same real set.

    With ``return_result=True``, a false result includes a counterexample from
    the symmetric difference and, when it can be determined cheaply, the failed
    implication direction.
    """

    left = _normalize_formula(lhs)
    right = _normalize_formula(rhs)
    universe = _normalize_variables(variables, sp.And(left, right))
    if sp.sstr(left) == sp.sstr(right):
        if return_result:
            return EquivalenceResult(True, left, right, universe, method="syntactic")
        return True
    try:
        if sp.simplify_logic(sp.Xor(left, right)) is sp.false:
            if return_result:
                return EquivalenceResult(True, left, right, universe, method="logic_simplify")
            return True
    except (TypeError, ValueError, NotImplementedError, AttributeError):
        pass
    if len(universe) == 1:
        try:
            same_set = bool(left.as_set() == right.as_set())
            if same_set:
                if return_result:
                    return EquivalenceResult(True, left, right, universe, method="sympy_set")
                return True
        except (TypeError, ValueError, NotImplementedError, AttributeError):
            pass
    left_not_right = sp.And(left, sp.Not(right))
    right_not_left = sp.And(right, sp.Not(left))
    lnr = is_satisfiable(
        left_not_right, universe, domain=domain, strategy=strategy, return_result=True
    )
    rnl = is_satisfiable(
        right_not_left, universe, domain=domain, strategy=strategy, return_result=True
    )
    equiv = not bool(lnr) and not bool(rnl)
    failed_direction = None
    witness = None
    if bool(lnr) and bool(rnl):
        failed_direction = "both"
        witness = lnr.witness or rnl.witness
    elif bool(lnr):
        failed_direction = "lhs_implies_rhs"
        witness = lnr.witness
    elif bool(rnl):
        failed_direction = "rhs_implies_lhs"
        witness = rnl.witness
    if return_result:
        return EquivalenceResult(
            equiv,
            left,
            right,
            universe,
            counterexample=witness,
            failed_direction=failed_direction,
            method="symmetric_difference",
            diagnostics={"lhs_not_rhs": lnr.diagnostics, "rhs_not_lhs": rnl.diagnostics},
        )
    return equiv


class CellOutput(tuple):
    """Tuple-compatible cells view with a ``.cells`` alias."""

    @property
    def cells(self):
        return tuple(self)


def _select_solution_output(result: SemialgebraicSolution, output: str | None) -> object:
    """Return a selected public view of a solution result."""

    if output is None:
        return result
    key = output.lower().replace("-", "_")
    if key in {"result", "solution", "object"}:
        return result
    if key in {"formula", "constraints"}:
        return result.formula
    if key in {"reduced_formula", "reduced", "best_formula", "reduce"}:
        # Public selector should preserve the solved semialgebraic set, including
        # closed boundaries. Cell formulas may intentionally describe only
        # open full-dimensional interiors, so use the stored formula here.
        return result.formula
    if key in {"piecewise", "indicator", "indicator_piecewise"}:
        return result.as_piecewise()
    if key in {"sample", "one_sample"}:
        return result.sample
    if key in {"samples", "points"}:
        return result.samples
    if key in {"components", "component"}:
        return result.components
    if key in {"cells", "cell"}:
        return CellOutput(result.cells)
    if key in {"cylindrical", "cylindrical_solution", "cylindrical_cells", "cad_cells"}:
        return result.cylindrical_solution
    if key in {"connectivity", "adjacency", "roadmap", "roadmap_graph", "components_graph"}:
        return result.connectivity
    if key in {"plot_data", "discretization", "discretized", "mesh_data"}:
        return result.discretize()
    if key in {"diagnostics", "explain", "explanation"}:
        return result.explain()
    if key in {
        "parameter_strata",
        "parameter_decomposition",
        "strata",
        "piecewise_solution",
        "parameterized_solution",
    }:
        return result.parameter_decomposition
    if key in {"conditions", "parameter_conditions", "solvability_conditions"}:
        if result.parameter_conditions is not None:
            return result.parameter_conditions
        if not result.parameters:
            return sp.true if result.satisfiable else sp.false
        return sp.false
    if key in {"satisfiable", "nonempty"}:
        return result.satisfiable
    if key in {"empty", "unsatisfiable"}:
        return result.empty
    raise ValueError(f"unsupported solve_semialgebraic output selector: {output!r}")


def _add_standard_solver_diagnostics(
    diagnostics: dict[str, object],
    *,
    method: str,
    variables: Sequence[sp.Symbol],
    projection_order: Sequence[sp.Symbol | str] | None,
    domain_normalization: object | None,
    metadata: Mapping[str, object] | None = None,
    parameter_decomposition: object | None = None,
    solved: object | None = None,
) -> dict[str, object]:
    """Populate common explainability keys used by solve results."""

    metadata = metadata or {}
    diagnostics.setdefault("backend", method)
    diagnostics.setdefault("normalization_steps", ())
    diagnostics.setdefault("removed_redundant_constraints", ())
    diagnostics.setdefault("unsupported_features", ())
    diagnostics["requested_method"] = method
    diagnostics["variable_order"] = tuple(sp.sstr(v) for v in variables)
    diagnostics["projection_order"] = tuple(
        sp.sstr(v) for v in _normalize_symbols(projection_order)
    )
    diagnostics["used_interval_decomposition"] = (
        bool(metadata.get("components")) and len(tuple(variables)) == 1
    )
    diagnostics["used_cad"] = bool(
        metadata.get("cells")
        or metadata.get("cylindrical_solution")
        or metadata.get("connectivity")
    ) or method in {"cad", "qe", "cylindrical"}
    diagnostics["used_qe"] = method in {"cad", "qe"} or solved is not None
    diagnostics["used_cylindrical_solution"] = metadata.get("cylindrical_solution") is not None
    diagnostics["used_connectivity"] = metadata.get("connectivity") is not None
    diagnostics["used_parameter_decomposition"] = parameter_decomposition is not None
    rewrites = (
        tuple(getattr(domain_normalization, "rewrites", ()))
        if domain_normalization is not None
        else ()
    )
    constraints = (
        tuple(sp.sstr(c) for c in getattr(domain_normalization, "domain_constraints", ()))
        if domain_normalization is not None
        else ()
    )
    active_domain_norm = bool(rewrites or constraints)
    diagnostics["domain_normalized"] = active_domain_norm
    if domain_normalization is not None:
        diagnostics["domain_rewrites"] = rewrites
        diagnostics["domain_constraints"] = constraints
        if active_domain_norm:
            diagnostics["normalization_steps"] = tuple(
                diagnostics.get("normalization_steps", ())
            ) + ("domain-sensitive-constraint-normalization",)
    return diagnostics


def solve_semialgebraic(
    constraints: FormulaLike | Iterable[FormulaLike],
    variables: Sequence[sp.Symbol | str] | None = None,
    *,
    parameters: Sequence[sp.Symbol | str] | None = None,
    domain: str = "reals",
    count: int = 1,
    samples: int | str | None = None,
    sample_mode: str | None = None,
    strategy: str | None = None,
    method: str = "auto",
    variable_order: Sequence[sp.Symbol | str] | None = None,
    projection_order: Sequence[sp.Symbol | str] | None = None,
    normalize_domains: bool = True,
    return_formula: bool = False,
    output: str | None = None,
) -> SemialgebraicSolution | sp.Expr | tuple[object, ...] | bool | None:
    """Reduce, sample, and summarize a semialgebraic system over the reals.

    ``output`` may be used as a convenience selector for common views of the
    solution. The default ``None`` preserves structured result-object behavior.
    Component- and cell-aware sampling is available through
    ``samples="per_component"``, ``samples="per_cell"``, or the equivalent
    ``sample_mode`` keyword. Supported selectors are ``"formula"``, ``"reduced_formula"``,
    ``"piecewise"``, ``"samples"``, ``"components"``, ``"cells"``,
    ``"cylindrical"``, and ``"conditions"``. The ``"conditions"`` selector returns
    the parameter-space condition under which the system is solvable.

    The returned ``SemialgebraicSolution`` includes best-effort metadata such as simplified
    constraints, parameter conditions, dimension, boundedness, compactness,
    exact 1D components, and 2D vertical-bound cells when those analyses are
    supported. Unsupported metadata is reported as ``None`` or an empty tuple
    rather than being guessed.
    """

    if domain.lower() not in {"real", "reals", "r", "rr"}:
        raise NotImplementedError("solve_semialgebraic currently supports only the real domain")
    expr_original = _normalize_formula(constraints)
    params = _normalize_symbols(parameters)
    vars_initial = _normalize_solve_variables(variables, expr_original, params)
    if variable_order is not None:
        ordered = tuple(sym for sym in _normalize_symbols(variable_order) if sym not in set(params))
        remaining = tuple(sym for sym in vars_initial if sym not in set(ordered))
        vars_initial = ordered + remaining
    method_key = method.lower().replace("-", "_")
    if method_key not in {
        "auto",
        "interval",
        "linear",
        "rur",
        "cad",
        "qe",
        "cylindrical",
        "sampling",
    }:
        raise ValueError(f"unsupported solve_semialgebraic method: {method!r}")
    if projection_order is not None and variable_order is None:
        ordered = tuple(
            sym for sym in _normalize_symbols(projection_order) if sym not in set(params)
        )
        remaining = tuple(sym for sym in vars_initial if sym not in set(ordered))
        vars_initial = ordered + remaining

    domain_normalization = None
    expr = expr_original
    if normalize_domains:
        try:
            from .domain_solve import normalize_domain_sensitive_constraints

            domain_normalization = normalize_domain_sensitive_constraints(
                expr_original, vars_initial
            )
            expr = domain_normalization.formula
        except Exception:
            expr = expr_original
    vars_ = _normalize_solve_variables(tuple(vars_initial), expr, params)
    sample_count, resolved_sample_mode = _normalize_sample_request(count, samples, sample_mode)
    if method_key == "interval" and len(vars_) != 1:
        raise NotImplementedError("method='interval' supports exactly one solve variable")

    parameter_conditions: sp.Expr | None = None
    parameter_decomposition: object | None = None
    if params:
        condition_input = _conjuncts(expr)
        parameter_formula = condition_input[0] if len(condition_input) == 1 else expr
        parameter_conditions = _fast_parameter_conditions(parameter_formula, vars_, params)
        if parameter_conditions is None:
            from ..parameters import solvability_conditions

            parameter_conditions = solvability_conditions(
                parameter_formula, vars_, params, domain=domain
            )
        try:
            from ..parameter_stratification import parameterized_cylindrical_decomposition

            parameter_decomposition = parameterized_cylindrical_decomposition(
                expr, vars_, params, domain=domain, specialize_fibers=True
            )
            if parameter_conditions is None:
                parameter_conditions = parameter_decomposition.parameter_condition
        except Exception:
            parameter_decomposition = None

    if expr is sp.true or expr == sp.true:
        reduced = sp.true
        samples_out = tuple(
            {var: sp.Integer(0) for var in vars_} for _ in range(1 if sample_count else 0)
        )
        meta = _collect_solution_metadata(reduced, vars_, strategy=strategy)
        result = SemialgebraicSolution(
            reduced,
            vars_,
            samples_out,
            True,
            "trivial",
            solution_capability_diagnostics(expr),
            parameters=params,
            simplified_constraints=(),
            parameter_conditions=parameter_conditions,
            parameter_decomposition=parameter_decomposition,
            dimension=meta["dimension"],
            bounded=meta["bounded"],
            closed=meta["closed"],
            compact=meta["compact"],
            components=meta["components"],
            cells=meta["cells"],
            cylindrical_solution=meta.get("cylindrical_solution"),
            connectivity=meta.get("connectivity"),
        )
        return result.formula if return_formula else _select_solution_output(result, output)
    if expr is sp.false or expr == sp.false:
        meta = _collect_solution_metadata(sp.false, vars_, strategy=strategy)
        result = SemialgebraicSolution(
            sp.false,
            vars_,
            (),
            False,
            "trivial",
            solution_capability_diagnostics(expr),
            parameters=params,
            simplified_constraints=(),
            parameter_conditions=parameter_conditions
            if parameter_conditions is not None
            else sp.false,
            dimension=meta["dimension"],
            bounded=meta["bounded"],
            closed=meta["closed"],
            compact=meta["compact"],
            components=meta["components"],
            cells=meta["cells"],
            cylindrical_solution=meta.get("cylindrical_solution"),
            connectivity=meta.get("connectivity"),
        )
        return result.formula if return_formula else _select_solution_output(result, output)

    condition_keys = {"conditions", "parameter_conditions", "solvability_conditions"}
    if (
        output is not None
        and output.lower().replace("-", "_") in condition_keys
        and not return_formula
    ):
        return parameter_conditions if parameter_conditions is not None else sp.true

    if params:
        satisfiable = parameter_conditions is not sp.false and parameter_conditions != sp.false
        simplified_constraints = _conjuncts(expr)
        result = SemialgebraicSolution(
            expr if satisfiable else sp.false,
            vars_,
            (),
            bool(satisfiable),
            "parameter_conditions",
            _add_standard_solver_diagnostics(
                solution_capability_diagnostics(
                    expr,
                    selected_output=output,
                    selected_sample_mode=resolved_sample_mode,
                    requested_sample_count=sample_count,
                    has_parameter_conditions=parameter_conditions is not None,
                    has_param_decomp=parameter_decomposition is not None,
                ),
                method="parameter_conditions",
                variables=vars_,
                projection_order=projection_order,
                domain_normalization=domain_normalization,
                metadata={},
                parameter_decomposition=parameter_decomposition,
            ),
            parameters=params,
            simplified_constraints=simplified_constraints,
            parameter_conditions=parameter_conditions,
            parameter_decomposition=parameter_decomposition,
            dimension=None,
            bounded=None,
            closed=None,
            compact=None,
            components=(),
            cells=(),
            cylindrical_solution=None,
            connectivity=None,
        )
        return result.formula if return_formula else _select_solution_output(result, output)

    selected_method: str | None = None
    rur_result = None
    if method_key in {"auto", "rur"}:
        rur_result = _try_rur_formula(
            expr, vars_, max_solutions=sample_count if sample_count else None
        )
        if method_key == "rur" and rur_result is None:
            raise NotImplementedError(
                "method='rur' supports finite zero-dimensional equality branches only"
            )
    if rur_result is not None and not rur_result.partial:
        assignments = tuple(dict(point) for point in rur_result.assignments)
        reduced = _finite_points_formula(assignments, vars_) if assignments else sp.false
        satisfiable = bool(assignments)
        selected_method = "rational_univariate"
        fast_method = selected_method
        fast_satisfiable = satisfiable
        solved = None
        explicit_cyl = None
    elif method_key in {"cad", "qe", "cylindrical", "sampling"}:
        reduced, fast_method, fast_satisfiable = expr, method_key, None
        solved = None
    else:
        reduced, fast_method, fast_satisfiable = _fast_solution_formula(expr, vars_)
        solved = None
    explicit_cyl = None
    if selected_method is None and fast_satisfiable is None and len(vars_) > 1:
        try:
            from ..cad.cells import extract_explicit_cylindrical_solution

            explicit_cyl = extract_explicit_cylindrical_solution(expr, vars_)
        except Exception:
            explicit_cyl = None
    if selected_method is None and fast_satisfiable is None and explicit_cyl is not None:
        reduced = expr
        satisfiable = True
        selected_method = "explicit_cylindrical_bounds"
    elif selected_method is None and fast_satisfiable is None:
        parsed = ParsedPrenexFormula(vars_, (), parse_formula(expr), expr)
        solved = reduce_formula(parsed, domain=domain, return_result=True, strategy=strategy)
        reduced = _safe_simplify_expr(solved.result)
        try:
            satisfiable = is_satisfiable(reduced, vars_, domain=domain, strategy=strategy)
        except PolynomialError:
            # Reconstructed CAD formulas may contain algebraic boundary functions.
            # A non-false CAD reconstruction already carries selected cells, so use
            # it as the satisfiability witness when polynomial parsing is not available.
            satisfiable = reduced is not sp.false and reduced != sp.false
        selected_method = getattr(solved, "method", "cad")
    elif selected_method is None:
        satisfiable = bool(fast_satisfiable)
        selected_method = fast_method

    selected_method = selected_method or fast_method
    simplified_constraints: tuple[sp.Expr, ...] = _conjuncts(reduced)
    # Record a simplified constraint tuple without making full redundancy removal
    # part of the critical solve path. The dedicated ``simplify_system`` API
    # remains available for heavier semantic cleanup.

    final_formula = sp.false if not satisfiable else reduced
    meta = _collect_solution_metadata(final_formula, vars_, strategy=strategy)
    if satisfiable:
        samples_out = _collect_structural_samples(
            final_formula,
            expr,
            vars_,
            meta,
            count=sample_count,
            mode=resolved_sample_mode,
            strategy=strategy,
        )
    else:
        samples_out = ()
    diagnostics = dict(getattr(solved, "metadata", {}) or {}) if solved is not None else {}
    diagnostics.update(solution_capability_diagnostics(expr))
    diagnostics["selected_output"] = output
    diagnostics["selected_sample_mode"] = resolved_sample_mode
    diagnostics["requested_sample_count"] = sample_count
    diagnostics["structural_sample_count"] = len(samples_out)
    diagnostics["simplified_constraint_count"] = len(simplified_constraints)
    diagnostics["has_parameter_conditions"] = parameter_conditions is not None
    diagnostics["has_parameter_decomposition"] = parameter_decomposition is not None
    diagnostics["used_rur"] = selected_method == "rational_univariate"
    if selected_method == "rational_univariate" and rur_result is not None:
        diagnostics["rur_solved_branches"] = rur_result.solved_branches
        diagnostics["rur_skipped_branches"] = rur_result.skipped_branches
        diagnostics["rur_notes"] = tuple(rur_result.notes)
    diagnostics = _add_standard_solver_diagnostics(
        diagnostics,
        method=method_key,
        variables=vars_,
        projection_order=projection_order,
        domain_normalization=domain_normalization,
        metadata=meta,
        parameter_decomposition=parameter_decomposition,
        solved=solved,
    )
    result = SemialgebraicSolution(
        final_formula,
        vars_,
        samples_out,
        satisfiable,
        selected_method,
        diagnostics,
        parameters=params,
        simplified_constraints=simplified_constraints,
        parameter_conditions=parameter_conditions,
        parameter_decomposition=parameter_decomposition,
        dimension=meta["dimension"],
        bounded=meta["bounded"],
        closed=meta["closed"],
        compact=meta["compact"],
        components=meta["components"],
        cells=meta["cells"],
        cylindrical_solution=meta.get("cylindrical_solution"),
        connectivity=meta.get("connectivity"),
    )
    return result.formula if return_formula else _select_solution_output(result, output)


def canonicalize_one_dimensional_formula(
    formula: FormulaLike | Iterable[FormulaLike], variable: sp.Symbol | str
) -> sp.Expr:
    """Return a canonical interval-union formula for supported 1D systems."""

    expr = _normalize_formula(formula)
    var = _as_real_symbol(variable)
    components = _one_dim_components(expr, var)
    if components is None:
        try:
            reduced = sp.reduce_inequalities(list(_conjuncts(expr)), var)
            components = _one_dim_components(reduced, var)
        except Exception:
            components = None
    if components is None:
        return _safe_simplify_expr(expr)
    return _components_formula(components)


__all__ = [
    "IntervalComponent",
    "SemialgebraicSolution",
    "EquivalenceResult",
    "ImplicationResult",
    "TautologyResult",
    "SatisfiabilityResult",
    "equivalent",
    "implies",
    "is_satisfiable",
    "is_tautology",
    "canonicalize_one_dimensional_formula",
    "solve_semialgebraic",
]
