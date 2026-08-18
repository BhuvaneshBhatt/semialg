from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence

import sympy as sp

from ...algebraic.rational_univariate import RationalUnivariateError, sign_of_algebraic_expression
from ...status import SolverStatus
from .atoms import _polynomial_degree, _to_negation_normal_form
from .eliminate import (
    _free_symbols_in_expr,
    _normalize_quantifier_name,
    can_use_quadratic_vs,
    eliminate_quadratic_variable,
)
from .substitution import _unique_polynomials
from .types import VirtualSubstitutionError, VirtualSubstitutionWitnessResult


def _truth_of_formula_at_assignment(
    formula: sp.Expr, assignment: Mapping[sp.Symbol, sp.Expr]
) -> bool | None:
    substituted = formula.subs(dict(assignment))
    try:
        simplified = sp.simplify(substituted)
    except (TypeError, ValueError, sp.SympifyError):
        simplified = substituted
    if simplified == sp.true or simplified is sp.true:
        return True
    if simplified == sp.false or simplified is sp.false:
        return False
    try:
        return bool(simplified)
    except TypeError:
        try:
            return bool(sp.N(simplified, 50))
        except (TypeError, ValueError, sp.SympifyError):
            return None


def _real_roots_of_low_degree_polynomial(
    polynomial: sp.Expr, variable: sp.Symbol
) -> tuple[sp.Expr, ...]:
    polynomial = sp.expand(polynomial)
    if polynomial == 0 or variable not in polynomial.free_symbols:
        return ()
    degree = _polynomial_degree(polynomial, variable)
    if degree <= 0:
        return ()
    if degree > 2:
        raise VirtualSubstitutionError(
            f"degree {degree} in {variable} exceeds the witness reconstruction fragment"
        )
    roots = sp.solve(sp.Eq(polynomial, 0), variable)
    real_roots: list[sp.Expr] = []
    for root in roots:
        root = sp.simplify(root)
        is_real = root.is_real
        if is_real is False:
            continue
        if is_real is True:
            real_roots.append(root)
            continue
        try:
            imag_part = abs(complex(sp.N(root, 50)).imag)
        except (TypeError, ValueError, sp.SympifyError):
            continue
        if imag_part < 1e-40:
            real_roots.append(root)
    return tuple(real_roots)


def _dedupe_by_expression(values: Iterable[sp.Expr]) -> tuple[sp.Expr, ...]:
    out: list[sp.Expr] = []
    seen: set[str] = set()
    for value in values:
        value = sp.simplify(value)
        key = sp.sstr(value)
        if key not in seen:
            out.append(value)
            seen.add(key)
    return tuple(out)


def _compare_real_algebraic_values(left: sp.Expr, right: sp.Expr) -> int:
    difference = sp.cancel(left - right)
    try:
        return sign_of_algebraic_expression(difference)
    except RationalUnivariateError:
        numeric_difference = sp.N(difference, 120)
        if numeric_difference > 0:
            return 1
        if numeric_difference < 0:
            return -1
        return 0


def _insert_ordered_real_value(values: list[sp.Expr], value: sp.Expr) -> None:
    for index, existing in enumerate(values):
        comparison = _compare_real_algebraic_values(value, existing)
        if comparison == 0:
            return
        if comparison < 0:
            values.insert(index, value)
            return
    values.append(value)


def _ordered_real_values(values: Sequence[sp.Expr]) -> list[sp.Expr]:
    ordered: list[sp.Expr] = []
    for value in values:
        _insert_ordered_real_value(ordered, value)
    return ordered


def _formula_samples(formula: sp.Expr, variable: sp.Symbol) -> tuple[sp.Expr, ...]:
    normalized = _to_negation_normal_form(formula)
    roots: list[sp.Expr] = []
    for polynomial in _unique_polynomials(normalized, variable):
        roots.extend(_real_roots_of_low_degree_polynomial(polynomial, variable))
    roots = list(_ordered_real_values(_dedupe_by_expression(roots)))

    candidates: list[sp.Expr] = []
    candidates.extend(roots)
    if not roots:
        candidates.append(sp.Integer(0))
        return _dedupe_by_expression(candidates)

    candidates.append(sp.simplify(roots[0] - 1))
    for left, right in zip(roots, roots[1:], strict=False):
        if sp.simplify(left - right) != 0:
            candidates.append(sp.simplify((left + right) / 2))
    candidates.append(sp.simplify(roots[-1] + 1))
    return _dedupe_by_expression(candidates)


def reconstruct_vs_value(
    formula: sp.Expr,
    variable: sp.Symbol,
    known_values: Mapping[sp.Symbol, sp.Expr],
) -> sp.Expr | None:
    """Find one concrete value of ``variable`` satisfying a quadratic stage.

    ``formula`` is the stage formula before ``variable`` was eliminated. All
    other free variables should already have values in ``known_values``. The
    reconstruction samples exact boundary roots and one representative from
    each open interval determined by those roots.
    """

    partially_evaluated = _to_negation_normal_form(formula.subs(dict(known_values)))
    if variable not in getattr(partially_evaluated, "free_symbols", set()):
        truth = _truth_of_formula_at_assignment(partially_evaluated, {})
        return sp.Integer(0) if truth else None
    for candidate in _formula_samples(partially_evaluated, variable):
        truth = _truth_of_formula_at_assignment(partially_evaluated, {variable: candidate})
        if truth is True:
            return sp.simplify(candidate)
    return None


def try_quadratic_virtual_substitution_witness(
    vars_: Sequence[sp.Symbol],
    quantifiers: Sequence[tuple[str, sp.Symbol]],
    matrix: sp.Expr,
    base_instance_finder,
    *,
    full: bool = True,
) -> VirtualSubstitutionWitnessResult | None:
    """Try to find a witness using quadratic virtual substitution.

    ``base_instance_finder`` is called as ``finder(reduced_formula, variables)``
    after VS has eliminated eligible existential variables. It should return a
    mapping for the remaining free variables, or ``None`` when no base witness
    is found.
    """

    normalized_quantifiers = tuple((_normalize_quantifier_name(q), sym) for q, sym in quantifiers)
    if not normalized_quantifiers:
        return None
    if any(quantifier != "exists" for quantifier, _ in normalized_quantifiers):
        return None

    current_formula = _to_negation_normal_form(matrix)
    stages: list[tuple[sp.Symbol, sp.Expr]] = []
    eliminated: list[sp.Symbol] = []
    notes: list[str] = []

    for _, variable in normalized_quantifiers:
        if variable not in getattr(current_formula, "free_symbols", set()):
            stages.append((variable, current_formula))
            eliminated.append(variable)
            notes.append(f"assigned vacuous existential variable {variable} during reconstruction")
            continue
        if not can_use_quadratic_vs(current_formula, variable):
            if full:
                return None
            notes.append(f"left {variable} unreconstructed; outside quadratic fragment")
            continue
        stages.append((variable, current_formula))
        try:
            current_formula = eliminate_quadratic_variable(
                current_formula, variable, simplify=False
            )
        except VirtualSubstitutionError:
            return None
        eliminated.append(variable)
        notes.append(
            f"eliminated existential variable {variable} and recorded reconstruction stage"
        )

    if not eliminated:
        return None
    remaining_symbols = tuple(
        sym for sym in _free_symbols_in_expr(current_formula) if sym not in set(eliminated)
    )
    base_instance = base_instance_finder(current_formula, remaining_symbols)
    if base_instance is None:
        return VirtualSubstitutionWitnessResult(
            instance=None,
            eliminated_variables=tuple(eliminated),
            reduced_formula=current_formula,
            reduced_variables=remaining_symbols,
            status=SolverStatus.UNSAT,
            notes=tuple(notes) + ("no witness found for reduced formula",),
        )

    instance: dict[sp.Symbol, sp.Expr] = dict(base_instance)
    for variable, stage_formula in reversed(stages):
        value = reconstruct_vs_value(stage_formula, variable, instance)
        if value is None:
            return VirtualSubstitutionWitnessResult(
                instance=None,
                eliminated_variables=tuple(eliminated),
                reduced_formula=current_formula,
                reduced_variables=remaining_symbols,
                status=SolverStatus.ERROR,
                notes=tuple(notes) + (f"could not reconstruct {variable}",),
            )
        instance[variable] = value

    ordered_instance = {sym: sp.simplify(instance[sym]) for sym in tuple(vars_) if sym in instance}
    if _truth_of_formula_at_assignment(matrix, ordered_instance) is not True:
        return VirtualSubstitutionWitnessResult(
            instance=None,
            eliminated_variables=tuple(eliminated),
            reduced_formula=current_formula,
            reduced_variables=remaining_symbols,
            status=SolverStatus.ERROR,
            notes=tuple(notes)
            + ("reconstructed witness did not validate against original formula",),
        )
    return VirtualSubstitutionWitnessResult(
        instance=ordered_instance,
        eliminated_variables=tuple(eliminated),
        reduced_formula=current_formula,
        reduced_variables=remaining_symbols,
        status=SolverStatus.SAT,
        notes=tuple(notes),
    )
