from __future__ import annotations

from collections.abc import Sequence

import sympy as sp

from ...status import CoverageStatus
from .atoms import (
    _and,
    _canonical_atom,
    _iter_atoms,
    _or,
    _polynomial_degree,
    _to_negation_normal_form,
)
from .substitution import (
    _root_candidates_from_polynomial,
    _unique_polynomials,
    substitute_infinity,
    substitute_perturbed_quadratic_root,
    substitute_quadratic_root,
)
from .types import (
    QuadraticVirtualSubstitutionResult,
    VirtualSubstitutionError,
    VirtualSubstitutionQEResult,
)


def _normalize_quantifier_name(name: str) -> str:
    lowered = name.lower()
    if lowered not in {"exists", "forall"}:
        raise VirtualSubstitutionError(f"unsupported quantifier: {name!r}")
    return lowered


def _free_symbols_in_expr(expr: sp.Expr) -> tuple[sp.Symbol, ...]:
    return tuple(sorted(getattr(expr, "free_symbols", set()), key=lambda sym: sym.name))


def can_use_quadratic_vs(formula: sp.Expr, variable: sp.Symbol) -> bool:
    """Return True when every atom is polynomial of degree at most two in variable."""

    try:
        normalized = _to_negation_normal_form(formula)
        for atom in _iter_atoms(normalized):
            polynomial, _ = _canonical_atom(atom)
            if _polynomial_degree(polynomial, variable) > 2:
                return False
        return True
    except VirtualSubstitutionError:
        return False


def try_quadratic_virtual_substitution_qe(
    vars_: Sequence[sp.Symbol],
    quantifiers: Sequence[tuple[str, sp.Symbol]],
    matrix: sp.Expr,
    *,
    full: bool = True,
    max_growth_factor: int = 2,
) -> VirtualSubstitutionQEResult | None:
    """Try a planner-level quadratic virtual-substitution pass before CAD.

    Existential blocks are eliminated directly. Pure universal blocks are
    eliminated by the identity ``ForAll[x, phi] == Not[Exists[x, Not[phi]]]``.
    Mixed alternating prefixes are intentionally left to CAD because correct
    prefix-order bookkeeping needs a full quantified-formula transformer rather
    than a local prepass.
    """

    normalized_quantifiers = tuple((_normalize_quantifier_name(q), sym) for q, sym in quantifiers)
    if not normalized_quantifiers:
        return None

    quantifier_names = {q for q, _ in normalized_quantifiers}
    if quantifier_names == {"forall"}:
        dual = try_quadratic_virtual_substitution_qe(
            vars_,
            tuple(("exists", sym) for _, sym in normalized_quantifiers),
            _to_negation_normal_form(matrix, negate=True),
            full=full,
            max_growth_factor=max_growth_factor,
        )
        if dual is None or dual.remaining_quantifiers:
            return None
        formula = _to_negation_normal_form(dual.formula, negate=True)
        pass
        result_symbols = _free_symbols_in_expr(formula)
        is_sentence = not result_symbols
        truth_value = None
        if is_sentence:
            simplified = sp.simplify(formula)
            if simplified == sp.true or simplified is sp.true:
                truth_value = True
                formula = sp.true
            elif simplified == sp.false or simplified is sp.false:
                truth_value = False
                formula = sp.false
        return VirtualSubstitutionQEResult(
            formula=formula,
            variables=tuple(vars_),
            free_variables=tuple(
                sym for sym in result_symbols if sym not in {v for _, v in normalized_quantifiers}
            ),
            quantified_variables=(),
            remaining_quantifiers=(),
            eliminated_variables=dual.eliminated_variables,
            is_sentence=is_sentence,
            truth_value=truth_value,
            status="complete",
            notes=(
                "used universal/existential duality for quadratic virtual substitution",
                *dual.notes,
            ),
        )
    if quantifier_names != {"exists"}:
        if full:
            return None

    current_formula = _to_negation_normal_form(matrix)
    current_size = max(1, int(sp.count_ops(current_formula, visual=False)))
    eliminated: list[sp.Symbol] = []
    remaining: list[tuple[str, sp.Symbol]] = []
    notes: list[str] = []

    for quantifier, variable in normalized_quantifiers:
        if quantifier != "exists":
            remaining.append((quantifier, variable))
            notes.append(f"left {quantifier} variable {variable} for CAD")
            continue
        if variable not in getattr(current_formula, "free_symbols", set()):
            eliminated.append(variable)
            notes.append(f"removed vacuous existential variable {variable}")
            continue
        if not can_use_quadratic_vs(current_formula, variable):
            remaining.append((quantifier, variable))
            notes.append(
                f"variable {variable} is outside the quadratic virtual-substitution fragment"
            )
            continue
        try:
            step = eliminate_exists_quadratic_variable(current_formula, variable, simplify=False)
        except VirtualSubstitutionError as exc:
            remaining.append((quantifier, variable))
            notes.append(f"virtual substitution failed for {variable}: {exc}")
            continue
        next_formula = step.formula
        next_size = max(1, int(sp.count_ops(next_formula, visual=False)))
        if not full and next_size > max_growth_factor * current_size:
            remaining.append((quantifier, variable))
            notes.append(f"skipped {variable}; formula growth exceeded limit")
            continue
        current_formula = next_formula
        current_size = next_size
        eliminated.append(variable)
        notes.append(
            f"eliminated existential variable {variable} by quadratic virtual substitution"
        )

    if not eliminated:
        return None
    if full and remaining:
        return None

    quantified_set = {var for _, var in normalized_quantifiers}
    free_variables = tuple(sym for sym in tuple(vars_) if sym not in quantified_set)
    result_symbols = _free_symbols_in_expr(current_formula)
    remaining_vars = {sym for _, sym in remaining}
    free_variables = tuple(sym for sym in free_variables if sym in result_symbols) + tuple(
        sym
        for sym in result_symbols
        if sym not in set(free_variables) and sym not in remaining_vars
    )
    quantified_variables = tuple(sym for _, sym in remaining)
    is_sentence = not result_symbols and not remaining
    truth_value = None
    formula_out = current_formula
    pass
    if is_sentence:
        simplified = sp.simplify(formula_out)
        if simplified == sp.true or simplified is sp.true:
            truth_value = True
            formula_out = sp.true
        elif simplified == sp.false or simplified is sp.false:
            truth_value = False
            formula_out = sp.false

    return VirtualSubstitutionQEResult(
        formula=formula_out,
        variables=tuple(vars_),
        free_variables=free_variables,
        quantified_variables=quantified_variables,
        remaining_quantifiers=tuple(remaining),
        eliminated_variables=tuple(eliminated),
        is_sentence=is_sentence,
        truth_value=truth_value,
        status=CoverageStatus.COMPLETE if not remaining else CoverageStatus.PARTIAL,
        notes=tuple(notes),
    )


def eliminate_quadratic_variable(
    formula: sp.Expr,
    variable: sp.Symbol,
    *,
    simplify: bool = True,
    return_result: bool = False,
) -> sp.Expr | QuadraticVirtualSubstitutionResult:
    """Eliminate one existential real variable from a quadratic formula.

    The input must be a quantifier-free SymPy Boolean formula built from
    polynomial relations, and every atom must have degree at most two in
    ``variable``. The returned formula is quantifier-free and does not contain
    ``variable``.
    """

    nnf_formula = _to_negation_normal_form(formula)
    if variable not in getattr(nnf_formula, "free_symbols", set()):
        result_formula = nnf_formula
    else:
        candidates: list[sp.Expr] = [substitute_infinity(nnf_formula, variable, -1)]
        for polynomial in _unique_polynomials(nnf_formula, variable):
            for guard, point in _root_candidates_from_polynomial(polynomial, variable):
                exact = substitute_quadratic_root(nnf_formula, variable, point)
                right = substitute_perturbed_quadratic_root(nnf_formula, variable, point, 1)
                candidates.append(_and(guard, exact))
                candidates.append(_and(guard, right))
        result_formula = _or(*candidates)

    if simplify:
        result_formula = (
            sp.factor(result_formula)
            if not isinstance(result_formula, (sp.And, sp.Or))
            else result_formula
        )
    if return_result:
        return QuadraticVirtualSubstitutionResult(
            formula=result_formula, eliminated_variable=variable
        )
    return result_formula


def eliminate_exists_quadratic_variable(
    formula: sp.Expr,
    variable: sp.Symbol,
    *,
    simplify: bool = True,
) -> QuadraticVirtualSubstitutionResult:
    """Return a structured result for one existential quadratic elimination."""

    return eliminate_quadratic_variable(formula, variable, simplify=simplify, return_result=True)  # type: ignore[return-value]
