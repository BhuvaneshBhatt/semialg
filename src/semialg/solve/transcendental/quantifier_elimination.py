from __future__ import annotations

from dataclasses import dataclass, field

import sympy as sp
from sympy import S

from .roots import decomp_univar_inequality, isolate_univar_roots
from .semantics import ResultSemantics
from .state import TransProblemState


@dataclass(frozen=True)
class QuantElimResult:
    variable: sp.Symbol
    quantifier: str
    resulting_formula: sp.Expr
    complete: bool = False
    method: str = "transcendental_real_quantifier_elimination"
    result_semantics: ResultSemantics = ResultSemantics.UNKNOWN
    validity_window: tuple[sp.Expr, sp.Expr] | None = None
    metadata: dict = field(default_factory=dict)


def extract_cand_form(state: TransProblemState, variable: sp.Symbol) -> sp.Expr:
    return sp.simplify(state.formula)


def _eliminate_exists_real(state: TransProblemState, variable: sp.Symbol) -> QuantElimResult | None:
    formula = extract_cand_form(state, variable)
    try:
        simplified = sp.simplify(sp.logic.boolalg.eliminate_implications(formula))
    except Exception:
        simplified = formula

    # First try direct logic-level existential reduction for purely univariate formulas.
    try:
        if simplified.free_symbols <= {variable}:
            sat = sp.satisfiable(simplified, use_lra_theory=True)
            if sat is False:
                return QuantElimResult(
                    variable,
                    "exists",
                    sp.false,
                    complete=True,
                    method="satisfiable_exists",
                    metadata={"reason": "unsat_univariate"},
                )
    except Exception:
        pass

    # Equation-driven existential reduction.  Do not feed inequalities or
    # arbitrary Boolean formulas into scalar root isolation.
    root_info = None
    if isinstance(simplified, sp.Equality):
        root_info = isolate_univar_roots(simplified, variable, domain=S.Reals)
    if root_info is not None and root_info.complete and root_info.roots:
        subs_truths = []
        for root in root_info.roots:
            try:
                val = bool(sp.simplify(simplified.subs(variable, root)))
            except Exception:
                val = False
            if val:
                subs_truths.append(root)
        if subs_truths:
            return QuantElimResult(
                variable,
                "exists",
                sp.true,
                complete=True,
                method="complete_exists_via_roots",
                result_semantics=ResultSemantics.EXACT,
                metadata={"witness_roots": tuple(subs_truths)},
            )

    dec = decomp_univar_inequality(simplified, variable, domain=S.Reals)
    if dec.true_intervals or dec.true_points:
        return QuantElimResult(
            variable,
            "exists",
            sp.true,
            complete=False,
            method="exists_via_numerical_interval_witness",
            result_semantics=ResultSemantics.WITNESS_SUBSET,
            validity_window=dec.validity_window,
            metadata={"intervals": dec.true_intervals, "points": dec.true_points},
        )
    if dec.method.startswith("numerical_"):
        return QuantElimResult(
            variable,
            "exists",
            sp.false,
            complete=False,
            method="exists_no_witness_in_numerical_window",
            result_semantics=ResultSemantics.WINDOW_NO_WITNESS,
            validity_window=dec.validity_window,
            metadata={"support_points": dec.support_points},
        )
    return None


def _eliminate_forall_real(state: TransProblemState, variable: sp.Symbol) -> QuantElimResult | None:
    formula = extract_cand_form(state, variable)
    dec = decomp_univar_inequality(formula, variable, domain=S.Reals)
    # A finite numerical window can neither prove nor refute a global universal
    # statement.  Return no reduction rather than a misleading Boolean formula.
    if dec.method.startswith("numerical_") or dec.method.startswith("empty_numerical_"):
        return None
    return None


def eliminate_lead_block(state: TransProblemState) -> QuantElimResult | None:
    if not state.quantifier_blocks:
        return None
    block = state.quantifier_blocks[0]
    if len(block.variables) != 1:
        return None
    variable = block.variables[0]
    if state.default_domain != state.default_domain.REALS:
        return None
    if state.formula.free_symbols - {variable}:
        # Current implementation handles only genuinely univariate quantified transcendental reductions.
        return None
    if block.quantifier == "exists":
        return _eliminate_exists_real(state, variable)
    if block.quantifier == "forall":
        return _eliminate_forall_real(state, variable)
    return None


__all__ = ["QuantElimResult", "eliminate_lead_block"]
