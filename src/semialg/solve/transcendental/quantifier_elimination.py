from __future__ import annotations

from dataclasses import dataclass, field

import sympy as sp
from sympy import S

from .roots import decomp_univar_inequality, isolate_univar_roots
from .state import TransProblemState


@dataclass(frozen=True)
class QuantElimResult:
    variable: sp.Symbol
    quantifier: str
    resulting_formula: sp.Expr
    complete: bool = False
    method: str = "transcendental_real_quantifier_elimination"
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

    # Equation-driven existential reduction.
    root_info = isolate_univar_roots(simplified, variable, domain=S.Reals)
    if root_info.complete and root_info.roots:
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
                metadata={"witness_roots": tuple(subs_truths)},
            )

    dec = decomp_univar_inequality(simplified, variable, domain=S.Reals)
    if dec.true_intervals or dec.true_points:
        return QuantElimResult(
            variable,
            "exists",
            sp.true,
            complete=bool(dec.true_intervals),
            method="exists_via_interval_decomposition",
            metadata={"intervals": dec.true_intervals, "points": dec.true_points},
        )
    if dec.method.startswith("certified_"):
        return QuantElimResult(
            variable,
            "exists",
            sp.false,
            complete=False,
            method="exists_no_witness_in_certified_window",
            metadata={"support_points": dec.support_points},
        )
    return None


def _eliminate_forall_real(state: TransProblemState, variable: sp.Symbol) -> QuantElimResult | None:
    formula = extract_cand_form(state, variable)
    dec = decomp_univar_inequality(formula, variable, domain=S.Reals)
    if dec.method.startswith("certified_"):
        # If decomposition covers a whole bounded fundamental/certified window and every tested interval is true,
        # we can certify only on that window. This is not global completeness, so we mark partial.
        if dec.true_intervals and len(dec.true_intervals) == 1:
            a, b = dec.true_intervals[0]
            return QuantElimResult(
                variable,
                "forall",
                sp.false,
                complete=False,
                method="forall_partial_window_only",
                metadata={"covered_interval": (a, b)},
            )
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
