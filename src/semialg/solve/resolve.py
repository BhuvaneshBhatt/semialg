from __future__ import annotations

import sympy as sp

from ..formula import ParsedPrenexFormula, parse_quant_form_text
from ..partial.qe import lazy_resolve_formula
from .domains import SolveDomain, normalize_domain
from .preprocess import semialgebraicize
from .reduce import reduce_formula
from .result import SolveResult


def _is_real_sentence(parsed: ParsedPrenexFormula) -> bool:
    quantified = {sym for _, sym in parsed.quantifiers}
    return set(parsed.vars).issubset(quantified) and bool(parsed.quantifiers)


def resolve_formula(
    parsed,
    config=None,
    *,
    domain: str | SolveDomain | None = None,
    assumptions=None,
    return_result: bool = False,
    strategy: str | None = "lazy",
):
    dom = normalize_domain(domain)
    strategy_name = (strategy or "lazy").lower()
    if (
        assumptions is None
        and dom is SolveDomain.REALS
        and strategy_name in {"lazy", "partial", "auto", "planner"}
        and _is_real_sentence(parsed)
    ):
        lazy = lazy_resolve_formula(parsed.vars, parsed.quantifiers, parsed.matrix)
        result = SolveResult(
            method="partial_cad_resolve",
            domain=dom,
            result=bool(lazy.truth_value),
            metadata={
                "lazy_result": lazy,
                "stats": lazy.stats,
                "witness": lazy.witness,
                "counterexample": lazy.counterexample,
            },
        )
        return result if return_result else result.result
    return reduce_formula(
        parsed,
        config=config,
        domain=domain,
        assumptions=assumptions,
        return_result=return_result,
        strategy=strategy,
    )


def resolve_text(
    text: str,
    *,
    symbols=None,
    variable_order=None,
    config=None,
    domain: str | SolveDomain | None = None,
    assumptions=None,
    use_preprocess: bool = True,
    return_result: bool = False,
    strategy: str | None = "lazy",
):
    if variable_order is not None:
        variable_order = [sp.Symbol(v) if isinstance(v, str) else v for v in variable_order]
    parsed = parse_quant_form_text(text, symbols=symbols, variable_order=variable_order)
    preprocess_changed = False
    if use_preprocess:
        prep = semialgebraicize(parsed.matrix)
        preprocess_changed = prep.changed
        if prep.changed:
            parsed = ParsedPrenexFormula(
                tuple(parsed.vars) + tuple(prep.aux_vars),
                tuple(parsed.quantifiers) + tuple(("exists", aux) for aux in prep.aux_vars),
                prep.formula,
                prep.sympy_expr,
            )
    solved = resolve_formula(
        parsed,
        config=config,
        domain=domain,
        assumptions=assumptions,
        return_result=True,
        strategy=strategy,
    )
    solved.preprocess_changed = preprocess_changed
    return solved if return_result else solved.result


__all__ = ["resolve_formula", "resolve_text"]
