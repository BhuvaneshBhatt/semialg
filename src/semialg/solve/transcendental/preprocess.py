from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field

import sympy as sp

from .families import classify_trans_fams, default_trans_handlers
from .state import QuantifierBlock, TransProblemState


@dataclass(frozen=True)
class FamilyReplacementStep:
    family_name: str
    original_matches: tuple[sp.Expr, ...]
    auxiliary_variables: tuple[sp.Symbol, ...]
    constraints_added: tuple[sp.Expr, ...]
    changed_formula: bool = True


@dataclass(frozen=True)
class QuantifierDispatchPlan:
    prior_univar_vars: tuple[sp.Symbol, ...]
    quantified_blocks: tuple[QuantifierBlock, ...]
    family_order: tuple[str, ...]
    metadata: dict = field(default_factory=dict)


@dataclass(frozen=True)
class TransPrepResult:
    state: TransProblemState
    steps: tuple[FamilyReplacementStep, ...]
    changed: bool = False
    quantifier_plan: QuantifierDispatchPlan | None = None


def _equation_atoms(formula: sp.Expr) -> tuple[sp.Expr, ...]:
    atoms = list(formula.args) if isinstance(formula, sp.And) else [formula]
    return tuple(a for a in atoms if isinstance(a, sp.Equality))


def simp_piecewise_subexprs(state: TransProblemState) -> TransProblemState:
    f = state.formula
    if not f.has(sp.Piecewise):
        return state
    try:
        f2 = sp.piecewise_fold(f)
        f2 = sp.simplify(f2)
    except Exception:
        f2 = f
    if f2 == f:
        return state
    return state.with_formula(f2, note="piecewise_fold")


def build_quantifier_plan(state: TransProblemState) -> QuantifierDispatchPlan:
    equation_atoms = _equation_atoms(state.formula)
    univariate = []
    for v in state.active_variable_order:
        count = 0
        for eq in equation_atoms:
            expr = sp.simplify(eq.lhs - eq.rhs)
            if expr.free_symbols <= {v}:
                count += 1
        if count:
            univariate.append((v, -count))
    prioritized = tuple(
        v
        for v, _ in sorted(
            univariate, key=lambda t: (t[1], state.active_variable_order.index(t[0]))
        )
    )
    families = classify_trans_fams(state.formula)
    family_order = tuple(
        det.family_name for det in sorted(families, key=lambda d: (-len(d.matches), d.family_name))
    )
    return QuantifierDispatchPlan(
        prior_univar_vars=prioritized,
        quantified_blocks=state.quantifier_blocks,
        family_order=family_order,
        metadata={
            "has_quantifiers": state.has_quantifiers,
            "equation_count": len(equation_atoms),
            "family_counts": {d.family_name: len(d.matches) for d in families},
        },
    )


def _make_aux_symbol(base: str, index: int, all_symbols: set[sp.Symbol]) -> sp.Symbol:
    candidate = sp.Symbol(f"{base}_{index}", real=True)
    while candidate in all_symbols:
        index += 1
        candidate = sp.Symbol(f"{base}_{index}", real=True)
    all_symbols.add(candidate)
    return candidate


def _replacement_constraint(family_name: str, aux: sp.Symbol, expr: sp.Expr) -> sp.Expr:
    if family_name in {
        "trigonometric",
        "hyperbolic",
        "exponential",
        "inverse",
        "productlog",
        "special",
    }:
        return sp.Eq(aux, expr)
    return sp.Eq(aux, expr)


def replace_function_auxilia(
    state: TransProblemState,
    family_name: str,
    *,
    only_in_equations: bool = False,
) -> tuple[TransProblemState, FamilyReplacementStep | None]:
    handlers = {h.family_name: h for h in default_trans_handlers()}
    if family_name not in handlers:
        return state, None
    detector = handlers[family_name].detector
    scope = _equation_atoms(state.formula) if only_in_equations else (state.formula,)
    matches = set()
    for expr in scope:
        matches.update(detector(expr))
    matches = tuple(sorted(matches, key=sp.default_sort_key))
    if not matches:
        return state, None

    used = set(state.all_variables) | {s for m in matches for s in m.free_symbols}
    repl = {}
    aux_vars = []
    constraints = []
    for index, m in enumerate(matches, start=1):
        aux = _make_aux_symbol(f"{family_name}_aux", index, used)
        repl[m] = aux
        aux_vars.append(aux)
        constraints.append(_replacement_constraint(family_name, aux, m))

    new_formula = state.formula.xreplace(repl)
    if constraints:
        new_formula = sp.And(new_formula, *constraints)
    new_qvars = tuple(v for v in state.quantified_variables if v not in repl)
    new_state = TransProblemState(
        formula=sp.simplify(new_formula),
        free_variables=state.free_variables,
        quantified_variables=new_qvars,
        parameter_variables=state.parameter_variables,
        quantifier_blocks=state.quantifier_blocks,
        variable_domains=state.variable_domains,
        default_domain=state.default_domain,
        variable_order=state.variable_order
        + tuple(v for v in aux_vars if v not in state.variable_order),
        notes=state.notes + (f"replaced_{family_name}_family",),
        metadata=dict(state.metadata),
    )
    return new_state, FamilyReplacementStep(
        family_name=family_name,
        original_matches=matches,
        auxiliary_variables=tuple(aux_vars),
        constraints_added=tuple(constraints),
        changed_formula=True,
    )


def prep_trans_problem(
    state: TransProblemState,
    *,
    families_to_replace: Sequence[str] = (
        "trigonometric",
        "hyperbolic",
        "exponential",
        "inverse",
        "productlog",
        "special",
    ),
    quantifier_aware: bool = True,
) -> TransPrepResult:
    current = simp_piecewise_subexprs(state)
    steps = []
    changed = current.formula != state.formula
    plan = build_quantifier_plan(current) if quantifier_aware else None
    family_order = (
        plan.family_order if plan is not None and plan.family_order else tuple(families_to_replace)
    )
    seen = set()
    ordered_families = [
        f for f in family_order if f in families_to_replace and not (f in seen or seen.add(f))
    ]
    ordered_families += [f for f in families_to_replace if f not in ordered_families]

    for family in ordered_families:
        current, step = replace_function_auxilia(
            current,
            family,
            only_in_equations=bool(current.has_quantifiers),
        )
        if step is not None:
            steps.append(step)
            changed = True

    return TransPrepResult(
        state=current,
        steps=tuple(steps),
        changed=changed,
        quantifier_plan=plan,
    )


__all__ = [
    "FamilyReplacementStep",
    "QuantifierDispatchPlan",
    "TransPrepResult",
    "simp_piecewise_subexprs",
    "build_quantifier_plan",
    "replace_function_auxilia",
    "prep_trans_problem",
]
