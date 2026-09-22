from __future__ import annotations

from collections.abc import Iterable, Sequence

import sympy as sp
from sympy.logic.boolalg import Boolean

from ._errors import EXACT_OPERATION_ERRORS as _RECOVERABLE_ERRORS
from .algebraic.rational_univariate import RationalUnivariateError
from .conditional import ParameterStratifiedResult, conditional_result
from .decision import implies, is_satisfiable
from .domain_solve import function_domain
from .normalization import normalize_formula, normalize_variables
from .parameters import solvability_conditions
from .reasoning_results import SignClassificationResult, SignProofResult

FormulaLike = sp.Expr | Boolean | bool


def _relation_formula(expr: sp.Expr, relation: str) -> sp.Expr:
    expression = sp.sympify(expr)
    if relation == "positive":
        return expression > 0
    if relation == "nonnegative":
        return expression >= 0
    if relation == "negative":
        return expression < 0
    if relation == "nonpositive":
        return expression <= 0
    if relation == "zero":
        return sp.Eq(expression, 0)
    if relation == "nonzero":
        return sp.Ne(expression, 0)
    raise ValueError(f"unknown sign proof relation: {relation!r}")


def _constant_sign_certificate(expr: sp.Expr, relation: str) -> bool | None:
    try:
        value = sp.simplify(expr)
    except _RECOVERABLE_ERRORS:
        value = expr
    if getattr(value, "free_symbols", set()):
        return None
    try:
        if relation == "positive":
            return bool(value > 0)
        if relation == "nonnegative":
            return bool(value >= 0)
        if relation == "negative":
            return bool(value < 0)
        if relation == "nonpositive":
            return bool(value <= 0)
        if relation == "zero":
            return bool(value == 0)
        if relation == "nonzero":
            return bool(value != 0)
    except _RECOVERABLE_ERRORS:
        return None
    return None


def _is_even_power_nonnegative(expr: sp.Expr) -> bool:
    if expr == 0:
        return True
    _base, exp = expr.as_base_exp()
    return bool(exp.is_integer and exp.is_even and exp.is_nonnegative)


def _is_obvious_square_product_nonnegative(expr: sp.Expr, variables: Sequence[sp.Symbol]) -> bool:
    try:
        factored = sp.factor(expr)
    except _RECOVERABLE_ERRORS:
        factored = expr
    if factored == 0:
        return True
    if _is_even_power_nonnegative(factored):
        return True
    coeff, factors = factored.as_coeff_mul()
    try:
        if not bool(coeff >= 0):
            return False
    except _RECOVERABLE_ERRORS:
        return False
    for factor in factors:
        if _is_even_power_nonnegative(factor):
            continue
        if getattr(factor, "free_symbols", set()):
            return False
        try:
            if not bool(factor >= 0):
                return False
        except _RECOVERABLE_ERRORS:
            return False
    return True


def _is_sum_of_squares_nonnegative(expr: sp.Expr, variables: Sequence[sp.Symbol]) -> bool:
    try:
        expanded = sp.expand(expr)
    except _RECOVERABLE_ERRORS:
        expanded = expr
    terms = sp.Add.make_args(expanded) if isinstance(expanded, sp.Add) else (expanded,)
    return bool(terms) and all(
        _is_obvious_square_product_nonnegative(term, variables) for term in terms
    )


def _zero_vector_counterexample(
    expr: sp.Expr, relation: str, variables: Sequence[sp.Symbol]
) -> dict[sp.Symbol, sp.Expr] | None:
    if relation not in {"positive", "negative"}:
        return None
    point = {var: sp.Integer(0) for var in variables}
    try:
        value = sp.simplify(sp.sympify(expr).subs(point))
    except _RECOVERABLE_ERRORS:
        return None
    if value == 0:
        return point
    return None


def _cheap_sign_proof(
    expr: sp.Expr, relation: str, assumptions: sp.Expr, variables: Sequence[sp.Symbol]
) -> tuple[bool, str, object] | None:
    if assumptions not in (sp.true, True):
        return None
    zero_counterexample = _zero_vector_counterexample(expr, relation, variables)
    if zero_counterexample is not None:
        return (False, "zero_counterexample", zero_counterexample)
    const = _constant_sign_certificate(expr, relation)
    if const is not None:
        return (const, "constant_sign", sp.simplify(expr))
    expression = sp.sympify(expr)
    if relation == "nonnegative" and _is_sum_of_squares_nonnegative(expression, variables):
        return (True, "sum_of_squares_or_even_powers", "syntactic_sum_of_squares")
    if relation == "nonpositive" and _is_sum_of_squares_nonnegative(-expression, variables):
        return (True, "negative_sum_of_squares_or_even_powers", "syntactic_sum_of_squares")
    target = (
        expression if relation == "positive" else -expression if relation == "negative" else None
    )
    if target is not None:
        try:
            const_term = (
                sp.Poly(target, *variables).coeff_monomial(1) if variables else sp.sympify(target)
            )
            remainder = sp.expand(target - const_term)
            if const_term.is_positive and _is_sum_of_squares_nonnegative(remainder, variables):
                return (
                    True,
                    "positive_constant_plus_squares",
                    {"constant": const_term, "remainder": remainder},
                )
        except _RECOVERABLE_ERRORS:
            pass
    return None


def _prove_sign(
    expr: sp.Expr,
    relation: str,
    variables: Sequence[sp.Symbol | str] | None = None,
    *,
    assumptions: FormulaLike | Iterable[FormulaLike] = True,
    strategy: str | None = None,
    parameters: Sequence[sp.Symbol | str] | None = None,
    return_result: bool = False,
) -> bool | SignProofResult | ParameterStratifiedResult:
    """Prove a sign relation by cheap structure checks, then exact implication.

    Structural certificates are preferred because they avoid CAD. If no cheap
    proof applies, the requested relation is proved as an implication from the
    assumptions; result mode retains the counterexample or implication proof.
    """
    expression = sp.sympify(expr)
    formula = _relation_formula(expression, relation)
    asm = normalize_formula(assumptions)
    vars_ = normalize_variables(
        variables,
        sp.And(asm, formula),
        append_context_symbols=variables is None and parameters is None,
    )
    if parameters is not None:
        normalized_parameters = normalize_variables(
            parameters, sp.And(asm, formula), append_context_symbols=False
        )
        if set(vars_) & set(normalized_parameters):
            raise ValueError("variables and parameters must be disjoint")
        violation_formula = sp.And(asm, sp.Not(formula))
        violation_condition = solvability_conditions(
            violation_formula, vars_, normalized_parameters
        )
        proof_condition = sp.Not(violation_condition)
        return conditional_result(
            normalized_parameters,
            ((proof_condition, True), (sp.Not(proof_condition), False)),
            method=f"{relation}_on_parameter_qe",
            diagnostics={"expression": sp.sstr(expression), "relation": relation},
        )
    cheap = _cheap_sign_proof(expression, relation, asm, vars_)
    if cheap is None and asm in (sp.true, True) and relation in {"nonnegative", "nonpositive"}:
        target = expression if relation == "nonnegative" else -expression
        try:
            from .polynomial_positivity import polynomial_nonnegative

            portfolio = polynomial_nonnegative(target, vars_, return_result=True)
        except (
            RationalUnivariateError,
            sp.PolynomialError,
            ValueError,
            TypeError,
            NotImplementedError,
        ):
            portfolio = None
        if portfolio is not None:
            cheap = (
                portfolio.decision,
                f"polynomial_portfolio:{portfolio.backend}",
                portfolio.certificate if portfolio.decision else portfolio.witness,
            )
    if cheap is not None:
        proven, method, certificate = cheap
        counterexample = None
        if return_result and not proven:
            implication = implies(asm, formula, vars_, strategy=strategy, return_result=True)
            counterexample = implication.counterexample
            if counterexample is not None:
                method = f"{method}+counterexample"
        if return_result:
            if (
                not proven
                and isinstance(certificate, dict)
                and all(var in certificate for var in vars_)
            ):
                counterexample = certificate
            return SignProofResult(
                bool(proven),
                relation,
                expression,
                asm,
                vars_,
                formula,
                counterexample=dict(counterexample) if counterexample else None,
                method=method,
                certificate=certificate,
                diagnostics={"strategy": strategy},
            )
        return bool(proven)
    implication = implies(asm, formula, vars_, strategy=strategy, return_result=True)
    proven = bool(implication)
    if return_result:
        return SignProofResult(
            proven,
            relation,
            expression,
            asm,
            vars_,
            formula,
            counterexample=dict(implication.counterexample) if implication.counterexample else None,
            method=implication.method,
            certificate=implication,
            diagnostics={"strategy": strategy, "implication_diagnostics": implication.diagnostics},
        )
    return proven


def prove_positive(
    expr: sp.Expr,
    variables: Sequence[sp.Symbol | str] | None = None,
    *,
    assumptions: FormulaLike | Iterable[FormulaLike] = True,
    strategy: str | None = None,
    parameters: Sequence[sp.Symbol | str] | None = None,
    return_result: bool = False,
) -> bool | SignProofResult | ParameterStratifiedResult:
    """Return whether the expression is certified positive on the stated domain."""
    return _prove_sign(
        expr,
        "positive",
        variables,
        assumptions=assumptions,
        strategy=strategy,
        parameters=parameters,
        return_result=return_result,
    )


def prove_nonnegative(
    expr: sp.Expr,
    variables: Sequence[sp.Symbol | str] | None = None,
    *,
    assumptions: FormulaLike | Iterable[FormulaLike] = True,
    strategy: str | None = None,
    parameters: Sequence[sp.Symbol | str] | None = None,
    return_result: bool = False,
) -> bool | SignProofResult | ParameterStratifiedResult:
    """Return whether the expression is certified nonnegative on the stated domain."""
    return _prove_sign(
        expr,
        "nonnegative",
        variables,
        assumptions=assumptions,
        strategy=strategy,
        parameters=parameters,
        return_result=return_result,
    )


def prove_negative(
    expr: sp.Expr,
    variables: Sequence[sp.Symbol | str] | None = None,
    *,
    assumptions: FormulaLike | Iterable[FormulaLike] = True,
    strategy: str | None = None,
    parameters: Sequence[sp.Symbol | str] | None = None,
    return_result: bool = False,
) -> bool | SignProofResult | ParameterStratifiedResult:
    """Return whether the expression is certified negative on the stated domain."""
    return _prove_sign(
        expr,
        "negative",
        variables,
        assumptions=assumptions,
        strategy=strategy,
        parameters=parameters,
        return_result=return_result,
    )


def prove_nonpositive(
    expr: sp.Expr,
    variables: Sequence[sp.Symbol | str] | None = None,
    *,
    assumptions: FormulaLike | Iterable[FormulaLike] = True,
    strategy: str | None = None,
    parameters: Sequence[sp.Symbol | str] | None = None,
    return_result: bool = False,
) -> bool | SignProofResult | ParameterStratifiedResult:
    """Return whether the expression is certified nonpositive on the stated domain."""
    return _prove_sign(
        expr,
        "nonpositive",
        variables,
        assumptions=assumptions,
        strategy=strategy,
        parameters=parameters,
        return_result=return_result,
    )


def prove_zero(
    expr: sp.Expr,
    variables: Sequence[sp.Symbol | str] | None = None,
    *,
    assumptions: FormulaLike | Iterable[FormulaLike] = True,
    strategy: str | None = None,
    parameters: Sequence[sp.Symbol | str] | None = None,
    return_result: bool = False,
) -> bool | SignProofResult | ParameterStratifiedResult:
    """Return whether the expression is identically zero on the stated domain."""
    return _prove_sign(
        expr,
        "zero",
        variables,
        assumptions=assumptions,
        strategy=strategy,
        parameters=parameters,
        return_result=return_result,
    )


def prove_nonzero(
    expr: sp.Expr,
    variables: Sequence[sp.Symbol | str] | None = None,
    *,
    assumptions: FormulaLike | Iterable[FormulaLike] = True,
    strategy: str | None = None,
    parameters: Sequence[sp.Symbol | str] | None = None,
    return_result: bool = False,
) -> bool | SignProofResult | ParameterStratifiedResult:
    """Return whether the expression is everywhere nonzero on the stated domain."""
    return _prove_sign(
        expr,
        "nonzero",
        variables,
        assumptions=assumptions,
        strategy=strategy,
        parameters=parameters,
        return_result=return_result,
    )


def function_sign(
    expr: sp.Expr,
    variables: Sequence[sp.Symbol | str] | None = None,
    *,
    assumptions: FormulaLike | Iterable[FormulaLike] = True,
    strategy: str | None = None,
    parameters: Sequence[sp.Symbol | str] | None = None,
    return_result: bool = False,
):
    """Classify the exact sign behavior of ``expr`` on a semialgebraic region.

    Canonical classifications are ``empty_domain``, ``positive``, ``negative``,
    ``zero``, ``nonnegative``, ``nonpositive``, ``nonzero``, and ``mixed``.
    The supplied assumptions are intersected automatically with the exact real
    domain recognized by :func:`function_domain`.  With parameters, a
    stratified classification over parameter space is returned.
    """
    expression = sp.sympify(expr)
    explicit_domain = normalize_formula(assumptions)
    combined = sp.And(explicit_domain, sp.Eq(expression, expression))
    normalized_variables = normalize_variables(
        variables,
        combined,
        append_context_symbols=variables is None and parameters is None,
    )
    normalized_parameters = (
        normalize_variables(parameters, combined, append_context_symbols=False)
        if parameters is not None
        else ()
    )
    natural_domain = function_domain(expression, (*normalized_variables, *normalized_parameters))
    domain = normalize_formula(sp.And(explicit_domain, natural_domain))
    if parameters is not None:
        conditions = {}
        relation_formulas = {
            "positive": expression > 0,
            "negative": expression < 0,
            "zero": sp.Eq(expression, 0),
            "nonnegative": expression >= 0,
            "nonpositive": expression <= 0,
            "nonzero": sp.Ne(expression, 0),
        }
        for name, relation_formula in relation_formulas.items():
            violating_parameters = solvability_conditions(
                sp.And(domain, sp.Not(relation_formula)),
                normalized_variables,
                normalized_parameters,
            )
            conditions[name] = sp.Not(violating_parameters)
        feasible_parameters = solvability_conditions(
            domain, normalized_variables, normalized_parameters
        )
        branches = [
            (sp.Not(feasible_parameters), "empty_domain"),
            (sp.And(feasible_parameters, conditions["positive"]), "positive"),
            (sp.And(feasible_parameters, conditions["negative"]), "negative"),
            (sp.And(feasible_parameters, conditions["zero"]), "zero"),
            (
                sp.And(
                    feasible_parameters,
                    conditions["nonnegative"],
                    sp.Not(conditions["positive"]),
                    sp.Not(conditions["zero"]),
                ),
                "nonnegative",
            ),
            (
                sp.And(
                    feasible_parameters,
                    conditions["nonpositive"],
                    sp.Not(conditions["negative"]),
                    sp.Not(conditions["zero"]),
                ),
                "nonpositive",
            ),
            (
                sp.And(
                    feasible_parameters,
                    conditions["nonzero"],
                    sp.Not(conditions["positive"]),
                    sp.Not(conditions["negative"]),
                ),
                "nonzero",
            ),
        ]
        covered = sp.Or(*(condition for condition, _ in branches))
        branches.append((sp.And(feasible_parameters, sp.Not(covered)), "mixed"))
        return conditional_result(
            normalized_parameters, branches, method="sign_classification_parameter_qe"
        )
    if not is_satisfiable(domain, normalized_variables):
        classification = "empty_domain"
        proofs = {}
    else:
        proofs = {}
        classification = "mixed"
        for relation in ("positive", "negative", "zero", "nonnegative", "nonpositive", "nonzero"):
            proof = _prove_sign(
                expression,
                relation,
                normalized_variables,
                assumptions=domain,
                strategy=strategy,
                return_result=True,
            )
            proofs[relation] = proof
            if proof.proven:
                classification = relation
                break
    result = SignClassificationResult(
        classification, expression, domain, tuple(normalized_variables), proofs
    )
    return result if return_result else classification


__all__ = [
    "SignClassificationResult",
    "SignProofResult",
    "prove_positive",
    "prove_nonnegative",
    "prove_negative",
    "prove_nonpositive",
    "prove_zero",
    "prove_nonzero",
    "function_sign",
]
