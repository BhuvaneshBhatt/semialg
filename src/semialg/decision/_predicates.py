from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence

import sympy as sp
from sympy.logic.boolalg import Boolean
from sympy.polys.polyerrors import CoercionFailed

from ..context import with_computation_context
from ..formula import parse_formula
from ..formulas.boolean import is_false_expr, is_true_expr
from ..incidence import decompose_conjunctive_formula
from ..inequality_reduction import reduce_conjunctive_inequalities
from ..normalization import normalize_formula
from ..qe import qe_by_complete_cad
from ..solve.planner import (
    affine_presolve,
    exact_linear_feasibility,
    profile_semialgebraic_system,
)
from ._inputs import (
    normalize_decision_variables as _normalize_variables,
)
from ._metadata import (
    one_dim_components as _one_dim_components,
)
from ._witnesses import find_validated_witness as _find_validated_witness
from .solution import (
    EquivalenceResult,
    ImplicationResult,
    SatisfiabilityResult,
    TautologyResult,
)

FormulaLike = sp.Expr | Boolean | bool

_RECOVERABLE_ERRORS = (
    ArithmeticError,
    TypeError,
    ValueError,
    NotImplementedError,
    sp.PolynomialError,
    CoercionFailed,
)


from .api import _make_quantified_sentence, _truth_from_qe_result, _try_rur_formula  # noqa: E402


@with_computation_context
def is_satisfiable(
    formula: FormulaLike | Iterable[FormulaLike],
    variables: Sequence[sp.Symbol | str] | None = None,
    *,
    domain: str = "reals",
    strategy: str | None = None,
    return_result: bool = False,
) -> bool | SatisfiabilityResult:
    """Return whether a real semialgebraic formula has a satisfying point.

    By default this returns a Boolean value. With
    ``return_result=True`` it returns a ``SatisfiabilityResult`` containing the
    normalized formula, variable order, backend method, and a validated witness
    when one is cheaply available.
    """

    if domain.lower() not in {"real", "reals", "r", "rr"}:
        raise NotImplementedError("is_satisfiable supports only the real domain")
    expr = normalize_formula(formula)
    vars_ = _normalize_variables(variables, expr)
    method = "trivial"
    witness: Mapping[sp.Symbol, sp.Expr] | None = None
    if is_true_expr(expr):
        sat = True
        witness = {var: sp.Integer(0) for var in vars_}
    elif is_false_expr(expr):
        sat = False
    else:
        components = decompose_conjunctive_formula(expr, vars_)
        if components:
            merged_witness: dict[sp.Symbol, sp.Expr] = {}
            sat = True
            methods: list[str] = []
            for component_formula, component_vars in components:
                component = is_satisfiable(
                    component_formula,
                    component_vars,
                    domain=domain,
                    strategy=strategy,
                    return_result=True,
                )
                methods.append(component.method)
                if not component.satisfiable:
                    sat = False
                    break
                if component.witness:
                    merged_witness.update(component.witness)
            method = "incidence_decomposition[" + ",".join(methods) + "]"
            witness = merged_witness if sat else None
        else:
            sat = None
        if sat is None:
            # Reuse the automatic solver's exact affine feasibility path before
            # finite algebraic solving or CAD.  This is especially important
            # for higher-dimensional polyhedra, where general CAD is avoidable.
            linear_profile = profile_semialgebraic_system(expr, vars_)
            if linear_profile.linear and linear_profile.conjunctive:
                try:
                    linear_presolve = affine_presolve(expr, vars_)
                    linear = exact_linear_feasibility(
                        linear_presolve.formula, linear_presolve.variables
                    )
                except _RECOVERABLE_ERRORS:
                    linear = None
                if linear is not None:
                    sat = bool(linear[0])
                    method = "linear_fourier_motzkin"
                    if sat:
                        witness = _find_validated_witness(expr, vars_, strategy=strategy)
        if sat is None:
            # Exact finite-system dispatch proves supported zero-dimensional
            # equality branches without constructing an ambient CAD.
            rur_result = _try_rur_formula(expr, vars_, max_solutions=1)
            if rur_result is not None and not rur_result.partial:
                sat = bool(rur_result.assignments)
                method = "rational_univariate"
                witness = dict(rur_result.assignments[0]) if rur_result.assignments else None
            else:
                # A validated witness can certify feasibility cheaply; failure to
                # find one carries no information about unsatisfiability.
                witness = _find_validated_witness(expr, vars_, strategy=strategy)
                if witness is not None:
                    sat = True
                    method = "validated_sample"
                else:
                    if len(vars_) == 1:
                        reduced = reduce_conjunctive_inequalities(expr, vars_[0])
                        comps = (
                            _one_dim_components(reduced, vars_[0]) if reduced is not None else None
                        )
                        if comps is not None:
                            sat = bool(comps)
                            method = "sympy_reduce_inequalities"
                            witness = {vars_[0]: comps[0].sample_point()} if comps else None
                        else:
                            result = qe_by_complete_cad(
                                vars_,
                                _make_quantified_sentence(expr, vars_),
                                parse_formula(expr),
                                return_result=True,
                            )
                            sat = _truth_from_qe_result(result)
                            method = getattr(result, "method", "complete_cad_qe")
                    else:
                        result = qe_by_complete_cad(
                            vars_,
                            _make_quantified_sentence(expr, vars_),
                            parse_formula(expr),
                            return_result=True,
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

    expr = normalize_formula(formula)
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

    premise = normalize_formula(assumptions)
    consequent = normalize_formula(conclusion)
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

    left = normalize_formula(lhs)
    right = normalize_formula(rhs)
    universe = _normalize_variables(variables, sp.And(left, right))
    if left == right:
        if return_result:
            return EquivalenceResult(True, left, right, universe, method="syntactic")
        return True
    # Generic Boolean simplification is only safe here for propositional
    # formulas.  Deep simplification of relational atoms can invoke expensive
    # algebraic equality checks before the exact semialgebraic solver runs.
    has_relations = bool(
        left.atoms(sp.core.relational.Relational) or right.atoms(sp.core.relational.Relational)
    )
    if not has_relations:
        try:
            if sp.simplify_logic(sp.Xor(left, right), deep=False) is sp.false:
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
