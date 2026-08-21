from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from itertools import product

import sympy as sp
from sympy import S
from sympy.core.relational import Equality as SymEquality

from ..algebraic.rational_univariate import RationalUnivariateError, solve_formula_with_rur
from ..decomposition.components import component_instances
from ..formula import Formula, ParsedPrenexFormula, parse_formula, parse_quant_form_text, to_sympy
from ..instances.real_fallbacks import find_real_witnesses, satisfies_formula
from ..partial.qe import lazy_find_inst_form
from ..qe.virtual_substitution import try_quadratic_virtual_substitution_witness
from ..status import SolverStatus
from .domains import SolveDomain, apply_assumptions, normalize_assumptions, normalize_domain
from .integer.diophantine import solve_int_methods
from .preprocess import semialgebraicize
from .result import SolveResult


@dataclass(frozen=True)
class InstanceResult:
    """Instances satisfying a formula over a requested domain.

    Exact samples are stored in ``instances``. Numeric approximations are kept
    separately so callers can display approximate points without discarding the
    exact algebraic or rational information returned by CAD.
    """

    instances: tuple[Mapping[sp.Symbol, sp.Expr], ...]
    approximate: tuple[Mapping[sp.Symbol, complex | float], ...]
    variables: tuple[sp.Symbol, ...]
    domain: SolveDomain
    status: str
    method: str
    exact: bool = True
    diagnostics: Mapping[str, object] = field(default_factory=dict)
    random_seed: int | None = None

    @property
    def found(self) -> bool:
        return bool(self.instances)

    def first(self) -> Mapping[sp.Symbol, sp.Expr] | None:
        return self.instances[0] if self.instances else None

    @property
    def result(self) -> Mapping[sp.Symbol, sp.Expr] | None:
        return self.first()

    def __len__(self) -> int:
        return len(self.instances)

    def __iter__(self):
        return iter(self.instances)


def _coerce_vars(variable_order):
    if variable_order is None:
        return None
    return tuple(sp.Symbol(v, real=True) if isinstance(v, str) else v for v in variable_order)


def _norm_vars(variables: Sequence[sp.Symbol | str] | None, expr: sp.Expr) -> tuple[sp.Symbol, ...]:
    if variables is None:
        return tuple(sorted(expr.free_symbols, key=lambda s: s.name))
    out: list[sp.Symbol] = []
    seen: set[sp.Symbol] = set()
    for var in variables:
        sym = sp.Symbol(var, real=True) if isinstance(var, str) else var
        if sym not in seen:
            out.append(sym)
            seen.add(sym)
    return tuple(out)


def _as_approx(inst: Mapping[sp.Symbol, sp.Expr]) -> Mapping[sp.Symbol, complex | float]:
    approx: dict[sp.Symbol, complex | float] = {}
    for var, value in inst.items():
        if value in (sp.true, True):
            approx[var] = 1.0
            continue
        if value in (sp.false, False):
            approx[var] = 0.0
            continue
        num = complex(sp.N(value, 30))
        approx[var] = num.real if abs(num.imag) < 1e-24 else num
    return approx


def _make_result(
    *,
    instances: Sequence[Mapping[sp.Symbol, sp.Expr]],
    variables: Sequence[sp.Symbol],
    domain: SolveDomain,
    method: str,
    diagnostics: Mapping[str, object] | None = None,
    exact: bool = True,
) -> InstanceResult:
    inst_tuple = tuple(dict(inst) for inst in instances)
    status = SolverStatus.SAT if inst_tuple else SolverStatus.UNSAT
    return InstanceResult(
        instances=inst_tuple,
        approximate=tuple(_as_approx(inst) for inst in inst_tuple),
        variables=tuple(variables),
        domain=domain,
        status=status,
        method=method,
        exact=exact,
        diagnostics=dict(diagnostics or {}),
    )


def _finite_set_instances(
    vars_: Sequence[sp.Symbol], solset, count: int
) -> list[dict[sp.Symbol, sp.Expr]]:
    out: list[dict[sp.Symbol, sp.Expr]] = []
    if isinstance(solset, sp.FiniteSet):
        for item in solset:
            values = item if isinstance(item, tuple) else (item,)
            out.append({var: sp.sympify(val) for var, val in zip(vars_, values, strict=True)})
            if len(out) >= count:
                break
    return out


def _real_instances(
    expr: sp.Expr,
    formula: Formula,
    vars_: Sequence[sp.Symbol],
    count: int,
    strategy: str,
    *,
    random_seed: int | None = None,
    strict: bool = True,
) -> InstanceResult:
    """Find real instances using fast fallbacks plus CAD components.

    The fallback pipeline tries cheap algebraic heuristics, bounded sampling,
    random sections, univariate decomposition, and strict candidate validation.
    CAD remains the component-aware backend and is used when requested or when
    the fast pipeline does not find a witness.
    """

    strategy_name = (strategy or "auto").lower()
    diagnostics: dict[str, object] = {}
    fallback_first = strategy_name in {"auto", "fallback", "heuristic", "heuristics", "random"}
    cad_allowed = strategy_name not in {"fallback", "heuristic", "heuristics", "random"}

    if strategy_name in {"auto", "rur", "rational-univariate", "rational_univariate"}:
        try:
            rur = solve_formula_with_rur(expr, vars_, real=True, max_solutions=count)
        except RationalUnivariateError as exc:
            diagnostics["rur_error"] = str(exc)
            rur = None
        if rur is not None:
            diagnostics["rur_result"] = rur
            if rur.assignments or (
                rur.complete
                and strategy_name in {"rur", "rational-univariate", "rational_univariate"}
            ):
                return _make_result(
                    instances=rur.assignments[:count],
                    variables=vars_,
                    domain=SolveDomain.REALS,
                    method="rational_univariate_instance",
                    diagnostics=diagnostics,
                )
            if rur.status == SolverStatus.UNSAT and not cad_allowed:
                return _make_result(
                    instances=[],
                    variables=vars_,
                    domain=SolveDomain.REALS,
                    method="rational_univariate_instance",
                    diagnostics=diagnostics,
                )

    if fallback_first:
        fb = find_real_witnesses(expr, vars_, count=count, seed=random_seed, strict=strict)
        diagnostics["fallback_status"] = fb.status
        diagnostics["fallback_method"] = fb.method
        diagnostics["fallback_attempts"] = tuple(
            {"name": attempt.name, "status": attempt.status, "details": dict(attempt.details)}
            for attempt in fb.attempts
        )
        if fb.instances:
            return _make_result(
                instances=fb.instances[:count],
                variables=vars_,
                domain=SolveDomain.REALS,
                method=f"{fb.method}_instance",
                diagnostics=diagnostics,
                exact=fb.exact,
            )
        if (fb.status == SolverStatus.UNSAT or fb.status == "unsat") and not cad_allowed:
            return _make_result(
                instances=[],
                variables=vars_,
                domain=SolveDomain.REALS,
                method=f"{fb.method}_instance",
                diagnostics=diagnostics,
                exact=fb.exact,
            )

    if cad_allowed:
        comps = component_instances(
            expr, vars_, strategy=strategy, max_components=count, return_result=True
        )
        cad_instances = [dict(inst) for inst in comps.instances[:count]]
        cad_instances = [
            inst for inst in cad_instances if satisfies_formula(expr, inst, strict=strict)
        ]
        diagnostics.update(
            {"component_count": len(comps), "component_diagnostics": comps.diagnostics}
        )
        if cad_instances:
            return _make_result(
                instances=cad_instances[:count],
                variables=vars_,
                domain=SolveDomain.REALS,
                method="component_cad_instance",
                diagnostics=diagnostics,
            )

    method = "component_cad_instance" if cad_allowed else "real_instance_fallback_pipeline_instance"
    return _make_result(
        instances=[],
        variables=vars_,
        domain=SolveDomain.REALS,
        method=method,
        diagnostics=diagnostics,
    )


def _quadratic_virtual_substitution_instance(
    parsed: ParsedPrenexFormula, count: int, strategy: str
) -> InstanceResult | None:
    """Find one real witness by eliminating existential quadratic variables.

    The VS witness backend reduces eligible quantified variables first, asks the
    existing real instance machinery for the reduced formula, and reconstructs
    eliminated variables from the recorded quadratic stages.
    """

    if count != 1 or not parsed.quantifiers:
        return None
    if any(str(quantifier).lower() != "exists" for quantifier, _ in parsed.quantifiers):
        return None

    def base_finder(reduced_formula: sp.Expr, reduced_variables: Sequence[sp.Symbol]):
        reduced_variables = tuple(reduced_variables)
        if not reduced_variables:
            try:
                truth = sp.simplify(reduced_formula)
            except (ValueError, TypeError, sp.SympifyError):
                truth = reduced_formula
            return {} if truth == sp.true or truth is sp.true else None

        candidate_lists: list[list[sp.Expr]] = []
        atoms = [atom for atom in reduced_formula.atoms(SymEquality)]
        for variable in reduced_variables:
            values: list[sp.Expr] = [sp.Integer(0), sp.Integer(1), sp.Integer(-1)]
            for atom in atoms:
                if isinstance(atom, sp.Equality):
                    expr = sp.expand(atom.lhs - atom.rhs)
                    if variable in expr.free_symbols and not (expr.free_symbols - {variable}):
                        try:
                            values.extend(sp.solve(sp.Eq(expr, 0), variable))
                        except (ValueError, TypeError, sp.SympifyError):
                            pass
            deduped: list[sp.Expr] = []
            seen: set[str] = set()
            for value in values:
                value = sp.simplify(value)
                key = sp.sstr(value)
                if key not in seen:
                    seen.add(key)
                    deduped.append(value)
            candidate_lists.append(deduped)
        for values in product(*candidate_lists):
            candidate = dict(zip(reduced_variables, values, strict=True))
            evaluated = reduced_formula.subs(candidate)
            try:
                evaluated = sp.simplify(evaluated)
            except (ValueError, TypeError, sp.SympifyError):
                pass
            if evaluated == sp.true or evaluated is sp.true:
                return candidate

        base = _real_instances(
            reduced_formula,
            parse_formula(reduced_formula),
            reduced_variables,
            1,
            strategy if strategy not in {"lazy", "partial"} else "auto",
            strict=True,
        )
        return base.first()

    witness = try_quadratic_virtual_substitution_witness(
        parsed.vars,
        parsed.quantifiers,
        parsed.matrix_expr,
        base_finder,
        full=True,
    )
    if witness is None:
        return None
    diagnostics = {
        "vs_witness_result": witness,
        "reduced_formula": witness.reduced_formula,
        "eliminated_variables": witness.eliminated_variables,
        "notes": witness.notes,
    }
    if witness.instance is None:
        return _make_result(
            instances=[],
            variables=parsed.vars,
            domain=SolveDomain.REALS,
            method="quadratic_virtual_substitution_instance",
            diagnostics=diagnostics,
        )
    return _make_result(
        instances=[witness.instance],
        variables=parsed.vars,
        domain=SolveDomain.REALS,
        method="quadratic_virtual_substitution_instance",
        diagnostics=diagnostics,
    )


def _lazy_real_instance(parsed: ParsedPrenexFormula, count: int) -> InstanceResult:
    lazy = lazy_find_inst_form(parsed.vars, parsed.matrix, quantifiers=parsed.quantifiers)
    instances = [lazy.instance] if lazy.instance is not None else []
    return _make_result(
        instances=instances[:count],
        variables=parsed.vars,
        domain=SolveDomain.REALS,
        method="partial_cad_instance",
        diagnostics={"lazy_result": lazy, "stats": lazy.stats, "found": lazy.found},
    )


def _complex_instances(expr: sp.Expr, vars_: Sequence[sp.Symbol], count: int) -> InstanceResult:
    instances: list[dict[sp.Symbol, sp.Expr]] = []
    if len(vars_) == 1:
        sol = sp.solveset(expr, vars_[0], domain=S.Complexes)
        instances = _finite_set_instances(vars_, sol, count)
    else:
        try:
            equations = [
                arg.lhs - arg.rhs if isinstance(arg, sp.Equality) else arg
                for arg in sp.And.make_args(expr)
            ]
            sol = sp.nonlinsolve(equations, vars_)
            instances = _finite_set_instances(vars_, sol, count)
        except (ValueError, TypeError, sp.SympifyError) as exc:
            return _make_result(
                instances=[],
                variables=vars_,
                domain=SolveDomain.COMPLEXES,
                method="sympy_complex_instance",
                diagnostics={"status_reason": f"complex solving failed: {exc}"},
                exact=False,
            )
    return _make_result(
        instances=instances,
        variables=vars_,
        domain=SolveDomain.COMPLEXES,
        method="sympy_complex_instance",
    )


def _integer_instances(expr: sp.Expr, vars_: Sequence[sp.Symbol], count: int) -> InstanceResult:
    instances: list[dict[sp.Symbol, sp.Expr]] = []
    specialized = solve_int_methods(expr, vars_)
    if specialized is not None and specialized.solutions:
        instances = [
            {var: sp.sympify(val) for var, val in zip(vars_, pt, strict=True)}
            for pt in specialized.solutions[:count]
        ]
    elif len(vars_) == 1:
        sol = sp.solveset(expr, vars_[0], domain=S.Integers)
        instances = _finite_set_instances(vars_, sol, count)
    return _make_result(
        instances=instances, variables=vars_, domain=SolveDomain.INTEGERS, method="integer_instance"
    )


def _boolean_instances(expr: sp.Expr, vars_: Sequence[sp.Symbol], count: int) -> InstanceResult:
    instances: list[dict[sp.Symbol, sp.Expr]] = []
    for values in product([False, True], repeat=len(vars_)):
        subs = dict(zip(vars_, values, strict=True))
        try:
            ok = bool(expr.subs(subs))
        except TypeError:
            ok = bool(sp.simplify(expr.subs(subs)))
        if ok:
            instances.append({var: sp.true if val else sp.false for var, val in subs.items()})
            if len(instances) >= count:
                break
    return _make_result(
        instances=instances,
        variables=vars_,
        domain=SolveDomain.BOOLEANS,
        method="boolean_bruteforce_instance",
        exact=True,
    )


def find_instance(
    formula: sp.Expr | Formula,
    variables: Sequence[sp.Symbol | str] | None = None,
    *,
    domain: str | SolveDomain | None = None,
    count: int = 1,
    strategy: str = "auto",
    assumptions: Iterable[sp.Expr] | sp.Expr | None = None,
    random_seed: int | None = None,
    exact: bool = True,
    strict: bool = False,
    return_result: bool = True,
):
    """Find satisfying assignments for a formula.

    For real quantifier-free formulas, CAD component extraction supplies one
    representative per connected component before any fallback sampling is used.
    Non-real domains currently use exact SymPy/specialized finite methods and
    return ``unknown``-style diagnostics through the result object when they
    cannot solve a requested class.
    """

    base_expr = (
        to_sympy(formula)
        if not isinstance(formula, (sp.Basic, sp.logic.boolalg.Boolean))
        else formula
    )
    expr = apply_assumptions(base_expr, assumptions)
    vars_ = _norm_vars(variables, expr)
    dom = normalize_domain(domain)
    if count <= 0:
        result = _make_result(
            instances=[], variables=vars_, domain=dom, method="empty_instance_request"
        )
    elif dom is SolveDomain.REALS:
        result = _real_instances(
            expr,
            parse_formula(expr),
            vars_,
            count,
            strategy,
            random_seed=random_seed,
            strict=strict,
        )
    elif dom is SolveDomain.COMPLEXES:
        result = _complex_instances(expr, vars_, count)
    elif dom is SolveDomain.INTEGERS:
        result = _integer_instances(expr, vars_, count)
    elif dom is SolveDomain.BOOLEANS:
        result = _boolean_instances(expr, vars_, count)
    else:
        result = InstanceResult(
            (),
            (),
            vars_,
            dom,
            "unknown",
            "unsupported_domain",
            diagnostics={"reason": f"unsupported domain {dom}"},
        )
    if random_seed is not None or assumptions is not None or not exact:
        diag = dict(result.diagnostics)
        diag.update(
            {
                "random_seed": random_seed,
                "assumptions": tuple(map(sp.sstr, normalize_assumptions(assumptions))),
                "exact_requested": exact,
                "strict": strict,
            }
        )
        result = InstanceResult(
            result.instances,
            result.approximate,
            result.variables,
            result.domain,
            result.status,
            result.method,
            exact=result.exact and exact,
            diagnostics=diag,
            random_seed=random_seed,
        )
    if return_result:
        return result
    if count == 1:
        return result.first()
    return result.instances


def find_instance_formula(
    parsed: ParsedPrenexFormula,
    config=None,
    *,
    domain: str | SolveDomain | None = None,
    max_instances: int | None = None,
    count: int | None = None,
    return_result: bool = True,
    strategy: str | None = "lazy",
    assumptions: Iterable[sp.Expr] | sp.Expr | None = None,
    random_seed: int | None = None,
    exact: bool = True,
    strict: bool = False,
):
    """Find a witness for a structured formula using exact domain-aware search strategies."""
    inst_count = count if count is not None else (max_instances if max_instances is not None else 1)
    dom = normalize_domain(domain)
    expr = apply_assumptions(to_sympy(parsed.matrix), assumptions)
    if dom is SolveDomain.REALS:
        strategy_name = (strategy or "lazy").lower()
        rur_result = None
        use_rur = all(str(q).lower() == "exists" for q, _ in parsed.quantifiers) and (
            bool(parsed.quantifiers)
            or strategy_name in {"rur", "rational-univariate", "rational_univariate"}
        )
        if use_rur:
            try:
                rur_result = solve_formula_with_rur(
                    expr, parsed.vars, real=True, max_solutions=inst_count
                )
            except RationalUnivariateError:
                rur_result = None
        if rur_result is not None and (
            rur_result.assignments
            or (
                rur_result.complete
                and strategy_name in {"rur", "rational-univariate", "rational_univariate"}
            )
        ):
            result = _make_result(
                instances=rur_result.assignments[:inst_count],
                variables=parsed.vars,
                domain=SolveDomain.REALS,
                method="rational_univariate_instance",
                diagnostics={"rur_result": rur_result},
            )
        else:
            vs_result = _quadratic_virtual_substitution_instance(parsed, inst_count, strategy_name)
            if vs_result is not None and (
                vs_result.found or vs_result.status == SolverStatus.UNSAT
            ):
                result = vs_result
            elif parsed.quantifiers or strategy_name in {"lazy", "partial"}:
                result = _lazy_real_instance(parsed, inst_count)
            else:
                result = _real_instances(
                    expr,
                    parse_formula(expr),
                    parsed.vars,
                    inst_count,
                    strategy_name,
                    random_seed=random_seed,
                    strict=strict,
                )
    elif parsed.quantifiers:
        result = InstanceResult(
            (),
            (),
            parsed.vars,
            dom,
            "unknown",
            "quantified_nonreal_instance",
            diagnostics={"reason": "quantified non-real instance finding is not implemented"},
        )
    elif dom is SolveDomain.COMPLEXES:
        result = _complex_instances(expr, parsed.vars, inst_count)
    elif dom is SolveDomain.INTEGERS:
        result = _integer_instances(expr, parsed.vars, inst_count)
    elif dom is SolveDomain.BOOLEANS:
        result = _boolean_instances(expr, parsed.vars, inst_count)
    else:
        result = InstanceResult((), (), parsed.vars, dom, "unknown", "unsupported_domain")
    if random_seed is not None or assumptions is not None or not exact:
        diag = dict(result.diagnostics)
        diag.update(
            {
                "random_seed": random_seed,
                "assumptions": tuple(map(sp.sstr, normalize_assumptions(assumptions))),
                "exact_requested": exact,
                "strict": strict,
            }
        )
        result = InstanceResult(
            result.instances,
            result.approximate,
            result.variables,
            result.domain,
            result.status,
            result.method,
            exact=result.exact and exact,
            diagnostics=diag,
            random_seed=random_seed,
        )
    if return_result:
        return result
    return result.first() if inst_count == 1 else result.instances


def find_instance_text(
    text: str,
    *,
    symbols=None,
    variable_order=None,
    variables=None,
    config=None,
    domain: str | SolveDomain | None = None,
    max_instances: int | None = None,
    count: int | None = None,
    use_preprocess: bool = True,
    return_result: bool = True,
    strategy: str | None = "lazy",
    assumptions: Iterable[sp.Expr] | sp.Expr | None = None,
    random_seed: int | None = None,
    exact: bool = True,
    strict: bool = False,
):
    """Find a witness for a textual formula after parsing and domain normalization."""
    variable_order = _coerce_vars(variable_order or variables)
    parse_symbols = dict(symbols or {})
    if variable_order is not None:
        for var in variable_order:
            parse_symbols.setdefault(var.name, var)
    parsed = parse_quant_form_text(text, symbols=parse_symbols, variable_order=variable_order)
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
    result = find_instance_formula(
        parsed,
        config=config,
        domain=domain,
        count=count,
        max_instances=max_instances,
        return_result=True,
        strategy=strategy,
        assumptions=assumptions,
        random_seed=random_seed,
        exact=exact,
        strict=strict,
    )
    diagnostics = dict(result.diagnostics)
    diagnostics["preprocess_changed"] = preprocess_changed
    result = InstanceResult(
        result.instances,
        result.approximate,
        result.variables,
        result.domain,
        result.status,
        result.method,
        exact=result.exact,
        diagnostics=diagnostics,
        random_seed=random_seed,
    )
    if return_result:
        return result
    inst_count = count if count is not None else (max_instances if max_instances is not None else 1)
    return result.first() if inst_count == 1 else result.instances


def instances_solve_result(result: InstanceResult) -> SolveResult:
    value = result.first() if len(result) == 1 else result.instances
    return SolveResult(
        method=result.method,
        domain=result.domain,
        result=value,
        metadata={"instance_result": result},
    )


__all__ = [
    "InstanceResult",
    "find_instance",
    "find_instance_formula",
    "find_instance_text",
    "instances_solve_result",
]
