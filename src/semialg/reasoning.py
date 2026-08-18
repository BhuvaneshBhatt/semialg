from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field

import sympy as sp
from sympy.logic.boolalg import Boolean

from .decision import equivalent, implies, is_satisfiable
from .formula import parse_formula
from .qe import qe_by_complete_cad
from .regions.operations import region_closure
from .symbolic_simplify import simplify_boole, simplify_piecewise

FormulaLike = sp.Expr | Boolean | bool


@dataclass(frozen=True)
class SimplifiedSystem:
    """Result of CAD/QE-backed constraint-system simplification."""

    formula: sp.Expr
    constraints: tuple[sp.Expr, ...]
    variables: tuple[sp.Symbol, ...]
    inconsistent: bool
    removed_redundant: tuple[sp.Expr, ...] = ()
    substitutions: dict[sp.Symbol, sp.Expr] = field(default_factory=dict)
    method: str = "cad_implication_redundancy"
    diagnostics: dict[str, object] = field(default_factory=dict)

    def __bool__(self) -> bool:
        return not self.inconsistent


@dataclass(frozen=True)
class AssumptionSimplificationResult:
    """Structured result for assumption-aware expression simplification."""

    expression: sp.Expr
    original: sp.Expr
    assumptions: sp.Expr
    variables: tuple[sp.Symbol, ...]
    conditions: tuple[sp.Expr, ...] = ()
    rewrites: tuple[str, ...] = ()
    diagnostics: dict[str, object] = field(default_factory=dict)

    def __bool__(self) -> bool:
        return self.expression is not sp.nan


@dataclass(frozen=True)
class SignProofResult:
    """Structured result for proving the sign of a real expression."""

    proven: bool
    relation: str
    expression: sp.Expr
    assumptions: sp.Expr
    variables: tuple[sp.Symbol, ...]
    formula: sp.Expr
    counterexample: dict[sp.Symbol, sp.Expr] | None = None
    method: str = "unknown"
    certificate: object | None = None
    diagnostics: dict[str, object] = field(default_factory=dict)

    def __bool__(self) -> bool:
        return self.proven


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
    return formula if isinstance(formula, (sp.Basic, Boolean)) else sp.sympify(formula)


def _normalize_variables(
    variables: Sequence[sp.Symbol | str] | None,
    expr: sp.Expr,
) -> tuple[sp.Symbol, ...]:
    out: list[sp.Symbol] = []
    seen: set[sp.Symbol] = set()
    if variables is not None:
        for var in variables:
            sym = _as_real_symbol(var)
            if sym not in seen:
                out.append(sym)
                seen.add(sym)
    for sym in sorted(expr.free_symbols, key=lambda item: item.name):
        if sym not in seen:
            out.append(sym)
            seen.add(sym)
    return tuple(out)


def _conjuncts(expr: sp.Expr) -> tuple[sp.Expr, ...]:
    if expr is sp.true or expr == sp.true:
        return ()
    if isinstance(expr, sp.And):
        out: list[sp.Expr] = []
        for arg in expr.args:
            out.extend(_conjuncts(arg))
        return tuple(out)
    return (expr,)


def _ordered_unique_exprs(exprs: Iterable[sp.Expr]) -> tuple[sp.Expr, ...]:
    out: list[sp.Expr] = []
    seen: set[str] = set()
    for expr in exprs:
        simplified = sp.simplify(expr) if not getattr(expr, "is_Relational", False) else expr
        key = sp.sstr(simplified)
        if key not in seen:
            out.append(simplified)
            seen.add(key)
    return tuple(out)


def _try_reduce_inequalities_1d(expr: sp.Expr, variables: Sequence[sp.Symbol]) -> sp.Expr | None:
    if len(variables) != 1:
        return None
    var = variables[0]
    try:
        reduced = sp.reduce_inequalities(list(_conjuncts(expr)), var)
    except Exception:
        return None
    if reduced is None:
        return None
    return reduced


def _simple_equality_substitutions(
    atoms: Sequence[sp.Expr], variables: Sequence[sp.Symbol]
) -> dict[sp.Symbol, sp.Expr]:
    substitutions: dict[sp.Symbol, sp.Expr] = {}
    variable_set = set(variables)
    for atom in atoms:
        if not isinstance(atom, sp.Equality):
            continue
        pairs = ((atom.lhs, atom.rhs), (atom.rhs, atom.lhs))
        for lhs, rhs in pairs:
            if (
                isinstance(lhs, sp.Symbol)
                and lhs in variable_set
                and lhs not in getattr(rhs, "free_symbols", set())
            ):
                if lhs not in substitutions:
                    substitutions[lhs] = sp.simplify(rhs)
                break
    # Keep only acyclic substitutions in user-variable order.
    clean: dict[sp.Symbol, sp.Expr] = {}
    for sym in variables:
        if sym not in substitutions:
            continue
        rhs = substitutions[sym].xreplace(clean)
        if sym not in getattr(rhs, "free_symbols", set()):
            clean[sym] = rhs
    return clean


def _apply_substitutions_to_constraints(
    atoms: Sequence[sp.Expr], substitutions: dict[sp.Symbol, sp.Expr], *, keep_definitions: bool
) -> tuple[sp.Expr, ...]:
    out: list[sp.Expr] = []
    for atom in atoms:
        if isinstance(atom, sp.Equality):
            if atom.lhs in substitutions and sp.simplify(atom.rhs - substitutions[atom.lhs]) == 0:
                if keep_definitions:
                    out.append(atom)
                continue
            if atom.rhs in substitutions and sp.simplify(atom.lhs - substitutions[atom.rhs]) == 0:
                if keep_definitions:
                    out.append(atom)
                continue
        try:
            new_atom = atom.subs(substitutions)
            if new_atom in (True, sp.true):
                continue
            if new_atom in (False, sp.false):
                out.append(sp.false)
            else:
                out.append(new_atom)
        except Exception:
            out.append(atom)
    return tuple(out)


def simplify_system(
    constraints: FormulaLike | Iterable[FormulaLike],
    variables: Sequence[sp.Symbol | str] | None = None,
    *,
    assumptions: FormulaLike | Iterable[FormulaLike] = True,
    return_result: bool = False,
    strategy: str | None = None,
    eliminate_equalities: bool = False,
    output: str = "formula",
) -> sp.Expr | tuple[sp.Expr, ...] | SimplifiedSystem:
    """Simplify a real semialgebraic system with CAD/QE-backed checks.

    The first implementation focuses on dependable semantic simplifications:
    contradiction detection, duplicate removal, redundancy removal by implication,
    and one-dimensional inequality normalization via SymPy's real inequality
    reducer. It intentionally returns a mathematically equivalent formula rather
    than promising a globally minimal or prettiest CAD cell union.
    """

    if output not in {"formula", "constraints", "result"}:
        raise ValueError("output must be 'formula', 'constraints', or 'result'")
    expr = _normalize_formula(constraints)
    asm = _normalize_formula(assumptions)
    vars_ = _normalize_variables(variables, sp.And(asm, expr))
    initial_atoms = _ordered_unique_exprs(_conjuncts(expr))
    substitutions = _simple_equality_substitutions(initial_atoms, vars_)
    if substitutions:
        substituted_atoms = _apply_substitutions_to_constraints(
            initial_atoms, substitutions, keep_definitions=not eliminate_equalities
        )
        expr = sp.And(*substituted_atoms) if substituted_atoms else sp.true
    combined = sp.And(asm, expr) if asm is not sp.true and asm != sp.true else expr

    explicit_cyl = None
    try:
        from .cad.cells import extract_explicit_cylindrical_solution

        explicit_cyl = extract_explicit_cylindrical_solution(combined, vars_)
    except Exception:
        explicit_cyl = None

    if combined is sp.false or combined == sp.false:
        result = SimplifiedSystem(
            sp.false,
            (),
            vars_,
            True,
            substitutions=substitutions,
            diagnostics={"reason": "inconsistent"},
        )
        return (
            result
            if (return_result or output == "result")
            else (() if output == "constraints" else sp.false)
        )
    if explicit_cyl is None and not is_satisfiable(combined, vars_, strategy=strategy):
        result = SimplifiedSystem(
            sp.false,
            (),
            vars_,
            True,
            substitutions=substitutions,
            diagnostics={"reason": "inconsistent"},
        )
        return (
            result
            if (return_result or output == "result")
            else (() if output == "constraints" else sp.false)
        )
    if explicit_cyl is not None and asm in (sp.true, True):
        formula = explicit_cyl.as_formula(closed=False)
        atoms = _ordered_unique_exprs(_conjuncts(formula))
        result = SimplifiedSystem(
            formula=formula,
            constraints=atoms,
            variables=vars_,
            inconsistent=False,
            removed_redundant=(),
            substitutions=substitutions,
            method="explicit_cylindrical_cell_formula",
            diagnostics={
                "used_cylindrical_cell_formula": True,
                "cell_count": len(explicit_cyl.cells),
            },
        )
        return (
            result
            if (return_result or output == "result")
            else (atoms if output == "constraints" else formula)
        )

    normalized_1d = _try_reduce_inequalities_1d(expr, vars_)
    if normalized_1d is not None and equivalent(expr, normalized_1d, vars_, strategy=strategy):
        expr = normalized_1d

    try:
        expr = simplify_boole(expr, vars_, assumptions=asm, strategy=strategy)
    except Exception:
        try:
            expr = sp.simplify_logic(expr, form="dnf")
        except Exception:
            pass

    atoms = _ordered_unique_exprs(_conjuncts(expr))
    if not atoms:
        result = SimplifiedSystem(
            sp.true, (), vars_, False, substitutions=substitutions, method="trivial"
        )
        if return_result or output == "result":
            return result
        return () if output == "constraints" else sp.true

    kept: list[sp.Expr] = list(atoms)
    removed: list[sp.Expr] = []
    changed = True
    while changed:
        changed = False
        for atom in list(kept):
            others = [other for other in kept if other is not atom]
            if not others:
                continue
            premise = sp.And(*(others + ([] if asm is sp.true or asm == sp.true else [asm])))
            try:
                if implies(premise, atom, vars_, strategy=strategy):
                    kept.remove(atom)
                    removed.append(atom)
                    changed = True
                    break
            except Exception:
                continue

    formula = sp.And(*kept) if kept else sp.true
    try:
        formula = simplify_boole(formula, vars_, assumptions=asm, strategy=strategy)
    except Exception:
        try:
            formula = sp.simplify_logic(formula, form="dnf")
        except Exception:
            pass
    cylindrical_formula_used = False
    if len(vars_) > 1:
        try:
            from .cad.cells import extract_cylindrical_solution

            cyl = extract_cylindrical_solution(formula, vars_, selected_only=True)
            # Use a finite cylindrical formula only when it is small enough to
            # improve readability rather than explode the output. This lets the
            # simplifier exploit arbitrary-dimensional nested CAD cells while
            # staying conservative for large decompositions.
            if cyl.cells and len(cyl.cells) <= 12:
                cyl_formula = cyl.as_formula(closed=False)
                if explicit_cyl is not None:
                    formula = cyl_formula
                    kept = list(_conjuncts(formula))
                    cylindrical_formula_used = True
                elif equivalent(formula, cyl_formula, vars_, strategy=strategy):
                    formula = cyl_formula
                    kept = list(_conjuncts(formula))
                    cylindrical_formula_used = True
        except Exception:
            pass

    result = SimplifiedSystem(
        formula=formula,
        constraints=tuple(kept),
        variables=vars_,
        inconsistent=False,
        removed_redundant=tuple(removed),
        substitutions=substitutions,
        diagnostics={
            "input_constraint_count": len(atoms),
            "kept_constraint_count": len(kept),
            "used_cylindrical_cell_formula": cylindrical_formula_used,
        },
    )
    if return_result or output == "result":
        return result
    if output == "constraints":
        return tuple(kept)
    return formula


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
    raise ValueError(f"unknown sign proof relation: {relation!r}")


def _constant_sign_certificate(expr: sp.Expr, relation: str) -> bool | None:
    try:
        value = sp.simplify(expr)
    except Exception:
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
    except Exception:
        return None
    return None


def _is_even_power_nonnegative(expr: sp.Expr) -> bool:
    if expr == 0:
        return True
    base, exp = expr.as_base_exp()
    return bool(exp.is_integer and exp.is_even and exp.is_nonnegative)


def _is_obvious_square_product_nonnegative(expr: sp.Expr, variables: Sequence[sp.Symbol]) -> bool:
    try:
        factored = sp.factor(expr)
    except Exception:
        factored = expr
    if factored == 0:
        return True
    if _is_even_power_nonnegative(factored):
        return True
    coeff, factors = factored.as_coeff_mul()
    try:
        if not bool(coeff >= 0):
            return False
    except Exception:
        return False
    for factor in factors:
        if _is_even_power_nonnegative(factor):
            continue
        if getattr(factor, "free_symbols", set()):
            return False
        try:
            if not bool(factor >= 0):
                return False
        except Exception:
            return False
    return True


def _is_sum_of_squares_nonnegative(expr: sp.Expr, variables: Sequence[sp.Symbol]) -> bool:
    try:
        expanded = sp.expand(expr)
    except Exception:
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
    except Exception:
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
        except Exception:
            pass
    return None


def _prove_sign(
    expr: sp.Expr,
    relation: str,
    variables: Sequence[sp.Symbol | str] | None = None,
    *,
    assumptions: FormulaLike | Iterable[FormulaLike] = True,
    strategy: str | None = None,
    return_result: bool = False,
) -> bool | SignProofResult:
    expression = sp.sympify(expr)
    formula = _relation_formula(expression, relation)
    asm = _normalize_formula(assumptions)
    vars_ = _normalize_variables(variables, sp.And(asm, formula))
    cheap = _cheap_sign_proof(expression, relation, asm, vars_)
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
    return_result: bool = False,
) -> bool | SignProofResult:
    return _prove_sign(
        expr,
        "positive",
        variables,
        assumptions=assumptions,
        strategy=strategy,
        return_result=return_result,
    )


def prove_nonnegative(
    expr: sp.Expr,
    variables: Sequence[sp.Symbol | str] | None = None,
    *,
    assumptions: FormulaLike | Iterable[FormulaLike] = True,
    strategy: str | None = None,
    return_result: bool = False,
) -> bool | SignProofResult:
    return _prove_sign(
        expr,
        "nonnegative",
        variables,
        assumptions=assumptions,
        strategy=strategy,
        return_result=return_result,
    )


def prove_negative(
    expr: sp.Expr,
    variables: Sequence[sp.Symbol | str] | None = None,
    *,
    assumptions: FormulaLike | Iterable[FormulaLike] = True,
    strategy: str | None = None,
    return_result: bool = False,
) -> bool | SignProofResult:
    return _prove_sign(
        expr,
        "negative",
        variables,
        assumptions=assumptions,
        strategy=strategy,
        return_result=return_result,
    )


def prove_nonpositive(
    expr: sp.Expr,
    variables: Sequence[sp.Symbol | str] | None = None,
    *,
    assumptions: FormulaLike | Iterable[FormulaLike] = True,
    strategy: str | None = None,
    return_result: bool = False,
) -> bool | SignProofResult:
    return _prove_sign(
        expr,
        "nonpositive",
        variables,
        assumptions=assumptions,
        strategy=strategy,
        return_result=return_result,
    )


def _radial_ball_radius(condition: sp.Expr, variables: Sequence[sp.Symbol]) -> sp.Expr | None:
    atoms = _conjuncts(condition)
    if len(atoms) != 1 or not getattr(atoms[0], "is_Relational", False):
        return None
    atom = atoms[0]
    if not isinstance(atom, (sp.StrictLessThan, sp.LessThan, sp.StrictGreaterThan, sp.GreaterThan)):
        return None
    expr = sp.expand(atom.lhs - atom.rhs)
    try:
        poly = sp.Poly(expr, *variables)
    except Exception:
        return None
    if not variables:
        return None
    coeff = sp.simplify(poly.coeff_monomial(variables[0] ** 2))
    if coeff == 0:
        return None
    for var in variables:
        if sp.simplify(poly.coeff_monomial(var**2) - coeff) != 0:
            return None
    allowed = {
        tuple(2 if i == j else 0 for i in range(len(variables))) for j in range(len(variables))
    }
    allowed.add(tuple(0 for _ in variables))
    if set(poly.monoms()) - allowed:
        return None
    radius_sq = sp.simplify(-poly.coeff_monomial(1) / coeff)
    if coeff.is_positive and isinstance(atom, (sp.StrictLessThan, sp.LessThan)):
        return radius_sq
    if coeff.is_negative and isinstance(atom, (sp.StrictGreaterThan, sp.GreaterThan)):
        return radius_sq
    return None


def region_subset(
    lhs: FormulaLike | Iterable[FormulaLike],
    rhs: FormulaLike | Iterable[FormulaLike],
    variables: Sequence[sp.Symbol | str] | None = None,
    *,
    strategy: str | None = None,
) -> bool:
    left = _normalize_formula(lhs)
    right = _normalize_formula(rhs)
    vars_ = _normalize_variables(variables, sp.And(left, right))
    left_radius = _radial_ball_radius(left, vars_)
    right_radius = _radial_ball_radius(right, vars_)
    if left_radius is not None and right_radius is not None:
        return bool(sp.simplify(left_radius <= right_radius))
    return implies(left, right, vars_, strategy=strategy)


def region_equal(
    lhs: FormulaLike | Iterable[FormulaLike],
    rhs: FormulaLike | Iterable[FormulaLike],
    variables: Sequence[sp.Symbol | str] | None = None,
    *,
    strategy: str | None = None,
) -> bool:
    left = _normalize_formula(lhs)
    right = _normalize_formula(rhs)
    vars_ = _normalize_variables(variables, sp.And(left, right))
    return equivalent(left, right, vars_, strategy=strategy)


def region_disjoint(
    lhs: FormulaLike | Iterable[FormulaLike],
    rhs: FormulaLike | Iterable[FormulaLike],
    variables: Sequence[sp.Symbol | str] | None = None,
    *,
    strategy: str | None = None,
) -> bool:
    left = _normalize_formula(lhs)
    right = _normalize_formula(rhs)
    vars_ = _normalize_variables(variables, sp.And(left, right))
    return not is_satisfiable(sp.And(left, right), vars_, strategy=strategy)


def _linear_box_bounds(
    condition: sp.Expr, variables: Sequence[sp.Symbol]
) -> dict[sp.Symbol, tuple[sp.Expr, sp.Expr]] | None:
    atoms = _conjuncts(condition)
    if not atoms:
        return None
    bounds = {var: (-sp.oo, sp.oo) for var in variables}
    for atom in atoms:
        if not getattr(atom, "is_Relational", False):
            return None
        if isinstance(atom, (sp.Equality, sp.Unequality)):
            return None
        lhs, rhs = atom.lhs, atom.rhs
        expr = sp.expand(lhs - rhs)
        involved = [var for var in variables if var in expr.free_symbols]
        if len(involved) != 1:
            return None
        var = involved[0]
        try:
            poly = sp.Poly(expr, var)
        except Exception:
            return None
        if poly.degree() > 1 or expr.free_symbols - {var}:
            return None
        coeff = sp.simplify(poly.coeff_monomial(var))
        rest = sp.simplify(poly.as_expr() - coeff * var)
        if coeff == 0:
            continue
        point = sp.simplify(-rest / coeff)
        relation = atom
        if coeff < 0:
            if isinstance(atom, sp.StrictGreaterThan):
                relation = sp.StrictLessThan(var, point)
            elif isinstance(atom, sp.GreaterThan):
                relation = sp.LessThan(var, point)
            elif isinstance(atom, sp.StrictLessThan):
                relation = sp.StrictGreaterThan(var, point)
            elif isinstance(atom, sp.LessThan):
                relation = sp.GreaterThan(var, point)
            else:
                return None
        else:
            if isinstance(atom, sp.StrictGreaterThan):
                relation = sp.StrictGreaterThan(var, point)
            elif isinstance(atom, sp.GreaterThan):
                relation = sp.GreaterThan(var, point)
            elif isinstance(atom, sp.StrictLessThan):
                relation = sp.StrictLessThan(var, point)
            elif isinstance(atom, sp.LessThan):
                relation = sp.LessThan(var, point)
            else:
                return None
        lo, hi = bounds[var]
        if isinstance(relation, (sp.StrictGreaterThan, sp.GreaterThan)):
            lo = point if lo == -sp.oo or bool(point > lo) else lo
        elif isinstance(relation, (sp.StrictLessThan, sp.LessThan)):
            hi = point if hi == sp.oo or bool(point < hi) else hi
        bounds[var] = (lo, hi)
    return bounds


def _radial_upper_bound(condition: sp.Expr, variables: Sequence[sp.Symbol]) -> bool | None:
    if len(variables) == 0:
        return True
    norm = sp.Add(*[var**2 for var in variables])
    for atom in _conjuncts(condition):
        if not getattr(atom, "is_Relational", False):
            continue
        if isinstance(atom, (sp.StrictLessThan, sp.LessThan, sp.StrictGreaterThan, sp.GreaterThan)):
            expr = sp.expand(atom.lhs - atom.rhs)
            if sp.simplify(expr - (norm - 1)) == 0 and isinstance(
                atom, (sp.StrictLessThan, sp.LessThan)
            ):
                return True
            if sp.simplify(expr + (norm - 1)) == 0 and isinstance(
                atom, (sp.StrictGreaterThan, sp.GreaterThan)
            ):
                return True
    return None


def _unbounded_by_qe(condition: sp.Expr, variables: Sequence[sp.Symbol]) -> bool | None:
    if not variables:
        return False
    radius = sp.Symbol("_semialg_radius_bound", real=True)
    norm_sq = sp.Add(*[var**2 for var in variables])
    # Unbounded iff for every positive radius threshold there is a point in the
    # region outside that squared-radius threshold.
    matrix = sp.Or(radius <= 0, sp.And(condition, norm_sq > radius))
    quantifiers = (("forall", radius),) + tuple(("exists", var) for var in variables)
    try:
        result = qe_by_complete_cad((radius, *variables), quantifiers, parse_formula(matrix))
    except Exception:
        return None
    if result.is_sentence:
        return bool(result.truth_value)
    reduced = sp.simplify(result.formula)
    if reduced is sp.true or reduced == sp.true:
        return True
    if reduced is sp.false or reduced == sp.false:
        return False
    return None


def region_bounded(
    region: FormulaLike | Iterable[FormulaLike],
    variables: Sequence[sp.Symbol | str] | None = None,
    *,
    strategy: str | None = None,
) -> bool:
    """Return whether a semialgebraic region is bounded over the reals.

    The implementation first recognizes common interval/box/radial cases, then
    falls back to a CAD/QE sentence expressing unboundedness.
    """

    expr = _normalize_formula(region)
    vars_ = _normalize_variables(variables, expr)
    if not is_satisfiable(expr, vars_, strategy=strategy):
        return True
    if len(vars_) == 0:
        return True
    if len(vars_) == 1:
        try:
            reduced = sp.reduce_inequalities(list(_conjuncts(expr)), vars_[0])
            bounds = _linear_box_bounds(reduced, vars_) or _linear_box_bounds(expr, vars_)
        except Exception:
            bounds = _linear_box_bounds(expr, vars_)
        if bounds is not None:
            lo, hi = bounds[vars_[0]]
            return lo != -sp.oo and hi != sp.oo
    box = _linear_box_bounds(expr, vars_)
    if box is not None:
        return all(lo != -sp.oo and hi != sp.oo for lo, hi in box.values())
    radial = _radial_upper_bound(expr, vars_)
    if radial is True:
        return True
    unbounded = _unbounded_by_qe(expr, vars_)
    if unbounded is not None:
        return not unbounded
    raise NotImplementedError("could not determine boundedness for this semialgebraic region")


def region_closed(
    region: FormulaLike | Iterable[FormulaLike],
    variables: Sequence[sp.Symbol | str] | None = None,
    *,
    strategy: str | None = None,
) -> bool:
    expr = _normalize_formula(region)
    vars_ = _normalize_variables(variables, expr)
    return equivalent(expr, region_closure(expr, vars_), vars_, strategy=strategy)


def region_compact(
    region: FormulaLike | Iterable[FormulaLike],
    variables: Sequence[sp.Symbol | str] | None = None,
    *,
    strategy: str | None = None,
) -> bool:
    expr = _normalize_formula(region)
    vars_ = _normalize_variables(variables, expr)
    return region_bounded(expr, vars_, strategy=strategy) and region_closed(
        expr, vars_, strategy=strategy
    )


def _dominates(
    candidate: sp.Expr,
    others: Sequence[sp.Expr],
    assumptions: sp.Expr,
    variables: Sequence[sp.Symbol],
    *,
    weak: bool,
    strategy: str | None,
) -> bool:
    rels = [
        (candidate >= other if weak else candidate <= other)
        for other in others
        if other != candidate
    ]
    return all(implies(assumptions, rel, variables, strategy=strategy) for rel in rels)


def _provable(
    condition: sp.Expr, assumptions: sp.Expr, variables: Sequence[sp.Symbol], strategy: str | None
) -> bool:
    """Return True when ``assumptions => condition`` can be proved."""

    if condition in (True, sp.true):
        return True
    if condition in (False, sp.false):
        return False
    try:
        return bool(implies(assumptions, condition, variables, strategy=strategy))
    except Exception:
        try:
            return bool(sp.simplify(condition) is sp.true or sp.simplify(condition) == sp.true)
        except Exception:
            return False


def _simplify_abs(
    expr: sp.Abs, assumptions: sp.Expr, variables: Sequence[sp.Symbol], strategy: str | None
) -> sp.Expr:
    arg = expr.args[0]
    if _provable(arg >= 0, assumptions, variables, strategy):
        return arg
    if _provable(arg <= 0, assumptions, variables, strategy):
        return -arg
    return expr


def _sqrt_square_base(base: sp.Expr) -> sp.Expr | None:
    """Return ``g`` when ``base`` is syntactically a square ``g**2``.

    This deliberately recognizes only identities that are valid over the real
    domain without changing the domain of the expression. The caller then
    simplifies ``Abs(g)`` under the active assumptions.
    """

    base = sp.factor(base)
    if isinstance(base, sp.Pow) and base.exp == 2:
        return base.base
    if isinstance(base, sp.Mul):
        coeff, factors = base.as_coeff_mul()
        square_root_coeff = None
        if coeff != 1:
            root = sp.sqrt(coeff)
            if root.is_Rational or (root.is_Integer if hasattr(root, "is_Integer") else False):
                square_root_coeff = root
            elif sp.simplify(root**2 - coeff) == 0 and root.is_real is not False:
                square_root_coeff = root
            else:
                return None
        pieces: list[sp.Expr] = []
        if square_root_coeff not in (None, 1):
            pieces.append(square_root_coeff)
        for factor in factors:
            if isinstance(factor, sp.Pow) and factor.exp.is_Integer and int(factor.exp) % 2 == 0:
                pieces.append(factor.base ** (int(factor.exp) // 2))
            else:
                return None
        return sp.Mul(*pieces) if pieces else sp.Integer(1)
    return None


def _simplify_power_under_assumptions(
    node: sp.Pow,
    assumptions: sp.Expr,
    variables: Sequence[sp.Symbol],
    strategy: str | None,
) -> sp.Expr | None:
    if node.exp == sp.Rational(1, 2):
        squared = _sqrt_square_base(node.base)
        if squared is not None:
            return _simplify_abs(sp.Abs(squared), assumptions, variables, strategy)
        return sp.sqrt(node.base)
    return None


def _simplify_log_under_assumptions(
    node: sp.log,
    assumptions: sp.Expr,
    variables: Sequence[sp.Symbol],
    strategy: str | None,
) -> sp.Expr | None:
    arg = node.args[0]
    if arg.func == sp.exp:
        inner = arg.args[0]
        if all(sym.is_real is not False for sym in inner.free_symbols):
            return inner
    if isinstance(arg, sp.Pow) and arg.exp == 2:
        base = arg.base
        if _provable(base > 0, assumptions, variables, strategy):
            return 2 * sp.log(base)
        if _provable(base < 0, assumptions, variables, strategy):
            return 2 * sp.log(-base)
    return None


def _proved_nonzero_denominator(
    denom: sp.Expr,
    assumptions: sp.Expr,
    variables: Sequence[sp.Symbol],
    strategy: str | None,
) -> bool:
    if denom == 1:
        return True
    if denom.is_number:
        try:
            return bool(denom != 0)
        except Exception:
            return False
    return _provable(sp.Ne(denom, 0), assumptions, variables, strategy)


def _safe_cancel_under_assumptions(
    expr: sp.Expr,
    assumptions: sp.Expr,
    variables: Sequence[sp.Symbol],
    strategy: str | None,
    conditions: list[sp.Expr],
    rewrites: list[str],
    *,
    allow_conditional: bool,
) -> sp.Expr:
    try:
        num, den = sp.fraction(sp.together(expr))
        cancelled = sp.cancel(expr)
        cnum, cden = sp.fraction(sp.together(cancelled))
    except Exception:
        return expr
    if cancelled == expr or den == cden:
        return expr
    # ``cancel`` is safe on the domain where the original denominator is nonzero.
    condition = sp.Ne(den, 0)
    if _proved_nonzero_denominator(den, assumptions, variables, strategy):
        rewrites.append("cancel_proved_nonzero_denominator")
        return cancelled
    if allow_conditional:
        conditions.append(condition)
        rewrites.append("cancel_with_side_condition")
        return cancelled
    return expr


def _safe_scalar_simplify(expr: sp.Expr) -> sp.Expr:
    """Apply only conservative scalar cleanup.

    Do not call ``sp.simplify`` here: for rational expressions it may cancel
    denominators and silently change the domain. Domain-changing cancellation
    is handled explicitly by ``_safe_cancel_under_assumptions``.
    """

    if isinstance(expr, Boolean):
        return expr
    try:
        return sp.factor_terms(expr)
    except Exception:
        return expr


def simplify_under_assumptions(
    expr: sp.Expr,
    assumptions: FormulaLike | Iterable[FormulaLike] = True,
    variables: Sequence[sp.Symbol | str] | None = None,
    *,
    strategy: str | None = None,
    return_conditions: bool = False,
    return_result: bool = False,
) -> sp.Expr | AssumptionSimplificationResult:
    """Simplify real expressions using provable assumptions.

    Supported rewrites include sign-sensitive ``Abs`` and ``sqrt(square)``
    simplification, branch-wise ``Piecewise`` simplification, ``Min``/``Max``
    dominance, real ``log(exp(x))``, positive-domain ``log(x**2)``, and safe
    rational cancellation when the original denominator is provably nonzero.

    By default, rewrites that need extra side conditions are not applied. With
    ``return_conditions=True`` or ``return_result=True``, such conditional
    rewrites may be returned together with their required side conditions.
    """

    expression = sp.sympify(expr)
    asm = _normalize_formula(assumptions)
    vars_ = _normalize_variables(
        variables, sp.And(asm, expression >= expression) if expression.free_symbols else asm
    )
    conditions: list[sp.Expr] = []
    rewrites: list[str] = []
    allow_conditional = bool(return_conditions or return_result)

    def rec(node: sp.Expr) -> sp.Expr:
        if isinstance(node, sp.Piecewise):
            rebuilt = sp.Piecewise(*[(rec(value), cond) for value, cond in node.args])
            result = simplify_piecewise(
                rebuilt, vars_, assumptions=asm, strategy=strategy, return_result=True
            )
            if getattr(result, "simplified_branch_values", 0):
                rewrites.append("piecewise_branch_values")
            return result.expression
        if isinstance(node, sp.Abs):
            simplified = _simplify_abs(sp.Abs(rec(node.args[0])), asm, vars_, strategy)
            if simplified != node:
                rewrites.append("abs_sign")
            return simplified
        if isinstance(node, sp.Max):
            args = tuple(rec(arg) for arg in node.args)
            for candidate in args:
                if all(
                    _provable(candidate >= other, asm, vars_, strategy)
                    for other in args
                    if other != candidate
                ):
                    rewrites.append("max_dominance")
                    return candidate
            return sp.Max(*args)
        if isinstance(node, sp.Min):
            args = tuple(rec(arg) for arg in node.args)
            for candidate in args:
                if all(
                    _provable(candidate <= other, asm, vars_, strategy)
                    for other in args
                    if other != candidate
                ):
                    rewrites.append("min_dominance")
                    return candidate
            return sp.Min(*args)
        if isinstance(node, sp.Pow):
            base = rec(node.base)
            rebuilt = sp.Pow(base, node.exp, evaluate=False)
            simplified_pow = _simplify_power_under_assumptions(rebuilt, asm, vars_, strategy)
            if simplified_pow is not None and simplified_pow != rebuilt:
                rewrites.append("sqrt_square")
                return simplified_pow
            return rebuilt
        if node.func == sp.log:
            arg = rec(node.args[0])
            rebuilt = sp.log(arg)
            simplified_log = _simplify_log_under_assumptions(rebuilt, asm, vars_, strategy)
            if simplified_log is not None:
                rewrites.append("log_domain")
                return simplified_log
            return rebuilt
        if not node.args:
            return node
        try:
            new_args = tuple(rec(arg) if isinstance(arg, sp.Expr) else arg for arg in node.args)
            rebuilt = node.func(*new_args)
            return _safe_scalar_simplify(rebuilt)
        except Exception:
            return node

    simplified = rec(expression)
    simplified = _safe_cancel_under_assumptions(
        simplified,
        asm,
        vars_,
        strategy,
        conditions,
        rewrites,
        allow_conditional=allow_conditional,
    )
    simplified = _safe_scalar_simplify(simplified)
    if return_conditions or return_result:
        return AssumptionSimplificationResult(
            expression=simplified,
            original=expression,
            assumptions=asm,
            variables=vars_,
            conditions=tuple(dict.fromkeys(conditions)),
            rewrites=tuple(dict.fromkeys(rewrites)),
            diagnostics={"conditional_rewrites_allowed": allow_conditional},
        )
    return simplified


__all__ = [
    "SimplifiedSystem",
    "AssumptionSimplificationResult",
    "SignProofResult",
    "simplify_system",
    "prove_positive",
    "prove_nonnegative",
    "prove_negative",
    "prove_nonpositive",
    "region_subset",
    "region_equal",
    "region_disjoint",
    "region_bounded",
    "region_closed",
    "region_compact",
    "simplify_under_assumptions",
]
