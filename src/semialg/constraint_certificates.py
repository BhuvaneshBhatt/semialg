"""Certified polynomial-constraint consequences for semialgebraic models."""

from __future__ import annotations

from dataclasses import dataclass

import sympy as sp

from .algebraic_geometry import IrreducibleAlgebraicComponent, irreducible_components
from .decision import implies
from .normalization import normalize_formula, normalize_problem_variables
from .semialgebraic_local_geometry import PolynomialConstraintSystem, polynomial_constraints
from .sos_certificates import SOSCertificate, search_sos_certificate, verify_sos_certificate


@dataclass(frozen=True)
class NonnegativeCombinationCertificate:
    target: sp.Expr
    premises: tuple[sp.Expr, ...]
    coefficients: tuple[sp.Expr, ...]
    residual: sp.Expr
    sos_certificate: SOSCertificate | None = None


@dataclass(frozen=True)
class ImpliedPolynomialInequality:
    polynomial: sp.Expr
    relation: str
    formula: sp.Expr
    method: str
    certificate: NonnegativeCombinationCertificate | SOSCertificate | None = None
    certified: bool = True


@dataclass(frozen=True)
class RedundantPolynomialInequality:
    clause_index: int
    constraint_index: int
    constraint: sp.Expr
    implication: ImpliedPolynomialInequality


@dataclass(frozen=True)
class ComponentConstraintDescription:
    component: IrreducibleAlgebraicComponent
    component_formula: sp.Expr
    model_on_component: sp.Expr
    constraints: PolynomialConstraintSystem


def _ge_polynomial(relation: sp.Expr, variables: tuple[sp.Symbol, ...]) -> sp.Expr:
    if isinstance(relation, (sp.GreaterThan, sp.StrictGreaterThan)):
        expr = sp.expand(relation.lhs - relation.rhs)
    elif isinstance(relation, (sp.LessThan, sp.StrictLessThan)):
        expr = sp.expand(relation.rhs - relation.lhs)
    else:
        raise TypeError("expected a polynomial inequality")
    sp.Poly(expr, *variables)
    return expr


def verify_nonnegative_combination_certificate(cert, variables) -> bool:
    """Replay an exact nonnegative-combination and optional SOS certificate.

    The verifier is a trust boundary: malformed typed payloads are rejected
    cleanly before polynomial replay begins.
    """
    if not isinstance(cert, NonnegativeCombinationCertificate):
        return False
    try:
        vars_ = tuple(variables)
        premises = tuple(map(sp.sympify, cert.premises))
        coefficients = tuple(map(sp.sympify, cert.coefficients))
        target = sp.sympify(cert.target)
        stored_residual = sp.sympify(cert.residual)
    except (TypeError, ValueError, sp.SympifyError):
        return False
    if any(not isinstance(variable, sp.Symbol) for variable in vars_):
        return False
    if len(premises) != len(coefficients):
        return False
    if any(coefficient.is_nonnegative is not True for coefficient in coefficients):
        return False
    try:
        residual = sp.expand(
            target
            - sum(
                coefficient * premise
                for coefficient, premise in zip(coefficients, premises, strict=True)
            )
        )
        if not sp.Poly(residual - stored_residual, *vars_).is_zero:
            return False
        if cert.sos_certificate is None:
            return bool(sp.Poly(residual, *vars_).is_zero)
    except (sp.PolynomialError, TypeError, ValueError):
        return False
    return verify_sos_certificate(residual, cert.sos_certificate, vars_)


def nonnegative_combination_certificate(
    target, premises, variables, *, allow_sos=True, sos_certificate=None
):
    """Find an exact certificate target = sum(lambda_i premise_i) + SOS, lambda_i >= 0."""
    vars_ = tuple(variables)
    target = sp.expand(sp.sympify(target))
    ps = tuple(sp.expand(sp.sympify(p)) for p in premises)
    lambdas = sp.symbols(f"_lambda0:{len(ps)}", nonnegative=True)
    remainder = sp.Poly(
        sp.expand(
            target
            - sum(
                coefficient * polynomial
                for coefficient, polynomial in zip(lambdas, ps, strict=True)
            )
        ),
        *vars_,
    )
    equations = tuple(remainder.coeff_monomial(m) for m in remainder.monoms())
    solset = sp.linsolve(equations, lambdas) if lambdas else ()
    for solution in solset:
        if any(v.free_symbols for v in solution) or any(
            v.is_nonnegative is not True for v in solution
        ):
            continue
        residual = sp.expand(target - sum(v * p for v, p in zip(solution, ps, strict=True)))
        cert = NonnegativeCombinationCertificate(target, ps, tuple(solution), residual)
        if verify_nonnegative_combination_certificate(cert, vars_):
            return cert
    if sos_certificate is not None:
        cert = NonnegativeCombinationCertificate(
            target, ps, tuple(sp.S.Zero for _ in ps), target, sos_certificate
        )
        return cert if verify_nonnegative_combination_certificate(cert, vars_) else None
    if allow_sos:
        search = search_sos_certificate(target, vars_)
        if search.certified:
            cert = NonnegativeCombinationCertificate(
                target, ps, tuple(sp.S.Zero for _ in ps), target, search.certificate
            )
            if verify_nonnegative_combination_certificate(cert, vars_):
                return cert
    return None


def implied_polynomial_inequality(region, inequality, variables=None, *, use_certificates=True):
    """Certify that ``region`` implies a polynomial inequality."""
    formula = normalize_formula(region)
    vars_ = normalize_problem_variables(variables, sp.And(formula, inequality))
    target = _ge_polynomial(inequality, vars_)
    if use_certificates and isinstance(inequality, (sp.GreaterThan, sp.LessThan)):
        system = polynomial_constraints(formula, vars_)
        if len(system.clauses) == 1:
            premises = []
            for c in system.clauses[0].constraints:
                if c.relation == "ge":
                    premises.append(c.polynomial)
                elif c.relation == "le":
                    premises.append(-c.polynomial)
            cert = nonnegative_combination_certificate(target, premises, vars_)
            if cert is not None:
                return ImpliedPolynomialInequality(
                    target, "ge", inequality, "nonnegative_combination_sos", cert
                )
    result = implies(formula, inequality, vars_, return_result=True)
    if not bool(result):
        return None
    relation = "gt" if isinstance(inequality, (sp.StrictGreaterThan, sp.StrictLessThan)) else "ge"
    return ImpliedPolynomialInequality(
        target, relation, inequality, "exact_semialgebraic_implication"
    )


def redundant_polynomial_inequalities(region, variables=None):
    """Return inequalities implied by the other constraints in their DNF clause."""
    system = polynomial_constraints(region, variables)
    out = []
    for ci, clause in enumerate(system.clauses):
        atoms = tuple(c.atom for c in clause.constraints)
        for i, c in enumerate(clause.constraints):
            if c.relation not in {"ge", "le", "gt", "lt"}:
                continue
            others = sp.And(*(a for j, a in enumerate(atoms) if j != i))
            proof = implied_polynomial_inequality(others, c.atom, system.variables)
            if proof is not None:
                out.append(RedundantPolynomialInequality(ci, i, c.atom, proof))
    return tuple(out)


def component_constraint_descriptions(region, variables=None, *, max_pieces=None):
    """Describe the real model separately on each certified irreducible closure component."""
    formula = normalize_formula(region)
    vars_ = normalize_problem_variables(variables, formula)
    equalities = tuple(sp.expand(a.lhs - a.rhs) for a in formula.atoms(sp.Equality))
    if not equalities:
        raise ValueError("component descriptions require polynomial equality constraints")
    components = irreducible_components(equalities, vars_, max_pieces=max_pieces)
    out = []
    for component in components:
        cf = sp.And(*(sp.Eq(eq, 0) for eq in component.equations))
        restricted = sp.And(formula, cf)
        out.append(
            ComponentConstraintDescription(
                component, cf, restricted, polynomial_constraints(restricted, vars_)
            )
        )
    return tuple(out)


__all__ = [
    "NonnegativeCombinationCertificate",
    "ImpliedPolynomialInequality",
    "RedundantPolynomialInequality",
    "ComponentConstraintDescription",
    "verify_nonnegative_combination_certificate",
    "nonnegative_combination_certificate",
    "implied_polynomial_inequality",
    "redundant_polynomial_inequalities",
    "component_constraint_descriptions",
]
