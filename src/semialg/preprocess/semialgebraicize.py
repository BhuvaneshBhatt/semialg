from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from enum import Enum

import sympy as sp

from ..formula import Formula, parse_formula, to_sympy
from .auxiliary import AuxiliaryDef, AuxiliaryFactory


class PowerPolicy(str, Enum):
    """Policies for interpreting rational powers during preprocessing."""

    REAL_PRINCIPAL = "real_principal"
    SEMIALGEBRAIC_REAL = "semialgebraic_real"
    STRICT_POLYNOMIAL = "strict_polynomial"


@dataclass(frozen=True)
class PreprocessResult:
    """Result of converting supported real semialgebraic syntax to polynomials."""

    formula: Formula
    sympy_expr: sp.Expr
    aux_vars: tuple[sp.Symbol, ...] = ()
    auxiliary_defs: tuple[AuxiliaryDef, ...] = ()
    substitutions: Mapping[sp.Expr, sp.Expr] = field(default_factory=dict)
    assumptions: tuple[sp.Expr, ...] = ()
    changed: bool = False
    notes: tuple[str, ...] = ()


def _replace_float_literals(expr: sp.Expr) -> tuple[sp.Expr, bool]:
    repl = {atom: sp.nsimplify(atom) for atom in expr.atoms(sp.Float)}
    if not repl:
        return expr, False
    return expr.xreplace(repl), True


def _branch_points_for_abs_arg(arg: sp.Expr) -> tuple[sp.Expr, ...]:
    if arg.is_polynomial(*sorted(arg.free_symbols, key=lambda s: s.name)):
        return (sp.Eq(sp.expand(arg), 0),)
    return tuple()


def _replace_abs(
    expr: sp.Expr, factory: AuxiliaryFactory
) -> tuple[sp.Expr, list[AuxiliaryDef], dict[sp.Expr, sp.Expr], list[str]]:
    aux_defs: list[AuxiliaryDef] = []
    substitutions: dict[sp.Expr, sp.Expr] = {}
    notes: list[str] = []
    current = expr
    # Replace innermost/simplest first; repeated xreplace is deterministic and
    # avoids leaving Abs inside later rational-power constraints.
    for atom in sorted(current.atoms(sp.Abs), key=lambda a: (len(sp.srepr(a)), sp.srepr(a))):
        if atom in substitutions:
            continue
        aux = factory.fresh("abs")
        arg = sp.expand(atom.args[0])
        constraints = (aux >= 0, sp.Eq(aux**2, sp.expand(arg**2)))
        aux_defs.append(
            AuxiliaryDef(
                symbol=aux,
                expression=atom,
                kind="abs",
                constraints=constraints,
                branch_points=_branch_points_for_abs_arg(arg),
            )
        )
        substitutions[atom] = aux
        current = current.xreplace({atom: aux})
        notes.append(f"Replaced {sp.sstr(atom)} by nonnegative auxiliary {sp.sstr(aux)}.")
    return current, aux_defs, substitutions, notes


def _is_supported_rational_power(atom: sp.Pow) -> bool:
    exp = sp.Rational(atom.exp)
    return exp.q > 1 and exp > 0


def _replace_rational_powers(
    expr: sp.Expr,
    factory: AuxiliaryFactory,
    *,
    policy: PowerPolicy,
) -> tuple[sp.Expr, list[AuxiliaryDef], dict[sp.Expr, sp.Expr], list[str]]:
    if policy == PowerPolicy.STRICT_POLYNOMIAL:
        return expr, [], {}, []
    aux_defs: list[AuxiliaryDef] = []
    substitutions: dict[sp.Expr, sp.Expr] = {}
    notes: list[str] = []
    current = expr

    # Largest powers first prevents replacing sqrt(base) inside base^(3/2) when
    # SymPy exposes both forms in an expression tree.
    atoms = sorted(
        [a for a in current.atoms(sp.Pow) if a.exp.is_Rational and _is_supported_rational_power(a)],
        key=lambda a: (-len(sp.srepr(a)), sp.srepr(a)),
    )
    for atom in atoms:
        if atom in substitutions:
            continue
        if not current.has(atom):
            continue
        exp = sp.Rational(atom.exp)
        base = sp.expand(atom.base)
        if base.is_number and base.is_nonnegative:
            continue

        aux = factory.fresh("pow")
        constraints: list[sp.Expr] = []
        if exp.q % 2 == 0:
            # Even roots are real only on a nonnegative base and use the
            # principal nonnegative branch.
            constraints.extend([base >= 0, aux >= 0])
        elif exp.p % 2 == 0:
            # Odd roots exist for all real bases. Even numerators make the
            # whole power nonnegative, e.g. x^(2/3).
            constraints.append(aux >= 0)
        constraints.append(sp.Eq(aux**exp.q, sp.expand(base**exp.p)))
        aux_defs.append(
            AuxiliaryDef(
                symbol=aux,
                expression=atom,
                kind="rational_power",
                constraints=tuple(constraints),
                branch_points=(sp.Eq(base, 0),),
            )
        )
        substitutions[atom] = aux
        current = current.xreplace({atom: aux})
        notes.append(
            f"Replaced rational power {sp.sstr(atom)} by auxiliary {sp.sstr(aux)} "
            f"with {exp.q}-th-power constraint."
        )
    return current, aux_defs, substitutions, notes


def semialgebraicize(
    formula: Formula | sp.Expr,
    *,
    variables: tuple[sp.Symbol, ...] | None = None,
    power_policy: PowerPolicy = PowerPolicy.REAL_PRINCIPAL,
) -> PreprocessResult:
    """Convert supported non-polynomial real semialgebraic syntax to polynomial form.

    Supported syntax is intentionally conservative: ``Abs(poly)`` and positive
    rational powers are replaced with existential auxiliaries plus polynomial
    constraints. The returned matrix is equivalent after existentially
    quantifying ``aux_vars``.
    """

    expr = (
        formula if isinstance(formula, (sp.Basic, sp.logic.boolalg.Boolean)) else to_sympy(formula)
    )
    expr = sp.sympify(expr)
    existing = set(expr.free_symbols) | set(variables or ())
    factory = AuxiliaryFactory(existing)
    notes: list[str] = []
    substitutions: dict[sp.Expr, sp.Expr] = {}
    aux_defs: list[AuxiliaryDef] = []
    changed = False

    current, did_float = _replace_float_literals(expr)
    if did_float:
        changed = True
        notes.append("Rationalized floating-point literals.")

    current, defs, repl, new_notes = _replace_abs(current, factory)
    aux_defs.extend(defs)
    substitutions.update(repl)
    notes.extend(new_notes)
    changed = changed or bool(defs)

    # Rational powers may have Abs bases; run after Abs replacement so all
    # generated constraints are polynomial when input bases are polynomial.
    current, defs, repl, new_notes = _replace_rational_powers(current, factory, policy=power_policy)
    aux_defs.extend(defs)
    substitutions.update(repl)
    notes.extend(new_notes)
    changed = changed or bool(defs)

    assumptions = tuple(cond for aux in aux_defs for cond in aux.constraints)
    matrix = sp.And(current, *assumptions) if assumptions else current
    # Avoid global simplify on Boolean formulas: SymPy may spend a long time
    # trying to prove equivalences between relational expressions. The
    # preprocessing pass only needs to build an equivalent polynomial matrix.
    return PreprocessResult(
        formula=parse_formula(matrix),
        sympy_expr=matrix,
        aux_vars=tuple(aux.symbol for aux in aux_defs),
        auxiliary_defs=tuple(aux_defs),
        substitutions=substitutions,
        assumptions=assumptions,
        changed=changed,
        notes=tuple(notes),
    )


__all__ = ["PowerPolicy", "PreprocessResult", "semialgebraicize"]
