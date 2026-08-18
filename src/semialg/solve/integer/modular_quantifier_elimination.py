from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field
from itertools import product

import sympy as sp

from .congruence import (
    ModularSolveResult,
    combine_mod_crt,
    eliminate_one_var,
    norm_mod_form,
    solve_mod_lin_sys,
    solve_mod_poly_sys,
    solve_quant_free_mod_sys,
    solve_quant_mod_sys,
)


@dataclass(frozen=True)
class ModularQuantifierBlock:
    quantifier: str
    variables: tuple[sp.Symbol, ...]


@dataclass(frozen=True)
class ModularQePlan:
    modulus: int
    component_moduli: tuple[int, ...]
    quantifier_blocks: tuple[ModularQuantifierBlock, ...]
    free_variables: tuple[sp.Symbol, ...]
    all_variables: tuple[sp.Symbol, ...]


@dataclass
class ModularQeResult:
    modulus: int
    free_variables: tuple[sp.Symbol, ...]
    quantifier_blocks: tuple[ModularQuantifierBlock, ...]
    formula: sp.Expr = sp.false
    points: list[tuple[int, ...]] = field(default_factory=list)
    method: str = "modular_quantifier_elimination"
    complete: bool = True
    metadata: dict = field(default_factory=dict)


def _coerce_modulus(modulus: int) -> int:
    modulus = int(modulus)
    if modulus <= 0:
        raise ValueError("modulus must be a positive integer")
    return modulus


def _factor_prime_powers(modulus: int) -> list[int]:
    fac = sp.factorint(_coerce_modulus(modulus))
    return [int(p**e) for p, e in sorted(fac.items())]


def _residue_formula_1d(variable: sp.Symbol, residues: Sequence[int], modulus: int) -> sp.Expr:
    residues = sorted(set(int(r) % modulus for r in residues))
    if len(residues) == 0:
        return sp.false
    if len(residues) == modulus:
        return sp.true
    if len(residues) == 1:
        return sp.Eq(variable, residues[0])
    if len(residues) == modulus - 1:
        missing = next(r for r in range(modulus) if r not in residues)
        return sp.Ne(variable, missing)
    return sp.Or(*[sp.Eq(variable, r) for r in residues])


def cart_prod_form(
    points: Sequence[Sequence[int]], variables: Sequence[sp.Symbol], modulus: int
) -> sp.Expr | None:
    variables = tuple(variables)
    pts = sorted(set(tuple(int(v) % modulus for v in pt) for pt in points))
    if not pts:
        return sp.false
    if not variables:
        return sp.true
    coordinate_sets = [sorted(set(pt[i] for pt in pts)) for i in range(len(variables))]
    total = 1
    for s in coordinate_sets:
        total *= len(s)
    if total != len(pts):
        return None
    if set(product(*coordinate_sets)) != set(pts):
        return None
    parts = [
        _residue_formula_1d(v, rs, modulus)
        for v, rs in zip(variables, coordinate_sets, strict=True)
    ]
    return sp.And(*parts) if parts else sp.true


def _factorized_formula(
    points: Sequence[Sequence[int]], variables: Sequence[sp.Symbol], modulus: int
) -> sp.Expr | None:
    variables = tuple(variables)
    pts = sorted(set(tuple(int(v) % modulus for v in pt) for pt in points))
    if not pts:
        return sp.false
    if len(variables) <= 1:
        return (
            _residue_formula_1d(variables[0], [pt[0] for pt in pts], modulus)
            if variables
            else sp.true
        )
    direct = cart_prod_form(pts, variables, modulus)
    if direct is not None:
        return direct

    # Try one-variable case split: x in S_i && formula_i(rest)
    best = None
    for i, var in enumerate(variables):
        buckets = {}
        for pt in pts:
            buckets.setdefault(pt[i], []).append(pt[:i] + pt[i + 1 :])
        score = len(buckets)
        if best is None or score < best[0]:
            best = (score, i, var, buckets)
    if best is None:
        return None
    _, idx, var, buckets = best
    rest_vars = variables[:idx] + variables[idx + 1 :]
    clauses = []
    for residue, subpts in sorted(buckets.items()):
        subf = _factorized_formula(subpts, rest_vars, modulus)
        if subf is None:
            subf = sp.Or(
                *[
                    sp.And(*[sp.Eq(v, a) for v, a in zip(rest_vars, pt, strict=True)])
                    for pt in sorted(set(tuple(map(int, spt)) for spt in subpts))
                ]
            )
        clauses.append(sp.And(sp.Eq(var, int(residue)), subf))
    return sp.Or(*clauses) if len(clauses) > 1 else clauses[0]


def _points_to_formula(
    points: Sequence[Sequence[int]], variables: Sequence[sp.Symbol], modulus: int
) -> sp.Expr:
    variables = tuple(variables)
    points = sorted(set(tuple(int(v) % modulus for v in pt) for pt in points))
    if not points:
        return sp.false
    symbolic = _factorized_formula(points, variables, modulus)
    if symbolic is not None:
        return sp.simplify(symbolic)
    clauses = [sp.And(*[sp.Eq(v, a) for v, a in zip(variables, pt, strict=True)]) for pt in points]
    return sp.simplify(sp.Or(*clauses) if len(clauses) > 1 else clauses[0])


def simp_proj_mod_form(
    formula: sp.Expr, free_variables: Sequence[sp.Symbol], modulus: int, *, max_points: int = 10000
) -> sp.Expr:
    free_variables = tuple(free_variables)
    result = solve_quant_free_mod_sys(formula, free_variables, modulus, max_points=max_points)
    return _points_to_formula(result.points, free_variables, modulus)


def norm_mod_quant_blocks(
    blocks: Sequence[ModularQuantifierBlock | tuple[str, Sequence[sp.Symbol]]],
) -> tuple[ModularQuantifierBlock, ...]:
    normalized = []
    seen: set[sp.Symbol] = set()
    for block in blocks:
        if isinstance(block, ModularQuantifierBlock):
            q = block.quantifier
            vs = tuple(block.variables)
        else:
            q, raw_vs = block
            vs = tuple(raw_vs)
        q = q.lower()
        if q not in {"exists", "forall"}:
            raise ValueError("quantifier must be 'exists' or 'forall'")
        if not vs:
            continue
        for v in vs:
            if v in seen:
                raise ValueError(f"duplicate quantified variable: {v}")
            seen.add(v)
        normalized.append(ModularQuantifierBlock(q, vs))
    return tuple(normalized)


def build_modular_qe_plan(
    modulus: int,
    free_variables: Sequence[sp.Symbol],
    quantified_variables: Sequence[sp.Symbol] | None = None,
    quantifier: str = "exists",
    *,
    quantifier_blocks: Sequence[ModularQuantifierBlock | tuple[str, Sequence[sp.Symbol]]]
    | None = None,
) -> ModularQePlan:
    modulus = _coerce_modulus(modulus)
    free_variables = tuple(free_variables)
    if quantifier_blocks is None:
        quant_blocks = norm_mod_quant_blocks([(quantifier, tuple(quantified_variables or ()))])
    else:
        quant_blocks = norm_mod_quant_blocks(quantifier_blocks)
    all_vars = free_variables + tuple(v for block in quant_blocks for v in block.variables)
    return ModularQePlan(
        modulus=modulus,
        component_moduli=tuple(_factor_prime_powers(modulus)),
        quantifier_blocks=quant_blocks,
        free_variables=free_variables,
        all_variables=all_vars,
    )


def _split_eqs_neqs(expr: sp.Expr):
    atoms = list(expr.args) if isinstance(expr, sp.And) else [expr]
    eqs = [a for a in atoms if isinstance(a, sp.Equality)]
    neqs = [a for a in atoms if isinstance(a, sp.Unequality)]
    others = [a for a in atoms if not isinstance(a, (sp.Equality, sp.Unequality))]
    return eqs, neqs, others


def supports_direct_lin_proj(expr: sp.Expr, variables: Sequence[sp.Symbol]) -> bool:
    eqs, _neqs, others = _split_eqs_neqs(expr)
    if others:
        return False
    if not eqs:
        return False
    return all(sp.Poly(sp.expand(eq.lhs - eq.rhs), *variables).total_degree() <= 1 for eq in eqs)


def _project_points(
    points: Sequence[Sequence[int]], all_vars: Sequence[sp.Symbol], keep_vars: Sequence[sp.Symbol]
) -> list[tuple[int, ...]]:
    keep_vars = tuple(keep_vars)
    idx = [tuple(all_vars).index(v) for v in keep_vars]
    return sorted(set(tuple(int(pt[i]) for i in idx) for pt in points))


def nonenum_block_elim(
    expr: sp.Expr,
    current_variables: Sequence[sp.Symbol],
    block: ModularQuantifierBlock,
    modulus: int,
    *,
    max_points: int,
) -> tuple[sp.Expr | None, dict]:
    current_variables = tuple(current_variables)
    remaining = tuple(v for v in current_variables if v not in set(block.variables))
    meta = {"strategy": None}

    if block.quantifier == "exists":
        try:
            if supports_direct_lin_proj(expr, current_variables):
                eqs, neqs, _others = _split_eqs_neqs(expr)
                lin = solve_mod_lin_sys(eqs, current_variables, modulus, max_points=max_points)
                good = []
                for pt in lin.points:
                    subst = dict(zip(current_variables, pt, strict=True))
                    if all(bool(sp.simplify(n.subs(subst))) for n in neqs):
                        good.append(pt)
                proj = _project_points(good, current_variables, remaining)
                meta["strategy"] = "linear_projection"
                meta["source_method"] = lin.method
                return _points_to_formula(proj, remaining, modulus), meta
        except Exception:
            pass

        try:
            eqs, neqs, others = _split_eqs_neqs(expr)
            if eqs and not others:
                poly = solve_mod_poly_sys(eqs, current_variables, modulus, max_points=max_points)
                good = []
                for pt in poly.points:
                    subst = dict(zip(current_variables, pt, strict=True))
                    if all(bool(sp.simplify(n.subs(subst))) for n in neqs):
                        good.append(pt)
                proj = _project_points(good, current_variables, remaining)
                meta["strategy"] = "quantified_groebner_projection"
                meta["source_method"] = poly.method
                return _points_to_formula(proj, remaining, modulus), meta
        except Exception:
            pass
    return None, meta


def _eliminate_single_form(
    formula: sp.Expr,
    current_variables: Sequence[sp.Symbol],
    block: ModularQuantifierBlock,
    modulus: int,
    *,
    max_points: int,
) -> tuple[sp.Expr, tuple[sp.Symbol, ...], dict]:
    current_variables = tuple(current_variables)
    working = norm_mod_form(formula, modulus)
    metadata = {"block": (block.quantifier, tuple(map(str, block.variables))), "steps": []}

    improved, improved_meta = nonenum_block_elim(
        working, current_variables, block, modulus, max_points=max_points
    )
    if improved is not None:
        remaining = tuple(v for v in current_variables if v not in set(block.variables))
        metadata["steps"].append(
            {"strategy": improved_meta.get("strategy"), "formula": sp.srepr(improved)}
        )
        return sp.simplify(improved), remaining, metadata

    for qvar in reversed(block.variables):
        if block.quantifier == "exists":
            working = eliminate_one_var(
                working, current_variables, qvar, modulus, max_points=max_points
            )
        else:
            neg = sp.Not(working)
            eliminated = eliminate_one_var(
                neg, current_variables, qvar, modulus, max_points=max_points
            )
            working = sp.Not(eliminated)
        working = norm_mod_form(sp.simplify(working), modulus)
        working = simp_proj_mod_form(
            working,
            tuple(v for v in current_variables if v != qvar),
            modulus,
            max_points=max_points,
        )
        current_variables = tuple(v for v in current_variables if v != qvar)
        metadata["steps"].append({"eliminated": str(qvar), "formula": sp.srepr(working)})
    return working, current_variables, metadata


def _apply_blocks_power(
    expr: sp.Expr,
    free_variables: Sequence[sp.Symbol],
    blocks: Sequence[ModularQuantifierBlock],
    modulus: int,
    *,
    max_points: int,
) -> ModularQeResult:
    free_variables = tuple(free_variables)
    blocks = tuple(blocks)
    working = norm_mod_form(expr, modulus)
    current = free_variables + tuple(v for block in blocks for v in block.variables)
    block_meta = []

    # Fast path for a single uniform quantifier block
    if len(blocks) == 1:
        block = blocks[0]
        try:
            base_result = solve_quant_mod_sys(
                working,
                free_variables,
                block.variables,
                modulus,
                quantifier=block.quantifier,
                max_points=max_points,
            )
            return ModularQeResult(
                modulus=modulus,
                free_variables=free_variables,
                quantifier_blocks=blocks,
                formula=_points_to_formula(base_result.points, free_variables, modulus),
                points=base_result.points,
                method=f"uniform_{block.quantifier}",
                complete=base_result.complete,
                metadata={"base_method": base_result.method, **base_result.metadata},
            )
        except Exception:
            pass

    for block in reversed(blocks):
        working, current, meta = _eliminate_single_form(
            working, current, block, modulus, max_points=max_points
        )
        block_meta.append(meta)

    points = solve_quant_free_mod_sys(
        working, free_variables, modulus, max_points=max_points
    ).points
    formula = _points_to_formula(points, free_variables, modulus)
    return ModularQeResult(
        modulus=modulus,
        free_variables=free_variables,
        quantifier_blocks=blocks,
        formula=formula,
        points=points,
        method="recursive_block_elimination",
        complete=True,
        metadata={
            "final_formula": sp.srepr(working),
            "block_elimination_trace": tuple(block_meta),
        },
    )


def solve_quant_decomp(
    expr: sp.Expr,
    free_variables: Sequence[sp.Symbol],
    quantified_variables: Sequence[sp.Symbol],
    modulus: int,
    *,
    quantifier: str = "exists",
    max_points: int = 10000,
) -> ModularSolveResult:
    result = solve_mod_qe_with_blocks(
        expr,
        free_variables,
        [(quantifier, tuple(quantified_variables))],
        modulus,
        max_points=max_points,
    )
    return ModularSolveResult(
        modulus=result.modulus,
        variables=result.free_variables,
        points=result.points,
        formula=result.formula,
        method=result.method,
        complete=result.complete,
        metadata=result.metadata,
    )


def solve_mod_qe_with_blocks(
    expr: sp.Expr,
    free_variables: Sequence[sp.Symbol],
    quantifier_blocks: Sequence[ModularQuantifierBlock | tuple[str, Sequence[sp.Symbol]]],
    modulus: int,
    *,
    max_points: int = 10000,
) -> ModularQeResult:
    plan = build_modular_qe_plan(modulus, free_variables, quantifier_blocks=quantifier_blocks)
    if len(plan.component_moduli) <= 1:
        result = _apply_blocks_power(
            expr, plan.free_variables, plan.quantifier_blocks, plan.modulus, max_points=max_points
        )
        result.metadata["component_moduli"] = plan.component_moduli
        return result

    parts = [
        _apply_blocks_power(
            expr, plan.free_variables, plan.quantifier_blocks, q, max_points=max_points
        )
        for q in plan.component_moduli
    ]
    crt_points = combine_mod_crt(
        [p.points for p in parts], plan.free_variables, plan.component_moduli
    )
    return ModularQeResult(
        modulus=plan.modulus,
        free_variables=plan.free_variables,
        quantifier_blocks=plan.quantifier_blocks,
        formula=_points_to_formula(crt_points, plan.free_variables, plan.modulus),
        points=crt_points,
        method="crt_block_modular_qe",
        complete=all(p.complete for p in parts),
        metadata={
            "component_moduli": plan.component_moduli,
            "component_methods": tuple(p.method for p in parts),
            "component_metadata": tuple(p.metadata for p in parts),
        },
    )


def eliminate_quant_m(
    expr: sp.Expr,
    free_variables: Sequence[sp.Symbol],
    quantified_variables: Sequence[sp.Symbol],
    modulus: int,
    *,
    quantifier: str = "exists",
    max_points: int = 10000,
) -> sp.Expr:
    result = solve_mod_qe_with_blocks(
        expr,
        free_variables,
        [(quantifier, tuple(quantified_variables))],
        modulus,
        max_points=max_points,
    )
    return result.formula


def eliminate_mod_form(
    expr: sp.Expr,
    free_variables: Sequence[sp.Symbol],
    quantifier_blocks: Sequence[ModularQuantifierBlock | tuple[str, Sequence[sp.Symbol]]],
    modulus: int,
    *,
    max_points: int = 10000,
) -> sp.Expr:
    result = solve_mod_qe_with_blocks(
        expr, free_variables, quantifier_blocks, modulus, max_points=max_points
    )
    return result.formula


__all__ = [
    "ModularQuantifierBlock",
    "ModularQePlan",
    "ModularQeResult",
    "norm_mod_quant_blocks",
    "build_modular_qe_plan",
    "simp_proj_mod_form",
    "solve_quant_decomp",
    "solve_mod_qe_with_blocks",
    "eliminate_quant_m",
    "eliminate_mod_form",
]
