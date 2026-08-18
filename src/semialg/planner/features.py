from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field

import sympy as sp

from ..formula import Formula, ParsedPrenexFormula, equational_constraints, formula_polynomials
from ..heuristics import suggest_var_polys


@dataclass(frozen=True)
class ProblemFeatures:
    variables: tuple[sp.Symbol, ...]
    quantified_variables: tuple[sp.Symbol, ...] = field(default_factory=tuple)
    free_variables: tuple[sp.Symbol, ...] = field(default_factory=tuple)
    num_polynomials: int = 0
    num_atoms: int = 0
    num_quantifiers: int = 0
    quantifier_alternations: int = 0
    max_total_degree: int = 0
    max_variable_degree: int = 0
    equality_count: int = 0
    inequality_count: int = 0
    has_disjunction: bool = False
    has_negation: bool = False
    has_ecs: bool = False
    is_univariate: bool = False
    is_pure_conjunction: bool = True
    suggested_variable_order: tuple[sp.Symbol, ...] = field(default_factory=tuple)

    @property
    def variable_count(self) -> int:
        return len(self.variables)

    @property
    def ec_density(self) -> float:
        return float(self.equality_count) / float(self.num_atoms or 1)

    @property
    def quantified_fraction(self) -> float:
        return float(self.num_quantifiers) / float(self.variable_count or 1)


def _formula_stats(formula: Formula) -> dict[str, int | bool]:
    from ..formula import And, Atom, BoolConst, Not, Or

    if isinstance(formula, BoolConst):
        return {
            "num_atoms": 0,
            "equality_count": 0,
            "inequality_count": 0,
            "has_disjunction": False,
            "has_negation": False,
            "is_pure_conjunction": True,
        }
    if isinstance(formula, Atom):
        return {
            "num_atoms": 1,
            "equality_count": int(formula.op == "="),
            "inequality_count": int(formula.op != "="),
            "has_disjunction": False,
            "has_negation": False,
            "is_pure_conjunction": True,
        }
    if isinstance(formula, Not):
        inner = _formula_stats(formula.arg)
        inner["has_negation"] = True
        inner["is_pure_conjunction"] = False
        return inner
    if isinstance(formula, And):
        out = {
            "num_atoms": 0,
            "equality_count": 0,
            "inequality_count": 0,
            "has_disjunction": False,
            "has_negation": False,
            "is_pure_conjunction": True,
        }
        for arg in formula.args:
            child = _formula_stats(arg)
            for key in ("num_atoms", "equality_count", "inequality_count"):
                out[key] += int(child[key])  # type: ignore[index]
            out["has_disjunction"] = bool(out["has_disjunction"] or child["has_disjunction"])
            out["has_negation"] = bool(out["has_negation"] or child["has_negation"])
            out["is_pure_conjunction"] = bool(
                out["is_pure_conjunction"] and child["is_pure_conjunction"]
            )
        return out
    if isinstance(formula, Or):
        out = {
            "num_atoms": 0,
            "equality_count": 0,
            "inequality_count": 0,
            "has_disjunction": True,
            "has_negation": False,
            "is_pure_conjunction": False,
        }
        for arg in formula.args:
            child = _formula_stats(arg)
            for key in ("num_atoms", "equality_count", "inequality_count"):
                out[key] += int(child[key])  # type: ignore[index]
            out["has_disjunction"] = True
            out["has_negation"] = bool(out["has_negation"] or child["has_negation"])
        return out
    raise TypeError(f"Unsupported formula node: {type(formula)!r}")


def _degree_stats(polys: Sequence[sp.Expr], variables: Sequence[sp.Symbol]) -> tuple[int, int]:
    max_total = 0
    max_var = 0
    for poly in polys:
        try:
            pobj = sp.Poly(poly, *variables) if variables else sp.Poly(poly)
            max_total = max(max_total, int(pobj.total_degree()))
            for var in variables:
                max_var = max(max_var, int(pobj.degree(var)))
        except Exception:
            syms = tuple(sorted(poly.free_symbols, key=lambda s: s.name))
            try:
                pobj = sp.Poly(poly, *syms) if syms else sp.Poly(poly)
                max_total = max(max_total, int(pobj.total_degree()))
                for var in syms:
                    max_var = max(max_var, int(pobj.degree(var)))
            except Exception:
                max_total = max(max_total, 1)
                max_var = max(max_var, 1)
    return max_total, max_var


def quantifier_alternations(quantifiers: Sequence[tuple[str, sp.Symbol]]) -> int:
    if not quantifiers:
        return 0
    count = 0
    prev = quantifiers[0][0]
    for q, _ in quantifiers[1:]:
        if q != prev:
            count += 1
            prev = q
    return count


def extract_problem_features(
    formula: Formula,
    *,
    variables: Iterable[sp.Symbol] | None = None,
    quantifiers: Sequence[tuple[str, sp.Symbol]] = (),
) -> ProblemFeatures:
    if isinstance(formula, tuple) and len(formula) == 2:
        formula = formula[1]
    polys = tuple(formula_polynomials(formula))
    if variables is None:
        variables = sorted(
            {sym for poly in polys for sym in poly.free_symbols}, key=lambda s: s.name
        )
    vars_tuple = tuple(variables)
    stats = _formula_stats(formula)
    max_total, max_var = _degree_stats(polys, vars_tuple)
    quantified = tuple(var for _, var in quantifiers)
    free = tuple(var for var in vars_tuple if var not in quantified)
    suggested = suggest_var_polys(polys, vars_tuple or None, strategy="degree")
    return ProblemFeatures(
        variables=vars_tuple,
        quantified_variables=quantified,
        free_variables=free,
        num_polynomials=len(polys),
        num_atoms=int(stats["num_atoms"]),
        num_quantifiers=len(quantifiers),
        quantifier_alternations=quantifier_alternations(quantifiers),
        max_total_degree=max_total,
        max_variable_degree=max_var,
        equality_count=int(stats["equality_count"]),
        inequality_count=int(stats["inequality_count"]),
        has_disjunction=bool(stats["has_disjunction"]),
        has_negation=bool(stats["has_negation"]),
        has_ecs=bool(equational_constraints(formula)) or int(stats["equality_count"]) > 0,
        is_univariate=len(vars_tuple) <= 1,
        is_pure_conjunction=bool(stats["is_pure_conjunction"]),
        suggested_variable_order=suggested,
    )


def extract_features_parsed(parsed: ParsedPrenexFormula) -> ProblemFeatures:
    return extract_problem_features(
        parsed.matrix, variables=parsed.vars, quantifiers=parsed.quantifiers
    )
