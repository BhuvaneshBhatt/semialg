from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass

import sympy as sp
from sympy.logic.boolalg import Boolean

from ..cad_algorithms.decomposition import CompleteCAD
from ..cad_algorithms.lifting.stack import CADCell
from ..cad_algorithms.projection.collins import ProjectionPolynomial
from ..formula import And, Atom, BoolConst, Formula, Not, Or
from ..normalization import normalize_formula, normalize_variables
from ..reconstruct.cylindrical import path_condition
from ..simplify.formula import simplify_qe_formula
from ..topology.operations import boundary_cells, cells_formula, interior_cells


@dataclass(frozen=True)
class ParametricCADBoundaryCause:
    """A polynomial whose zero set can carry nongeneric CAD behavior.

    The source string is small and stable so callers and tests can
    distinguish input boundaries from projection-generated degeneracy loci.
    """

    polynomial: sp.Poly
    source: str
    variables: tuple[sp.Symbol, ...]
    parents: tuple[str, ...] = ()

    def equation(self) -> sp.Equality:
        return sp.Eq(sp.expand(self.polynomial.as_expr()), 0)


def _normalize_poly(expr: sp.Expr, variables: Sequence[sp.Symbol]) -> sp.Poly | None:
    if not variables:
        return None
    try:
        poly = sp.Poly(sp.expand(expr), *variables)
    except (ArithmeticError, TypeError, ValueError, NotImplementedError, sp.PolynomialError):
        return None
    if poly.is_zero or poly.total_degree() == 0:
        return None
    primitive = poly.primitive()[1]
    if primitive.LC().could_extract_minus_sign():
        primitive = -primitive
    return primitive


def _atom_boundary_exprs(formula: Formula) -> tuple[sp.Expr, ...]:
    if isinstance(formula, Atom):
        if formula.op in {"=", "<", "<=", ">", ">=", "!="}:
            return (sp.expand(formula.expr),)
        return tuple()
    if isinstance(formula, BoolConst):
        return tuple()
    if isinstance(formula, (And, Or)):
        out: list[sp.Expr] = []
        for arg in formula.args:
            out.extend(_atom_boundary_exprs(arg))
        return tuple(out)
    if isinstance(formula, Not):
        return _atom_boundary_exprs(formula.arg)
    return tuple()


def input_boundary_causes(
    formula: Formula, variables: Sequence[sp.Symbol]
) -> tuple[ParametricCADBoundaryCause, ...]:
    """Return atom-boundary polynomials from the user formula.

    These are the most important exceptional
    equations: inequalities contribute their boundary, equalities contribute the
    algebraic variety, and disequalities contribute their deleted locus.
    """

    seen: set[sp.Expr] = set()
    causes: list[ParametricCADBoundaryCause] = []
    vars_tuple = tuple(variables)
    for expr in _atom_boundary_exprs(formula):
        poly = _normalize_poly(expr, vars_tuple)
        if poly is None:
            continue
        key = sp.expand(poly.as_expr())
        if key in seen:
            continue
        seen.add(key)
        causes.append(ParametricCADBoundaryCause(poly, "input_boundary", vars_tuple))
    return tuple(causes)


def projection_causes(
    cad: CompleteCAD, variables: Sequence[sp.Symbol]
) -> tuple[ParametricCADBoundaryCause, ...]:
    """Return discriminant/resultant/nullification-style projection causes."""

    vars_tuple = tuple(variables)
    seen: set[tuple[str, sp.Expr]] = set()
    causes: list[ParametricCADBoundaryCause] = []
    for level in cad.tower.levels:
        for entry in level.entries:
            source = _projection_source(entry)
            if source is None:
                continue
            poly = _normalize_poly(entry.poly.as_expr(), vars_tuple[: level.level])
            if poly is None:
                continue
            key = (source, sp.expand(poly.as_expr()))
            if key in seen:
                continue
            seen.add(key)
            causes.append(
                ParametricCADBoundaryCause(
                    polynomial=poly,
                    source=source,
                    variables=tuple(vars_tuple[: level.level]),
                    parents=entry.parents,
                )
            )
    return tuple(causes)


def _projection_source(entry: ProjectionPolynomial) -> str | None:
    source = entry.source
    if "discriminant" in source:
        return "discriminant"
    if "resultant" in source:
        return "resultant"
    if "coefficient" in source or "content" in source:
        return "nullification"
    return None


def exceptional_formula_from_causes(causes: Iterable[ParametricCADBoundaryCause]) -> sp.Expr:
    equations = [cause.equation() for cause in causes]
    if not equations:
        return sp.false
    return simplify_qe_formula(sp.Or(*equations), implication_minimize=False)


def relevant_causes_for_cells(
    causes: Iterable[ParametricCADBoundaryCause],
    exceptional_formula: sp.Expr,
) -> tuple[ParametricCADBoundaryCause, ...]:
    """Filter causes to those plausibly useful for the visible exception set.

    A full semantic minimal hitting set would require an additional QE pass. We
    keep the predicate conservative: when the exceptional cell
    formula is nonempty, retain the atom-boundary causes first, and let callers
    include projection causes in diagnostics for auditability.
    """

    if exceptional_formula is sp.false or exceptional_formula == sp.false:
        return tuple()
    kept = [cause for cause in causes if cause.source == "input_boundary"]
    return tuple(kept or tuple(causes))


@dataclass(frozen=True)
class ParametricCADSplit:
    """Cell-level split used by generic CAD output."""

    generic_cells: tuple[CADCell, ...]
    exceptional_cells: tuple[CADCell, ...]
    exceptional_polys: tuple[sp.Poly, ...]
    exceptional_causes: tuple[ParametricCADBoundaryCause, ...]
    generic_formula: sp.Expr
    exceptional_formula: sp.Expr


def _generic_cell_formula(
    cells: Sequence[CADCell],
    variables: Sequence[sp.Symbol],
    cells_by_level,
) -> sp.Expr:
    pieces = [path_condition(cell, variables, cells_by_level, closed=False) for cell in cells]
    kept = [piece for piece in pieces if piece is not sp.false and piece != sp.false]
    if not kept:
        return sp.false
    return simplify_qe_formula(sp.Or(*kept), implication_minimize=False)


def parametric_split_from_cells(
    selected_cells: Sequence[CADCell],
    cad: CompleteCAD,
    variables: Sequence[sp.Symbol],
    causes: Sequence[ParametricCADBoundaryCause] = (),
) -> ParametricCADSplit:
    """Split a selected CAD cell set into full-dimensional and exceptional parts."""

    vars_tuple = tuple(variables)
    generic_cells = interior_cells(selected_cells, cad, vars_tuple)
    exceptional_cells = boundary_cells(selected_cells, cad, vars_tuple)
    generic_formula = _generic_cell_formula(generic_cells, vars_tuple, cad.cells_by_level)
    cell_exceptional = cells_formula(exceptional_cells, vars_tuple, cad.cells_by_level)
    relevant = relevant_causes_for_cells(causes, cell_exceptional)
    cause_exceptional = exceptional_formula_from_causes(relevant)
    exceptional_formula = cause_exceptional if cause_exceptional != sp.false else cell_exceptional
    return ParametricCADSplit(
        generic_cells=tuple(generic_cells),
        exceptional_cells=tuple(exceptional_cells),
        exceptional_polys=tuple(cause.polynomial for cause in relevant),
        exceptional_causes=tuple(relevant),
        generic_formula=simplify_qe_formula(generic_formula, implication_minimize=False),
        exceptional_formula=simplify_qe_formula(exceptional_formula, implication_minimize=False),
    )


FormulaLike = sp.Expr | Boolean | bool


@dataclass(frozen=True)
class ParametricCADBoundaryPolynomial:
    polynomial: sp.Expr
    source: str
    variable: sp.Symbol | None = None
    parents: tuple[sp.Expr, ...] = ()

    def branches(self) -> tuple[sp.Expr, sp.Expr, sp.Expr]:
        p = sp.expand(self.polynomial)
        return (p < 0, sp.Eq(p, 0), p > 0)


@dataclass(frozen=True)
class ParametricCADBoundaryAnalysis:
    variables: tuple[sp.Symbol, ...]
    parameters: tuple[sp.Symbol, ...]
    causes: tuple[ParametricCADBoundaryPolynomial, ...]

    @property
    def exceptional_condition(self) -> sp.Expr:
        if not self.causes:
            return sp.false
        return sp.Or(*(sp.Eq(c.polynomial, 0) for c in self.causes))

    @property
    def branching_polynomials(self) -> tuple[sp.Expr, ...]:
        return tuple(c.polynomial for c in self.causes)


def analyze_parametric_boundaries(
    formula: FormulaLike | Iterable[FormulaLike],
    variables: Sequence[sp.Symbol | str],
    parameters: Sequence[sp.Symbol | str],
) -> ParametricCADBoundaryAnalysis:
    """Discover parameter polynomials where algebraic problem type can change.

    Causes include leading coefficients (degree drops), coefficient sign
    boundaries, discriminants (multiplicity/root-count changes), and pairwise
    resultants (root collisions/common factors).  Only parameter-only
    polynomials are retained, so every reported branch is a valid parameter
    stratum boundary.
    """

    expr = normalize_formula(formula)
    params = normalize_variables(parameters, expr, append_context_symbols=False)
    vars_ = normalize_variables(variables, expr, append_context_symbols=False, exclude=params)

    def relational_atoms(node: sp.Expr) -> tuple[sp.Expr, ...]:
        if getattr(node, "is_Relational", False):
            return (node,)
        if isinstance(node, (sp.And, sp.Or)):
            return tuple(atom for arg in node.args for atom in relational_atoms(arg))
        if isinstance(node, sp.Not):
            return relational_atoms(node.args[0])
        return tuple()

    polys = [sp.expand(atom.lhs - atom.rhs) for atom in relational_atoms(expr)]
    causes: list[ParametricCADBoundaryPolynomial] = []
    seen: set[tuple[str, str]] = set()

    def add(poly_expr: sp.Expr, source: str, variable=None, parents=()):
        candidate = sp.expand(poly_expr)
        if candidate == 0 or candidate.free_symbols - set(params):
            return
        try:
            pp = sp.Poly(candidate, *params) if params else sp.Poly(candidate)
            if pp.total_degree() == 0:
                return
            _, primitive = pp.primitive()
            factors = primitive.factor_list()[1]
            factor_polys = [factor for factor, _multiplicity in factors] or [primitive]
        except (sp.PolynomialError, TypeError, ValueError, NotImplementedError):
            return
        for factor_poly in factor_polys:
            factor_expr = sp.expand(factor_poly.as_expr())
            if factor_expr.could_extract_minus_sign():
                factor_expr = -factor_expr
            key = (source, sp.srepr(factor_expr))
            if key not in seen:
                seen.add(key)
                causes.append(
                    ParametricCADBoundaryPolynomial(factor_expr, source, variable, tuple(parents))
                )

    for expr_poly in polys:
        for var in vars_:
            try:
                p = sp.Poly(expr_poly, var)
            except (sp.PolynomialError, TypeError, ValueError):
                continue
            if p.degree() <= 0:
                continue
            add(p.LC(), "degree_drop", var, (expr_poly,))
            for coeff in p.all_coeffs():
                add(coeff, "coefficient_sign", var, (expr_poly,))
            if p.degree() >= 2:
                try:
                    add(sp.discriminant(p.as_expr(), var), "discriminant", var, (expr_poly,))
                except (sp.PolynomialError, TypeError, ValueError):
                    pass

    for i, left in enumerate(polys):
        for right in polys[i + 1 :]:
            for var in vars_:
                if var not in left.free_symbols or var not in right.free_symbols:
                    continue
                try:
                    from ..algebraic.modular import modular_resultant_qq

                    accelerated = modular_resultant_qq(left, right, var)
                    resultant = (
                        accelerated.resultant
                        if accelerated is not None
                        else sp.resultant(left, right, var)
                    )
                    add(resultant, "resultant", var, (left, right))
                except (sp.PolynomialError, TypeError, ValueError):
                    pass
    causes.sort(
        key=lambda c: (
            c.source,
            "" if c.variable is None else c.variable.name,
            sp.srepr(c.polynomial),
        )
    )
    return ParametricCADBoundaryAnalysis(vars_, params, tuple(causes))
