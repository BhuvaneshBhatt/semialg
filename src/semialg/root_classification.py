from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field

import sympy as sp

from .sampling import sample_point


@dataclass(frozen=True)
class RootClassificationCell:
    """One parameter-space condition with a constant real-root count."""

    condition: sp.Expr
    root_count: sp.Expr
    multiplicity_pattern: tuple[int, ...] = ()
    sample: Mapping[sp.Symbol, sp.Expr] | None = None
    diagnostics: Mapping[str, object] = field(default_factory=dict)


@dataclass(frozen=True)
class RootClassificationResult:
    """Real-root classification for a univariate polynomial family."""

    polynomial: sp.Expr
    variable: sp.Symbol
    parameters: tuple[sp.Symbol, ...]
    cells: tuple[RootClassificationCell, ...]
    method: str
    diagnostics: Mapping[str, object] = field(default_factory=dict)

    @property
    def generic_root_count(self) -> sp.Expr | None:
        return self.cells[0].root_count if self.cells else None

    @property
    def generic_multiplicity_pattern(self) -> tuple[int, ...]:
        return self.cells[0].multiplicity_pattern if self.cells else ()


def _normalize_parameters(
    parameters: Sequence[sp.Symbol | str] | None,
    expr: sp.Expr,
    variable: sp.Symbol,
) -> tuple[sp.Symbol, ...]:
    if parameters is not None:
        out: list[sp.Symbol] = []
        seen: set[sp.Symbol] = set()
        for param in parameters:
            sym = sp.Symbol(param, real=True) if isinstance(param, str) else param
            if sym != variable and sym not in seen:
                out.append(sym)
                seen.add(sym)
        return tuple(out)
    return tuple(sorted(expr.free_symbols - {variable}, key=lambda sym: sym.name))


def _multiplicity_pattern(poly: sp.Poly) -> tuple[int, ...]:
    roots = sp.roots(poly.as_expr(), poly.gens[0])
    real_mults: list[int] = []
    for root, mult in roots.items():
        if root.is_real is True or bool(sp.N(sp.im(root), 80) == 0):
            real_mults.append(int(mult))
    if not real_mults:
        try:
            real_mults = [1 for _ in sp.real_roots(poly.as_expr())]
        except Exception:
            real_mults = []
    return tuple(sorted(real_mults))


def _unparameterized(expr: sp.Expr, variable: sp.Symbol) -> RootClassificationResult:
    poly = sp.Poly(expr, variable)
    pattern = _multiplicity_pattern(poly)
    cell = RootClassificationCell(sp.true, sp.Integer(len(pattern)), pattern, sample={})
    return RootClassificationResult(
        sp.expand(expr),
        variable,
        (),
        (cell,),
        "univariate_exact_roots",
        {"degree": poly.degree()},
    )


def _sample_for(
    condition: sp.Expr, parameters: Sequence[sp.Symbol]
) -> Mapping[sp.Symbol, sp.Expr] | None:
    if not parameters:
        return {}
    if condition is sp.true or condition is True:
        return {param: sp.Integer(0) for param in parameters}
    if condition is sp.false or condition is False:
        return None
    return sample_point(condition, parameters, strategy="fallback", strict=False)


def _linear_family(
    poly: sp.Poly, variable: sp.Symbol, parameters: tuple[sp.Symbol, ...]
) -> RootClassificationResult:
    coeffs = poly.all_coeffs()
    leading = sp.factor(coeffs[0])
    constant = sp.factor(coeffs[1])
    conditions = (
        (sp.Ne(leading, 0), sp.Integer(1), (1,)),
        (sp.And(sp.Eq(leading, 0), sp.Eq(constant, 0)), sp.oo, ()),
        (sp.And(sp.Eq(leading, 0), sp.Ne(constant, 0)), sp.Integer(0), ()),
    )
    cells = tuple(
        RootClassificationCell(
            condition,
            count,
            pattern,
            _sample_for(condition, parameters),
            {"degree_case": "linear"},
        )
        for condition, count, pattern in conditions
    )
    return RootClassificationResult(
        sp.expand(poly.as_expr()),
        variable,
        parameters,
        cells,
        "linear_parameter_classification",
        {"leading_coefficient": leading, "constant_coefficient": constant},
    )


def _quadratic_family(
    poly: sp.Poly, variable: sp.Symbol, parameters: tuple[sp.Symbol, ...]
) -> RootClassificationResult:
    coeffs = poly.all_coeffs()
    leading, middle, constant = (sp.factor(c) for c in coeffs)
    disc = sp.factor(middle**2 - 4 * leading * constant)
    # The nondegenerate quadratic strata are first because they are the generic
    # cases used most often by callers. Degenerate linear strata keep the API
    # mathematically complete for parameter values with a zero leading term.
    conditions = [
        (
            sp.And(sp.Ne(leading, 0), sp.Gt(disc, 0)),
            sp.Integer(2),
            (1, 1),
            "quadratic_discriminant_positive",
        ),
        (
            sp.And(sp.Ne(leading, 0), sp.Eq(disc, 0)),
            sp.Integer(1),
            (2,),
            "quadratic_discriminant_zero",
        ),
        (
            sp.And(sp.Ne(leading, 0), sp.Lt(disc, 0)),
            sp.Integer(0),
            (),
            "quadratic_discriminant_negative",
        ),
        (sp.And(sp.Eq(leading, 0), sp.Ne(middle, 0)), sp.Integer(1), (1,), "linear_degenerate"),
        (
            sp.And(sp.Eq(leading, 0), sp.Eq(middle, 0), sp.Eq(constant, 0)),
            sp.oo,
            (),
            "zero_polynomial",
        ),
        (
            sp.And(sp.Eq(leading, 0), sp.Eq(middle, 0), sp.Ne(constant, 0)),
            sp.Integer(0),
            (),
            "nonzero_constant",
        ),
    ]
    # For monic quadratic families, keep the printed public conditions concise.
    if sp.simplify(leading - 1) == 0:
        conditions = [
            (sp.Gt(disc, 0), sp.Integer(2), (1, 1), "quadratic_discriminant_positive"),
            (sp.Eq(disc, 0), sp.Integer(1), (2,), "quadratic_discriminant_zero"),
            (sp.Lt(disc, 0), sp.Integer(0), (), "quadratic_discriminant_negative"),
        ]
    cells = tuple(
        RootClassificationCell(
            condition,
            count,
            pattern,
            _sample_for(condition, parameters),
            {"case": case, "discriminant": disc},
        )
        for condition, count, pattern, case in conditions
    )
    return RootClassificationResult(
        sp.expand(poly.as_expr()),
        variable,
        parameters,
        cells,
        "quadratic_discriminant_classification",
        {"degree": 2, "discriminant": disc, "leading_coefficient": leading},
    )


def _sampled_parameter_cells(
    poly: sp.Poly,
    variable: sp.Symbol,
    parameters: tuple[sp.Symbol, ...],
) -> RootClassificationResult:
    """Fallback classification by discriminant sign cells for unsupported degrees.

    This is intentionally modest: it gives callers useful parameter conditions
    around the discriminant without claiming a complete real-root-classification
    algorithm for arbitrary high-degree families.
    """

    disc = sp.factor(sp.discriminant(poly.as_expr(), variable))
    cells: list[RootClassificationCell] = []
    for condition in (sp.Gt(disc, 0), sp.Eq(disc, 0), sp.Lt(disc, 0)):
        sample = _sample_for(condition, parameters)
        count: sp.Expr = sp.Integer(-1)
        pattern: tuple[int, ...] = ()
        if sample is not None:
            specialized = sp.Poly(poly.as_expr().subs(sample), variable)
            try:
                pattern = _multiplicity_pattern(specialized)
                count = sp.Integer(len(pattern))
            except Exception:
                count = sp.Integer(-1)
        cells.append(
            RootClassificationCell(
                condition,
                count,
                pattern,
                sample,
                {"case": "sampled_discriminant_cell", "complete": False, "discriminant": disc},
            )
        )
    return RootClassificationResult(
        sp.expand(poly.as_expr()),
        variable,
        parameters,
        tuple(cells),
        "sampled_discriminant_classification",
        {"complete": False, "degree": poly.degree(), "discriminant": disc},
    )


def classify_real_roots(
    polynomial: sp.Poly | sp.Expr,
    variable: sp.Symbol | str,
    *,
    parameters: Sequence[sp.Symbol | str] | None = None,
) -> RootClassificationResult:
    """Classify real roots of a univariate polynomial or polynomial family.

    The current public implementation is exact for unparameterized polynomials,
    complete for linear and quadratic parameter families, and provides a marked
    sampled discriminant stratification for higher-degree parameter families.
    """

    var = sp.Symbol(variable, real=True) if isinstance(variable, str) else variable
    expr = polynomial.as_expr() if isinstance(polynomial, sp.Poly) else sp.sympify(polynomial)
    params = _normalize_parameters(parameters, expr, var)
    poly = sp.Poly(expr, var)
    if not params:
        return _unparameterized(expr, var)
    degree = poly.degree()
    if degree == 0:
        condition = sp.Eq(poly.as_expr(), 0)
        nonzero = sp.Ne(poly.as_expr(), 0)
        cells = (
            RootClassificationCell(
                condition,
                sp.oo,
                (),
                _sample_for(condition, params),
                {"degree_case": "zero_constant"},
            ),
            RootClassificationCell(
                nonzero,
                sp.Integer(0),
                (),
                _sample_for(nonzero, params),
                {"degree_case": "nonzero_constant"},
            ),
        )
        return RootClassificationResult(
            sp.expand(expr), var, params, cells, "constant_parameter_classification", {"degree": 0}
        )
    if degree == 1:
        return _linear_family(poly, var, params)
    if degree == 2:
        return _quadratic_family(poly, var, params)
    return _sampled_parameter_cells(poly, var, params)


__all__ = ["RootClassificationCell", "RootClassificationResult", "classify_real_roots"]
