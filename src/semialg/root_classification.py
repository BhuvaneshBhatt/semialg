from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field

import sympy as sp

from .algebraic.samples import sample_to_expr
from .algebraic.signs import sign_at_sample
from .cad_algorithms.decomposition import decomp_collins_complete
from .normalization import normalize_parameters
from .sampling import sample_point
from .structural_keys import symbol_identity_key
from .symbol_resolution import resolve_symbol

_RECOVERABLE_ERRORS = (
    ArithmeticError,
    TypeError,
    ValueError,
    NotImplementedError,
    sp.PolynomialError,
)


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
        return normalize_parameters(parameters, expr, exclude=(variable,))
    return tuple(sorted(expr.free_symbols - {variable}, key=symbol_identity_key))


def _multiplicity_pattern(poly: sp.Poly) -> tuple[int, ...]:
    """Return real-root multiplicities using exact root isolation only."""

    try:
        roots = list(sp.real_roots(poly.as_expr()))
    except (NotImplementedError, sp.PolynomialError, ValueError):
        return ()
    if not roots:
        return ()
    groups: list[tuple[sp.Expr, int]] = []
    for root in roots:
        for index, (known, multiplicity) in enumerate(groups):
            if sp.simplify(root - known) == 0:
                groups[index] = (known, multiplicity + 1)
                break
        else:
            groups.append((root, 1))
    return tuple(sorted(multiplicity for _, multiplicity in groups))


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
    poly: sp.Poly,
    variable: sp.Symbol,
    parameters: tuple[sp.Symbol, ...],
    *,
    sample_cells: bool = True,
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
            _sample_for(condition, parameters) if sample_cells else None,
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
    poly: sp.Poly,
    variable: sp.Symbol,
    parameters: tuple[sp.Symbol, ...],
    *,
    sample_cells: bool = True,
) -> RootClassificationResult:
    """Classify a quadratic parameter family into exact real-root strata."""
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
            _sample_for(condition, parameters) if sample_cells else None,
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


def _classify_low_degree(
    poly: sp.Poly,
    variable: sp.Symbol,
    parameters: tuple[sp.Symbol, ...],
    *,
    sample_cells: bool = False,
) -> RootClassificationResult:
    """Dispatch a structurally lower-degree parameter family without wasted samples."""

    degree = poly.degree()
    if degree == 0:
        value = sp.factor(poly.as_expr())
        zero = sp.Eq(value, 0)
        nonzero = sp.Ne(value, 0)
        cells = (
            RootClassificationCell(
                zero,
                sp.oo,
                (),
                _sample_for(zero, parameters) if sample_cells else None,
                {"degree_case": "zero_constant"},
            ),
            RootClassificationCell(
                nonzero,
                sp.Integer(0),
                (),
                _sample_for(nonzero, parameters) if sample_cells else None,
                {"degree_case": "nonzero_constant"},
            ),
        )
        return RootClassificationResult(
            sp.expand(value),
            variable,
            parameters,
            cells,
            "constant_parameter_classification",
            {"degree": 0},
        )
    if degree == 1:
        return _linear_family(poly, variable, parameters, sample_cells=sample_cells)
    if degree == 2:
        return _quadratic_family(poly, variable, parameters, sample_cells=sample_cells)
    if degree == 3:
        return _cubic_family(poly, variable, parameters)
    raise ValueError("low-degree classifier received a polynomial of degree greater than three")


def _guard_cells(
    cells: tuple[RootClassificationCell, ...],
    guard: sp.Expr,
    parameters: tuple[sp.Symbol, ...],
    *,
    case_prefix: str,
) -> tuple[RootClassificationCell, ...]:
    """Restrict classification cells to an exact parameter guard."""

    guarded: list[RootClassificationCell] = []
    for cell in cells:
        condition = sp.And(guard, cell.condition)
        diagnostics = dict(cell.diagnostics)
        diagnostics["degree_guard"] = guard
        diagnostics["case"] = (
            f"{case_prefix}:{diagnostics.get('case', diagnostics.get('degree_case', 'cell'))}"
        )
        guarded.append(
            RootClassificationCell(
                condition,
                cell.root_count,
                cell.multiplicity_pattern,
                None,
                diagnostics,
            )
        )
    return tuple(guarded)


def _cubic_family(
    poly: sp.Poly, variable: sp.Symbol, parameters: tuple[sp.Symbol, ...]
) -> RootClassificationResult:
    """Classify a cubic family exactly, including leading-coefficient degree drops."""

    leading, quad, linear, constant = (sp.factor(c) for c in poly.all_coeffs())
    disc = sp.factor(sp.discriminant(poly.as_expr(), variable))
    triple_test = sp.factor(quad**2 - 3 * leading * linear)
    nondegenerate = sp.Ne(leading, 0)

    conditions = [
        (
            sp.And(nondegenerate, sp.Gt(disc, 0)),
            sp.Integer(3),
            (1, 1, 1),
            "cubic_discriminant_positive",
        ),
        (
            sp.And(nondegenerate, sp.Eq(disc, 0), sp.Eq(triple_test, 0)),
            sp.Integer(1),
            (3,),
            "cubic_triple_root",
        ),
        (
            sp.And(nondegenerate, sp.Eq(disc, 0), sp.Ne(triple_test, 0)),
            sp.Integer(2),
            (1, 2),
            "cubic_double_and_simple_root",
        ),
        (
            sp.And(nondegenerate, sp.Lt(disc, 0)),
            sp.Integer(1),
            (1,),
            "cubic_discriminant_negative",
        ),
    ]
    cells = [
        RootClassificationCell(
            condition,
            count,
            pattern,
            None,
            {"case": case, "discriminant": disc},
        )
        for condition, count, pattern, case in conditions
    ]

    degenerate_expr = sp.expand(quad * variable**2 + linear * variable + constant)
    degenerate = _classify_low_degree(
        sp.Poly(degenerate_expr, variable),
        variable,
        parameters,
        sample_cells=False,
    )
    cells.extend(
        _guard_cells(
            degenerate.cells,
            sp.Eq(leading, 0),
            parameters,
            case_prefix="cubic_degree_drop",
        )
    )
    return RootClassificationResult(
        sp.expand(poly.as_expr()),
        variable,
        parameters,
        tuple(cells),
        "cubic_discriminant_classification",
        {
            "degree": 3,
            "discriminant": disc,
            "leading_coefficient": leading,
            "triple_root_test": triple_test,
        },
    )


def _quartic_family(
    poly: sp.Poly, variable: sp.Symbol, parameters: tuple[sp.Symbol, ...]
) -> RootClassificationResult:
    """Classify nondegenerate quartic strata exactly and preserve unknown multiple-root strata.

    For a real quartic with nonzero leading coefficient and nonzero
    discriminant, the discriminant together with the classical ``P`` and ``D``
    invariants determines the number of distinct real roots:

    - discriminant < 0: exactly two real roots;
    - discriminant > 0, P < 0, D < 0: four real roots;
    - discriminant > 0 otherwise: no real roots.

    The discriminant-zero locus contains several multiple-root configurations
    and is deliberately left as an unknown-count stratum here. Leading-
    coefficient degree drops recurse through the complete cubic classifier.
    """

    leading, cubic, quad, linear, constant = (sp.factor(c) for c in poly.all_coeffs())
    expr = sp.expand(poly.as_expr())
    disc = sp.factor(sp.discriminant(expr, variable))
    p_inv = sp.factor(8 * leading * quad - 3 * cubic**2)
    d_inv = sp.factor(
        64 * leading**3 * constant
        - 16 * leading**2 * quad**2
        + 16 * leading * cubic**2 * quad
        - 16 * leading**2 * cubic * linear
        - 3 * cubic**4
    )
    nondegenerate = sp.Ne(leading, 0)

    conditions = [
        (
            sp.And(nondegenerate, sp.Lt(disc, 0)),
            sp.Integer(2),
            (1, 1),
            "quartic_discriminant_negative",
        ),
        (
            sp.And(nondegenerate, sp.Gt(disc, 0), sp.Lt(p_inv, 0), sp.Lt(d_inv, 0)),
            sp.Integer(4),
            (1, 1, 1, 1),
            "quartic_four_real_roots",
        ),
        (
            sp.And(
                nondegenerate,
                sp.Gt(disc, 0),
                sp.Or(sp.Ge(p_inv, 0), sp.Ge(d_inv, 0)),
            ),
            sp.Integer(0),
            (),
            "quartic_no_real_roots",
        ),
        (
            sp.And(nondegenerate, sp.Eq(disc, 0)),
            sp.Integer(-1),
            (),
            "quartic_multiple_root_unknown",
        ),
    ]
    cells = [
        RootClassificationCell(
            condition,
            count,
            pattern,
            None,
            {
                "case": case,
                "discriminant": disc,
                "P": p_inv,
                "D": d_inv,
                "complete": count != -1,
            },
        )
        for condition, count, pattern, case in conditions
    ]

    degenerate_expr = sp.expand(
        cubic * variable**3 + quad * variable**2 + linear * variable + constant
    )
    degenerate = _classify_low_degree(
        sp.Poly(degenerate_expr, variable),
        variable,
        parameters,
        sample_cells=False,
    )
    cells.extend(
        _guard_cells(
            degenerate.cells,
            sp.Eq(leading, 0),
            parameters,
            case_prefix="quartic_degree_drop",
        )
    )
    complete = all(cell.root_count != -1 for cell in cells)
    return RootClassificationResult(
        expr,
        variable,
        parameters,
        tuple(cells),
        "quartic_invariant_classification",
        {
            "complete": complete,
            "degree": 4,
            "discriminant": disc,
            "leading_coefficient": leading,
            "P": p_inv,
            "D": d_inv,
        },
    )


def _sturm_signature_polynomials(
    poly: sp.Poly, variable: sp.Symbol, parameters: tuple[sp.Symbol, ...]
) -> tuple[sp.Expr, ...]:
    """Return parameter polynomials whose signs determine the specialized Sturm profile.

    We use coefficients of the subresultant PRS of ``p`` and ``p'`` rather
    than a quotient-based Sturm chain.  Subresultants specialize polynomially
    in the parameters and encode degree drops/gcd changes without introducing
    parameter denominators.  A sign-invariant CAD for these coefficients is
    therefore a certified stratification on which the distinct real-root count
    is constant.
    """

    expr = sp.expand(poly.as_expr())
    derivative = sp.diff(expr, variable)
    try:
        chain = sp.subresultants(expr, derivative, variable)
    except _RECOVERABLE_ERRORS:
        chain = [expr, derivative]
    candidates: list[sp.Expr] = []
    for item in chain:
        try:
            item_poly = sp.Poly(sp.expand(item), variable)
        except _RECOVERABLE_ERRORS:
            continue
        candidates.extend(sp.factor(coeff) for coeff in item_poly.all_coeffs())
    # Input coefficients are retained explicitly even if SymPy shortens a
    # degenerate subresultant chain.  Constants do not stratify parameter space.
    candidates.extend(sp.factor(coeff) for coeff in poly.all_coeffs())
    unique: list[sp.Expr] = []
    seen: set[sp.Expr] = set()
    parameter_set = set(parameters)
    for candidate in candidates:
        candidate = sp.factor(candidate)
        if candidate == 0 or not candidate.free_symbols:
            continue
        if not candidate.free_symbols <= parameter_set:
            continue
        if candidate not in seen:
            seen.add(candidate)
            unique.append(candidate)
    return tuple(unique)


def _signature_condition(
    signature_polynomials: tuple[sp.Expr, ...], signs: tuple[int, ...]
) -> sp.Expr:
    atoms: list[sp.Expr] = []
    for polynomial, sign in zip(signature_polynomials, signs, strict=True):
        if sign < 0:
            atoms.append(polynomial < 0)
        elif sign > 0:
            atoms.append(polynomial > 0)
        else:
            atoms.append(sp.Eq(polynomial, 0))
    return sp.And(*atoms) if atoms else sp.true


def _exact_specialized_root_data(
    expr: sp.Expr,
    variable: sp.Symbol,
    parameters: tuple[sp.Symbol, ...],
    sample: tuple[object, ...],
) -> tuple[sp.Expr, tuple[int, ...], Mapping[sp.Symbol, sp.Expr]]:
    assignment = {
        parameter: sample_to_expr(value)
        for parameter, value in zip(parameters, sample, strict=True)
    }
    specialized = sp.expand(expr.subs(assignment))
    try:
        specialized_poly = sp.Poly(specialized, variable, extension=True)
    except _RECOVERABLE_ERRORS:
        specialized_poly = sp.Poly(specialized, variable)
    if specialized_poly.is_zero:
        return sp.oo, (), assignment
    if specialized_poly.degree() <= 0:
        return sp.Integer(0), (), assignment
    try:
        count = sp.Integer(specialized_poly.count_roots(-sp.oo, sp.oo))
    except _RECOVERABLE_ERRORS:
        pattern = _multiplicity_pattern(specialized_poly)
        return sp.Integer(len(pattern)), pattern, assignment
    pattern = _multiplicity_pattern(specialized_poly)
    return count, pattern, assignment


def _sturm_parameter_stratification(
    poly: sp.Poly, variable: sp.Symbol, parameters: tuple[sp.Symbol, ...]
) -> RootClassificationResult | None:
    """Classify an arbitrary-degree family by subresultant/Sturm sign strata.

    The method constructs a sign-invariant CAD in parameter space for the
    coefficient data of the subresultant PRS.  Root counts are evaluated at one
    exact sample per CAD cell; because the entire Sturm specialization profile
    is sign invariant on the cell, the sampled distinct-root count is a
    cell-wide certificate, not a heuristic extrapolation.
    """

    signatures = _sturm_signature_polynomials(poly, variable, parameters)
    if not signatures:
        return None
    try:
        signature_polys = tuple(sp.Poly(item, *parameters) for item in signatures)
        cad = decomp_collins_complete(signature_polys, parameters)
    except _RECOVERABLE_ERRORS:
        return None

    by_signature: dict[
        tuple[int, ...], tuple[sp.Expr, tuple[int, ...], Mapping[sp.Symbol, sp.Expr]]
    ] = {}
    conflicts: set[tuple[int, ...]] = set()
    for cell in cad.cells:
        try:
            signs = tuple(sign_at_sample(item, cell.sample) for item in signature_polys)
            count, pattern, assignment = _exact_specialized_root_data(
                poly.as_expr(), variable, parameters, cell.sample
            )
        except _RECOVERABLE_ERRORS:
            continue
        previous = by_signature.get(signs)
        if previous is not None and previous[0] != count:
            conflicts.add(signs)
        else:
            by_signature.setdefault(signs, (count, pattern, assignment))

    if not by_signature:
        return None

    cells: list[RootClassificationCell] = []
    for signs, (count, pattern, assignment) in sorted(
        by_signature.items(), key=lambda item: item[0]
    ):
        certified_count = sp.Integer(-1) if signs in conflicts else count
        condition = _signature_condition(signatures, signs)
        cells.append(
            RootClassificationCell(
                condition,
                certified_count,
                pattern if certified_count != -1 else (),
                assignment,
                {
                    "case": "subresultant_sturm_signature",
                    "signature": signs,
                    "complete": certified_count != -1,
                },
            )
        )
    return RootClassificationResult(
        sp.expand(poly.as_expr()),
        variable,
        parameters,
        tuple(cells),
        "subresultant_sturm_parameter_stratification",
        {
            "complete": not conflicts,
            "degree": poly.degree(),
            "signature_polynomials": signatures,
            "parameter_cad_cells": len(cad.cells),
            "conflicting_signatures": tuple(sorted(conflicts)),
        },
    )


def _sampled_parameter_cells(
    poly: sp.Poly,
    variable: sp.Symbol,
    parameters: tuple[sp.Symbol, ...],
) -> RootClassificationResult:
    """Fallback classification by discriminant sign cells for unsupported degrees.

    This is intentionally conservative: discriminant cells are useful diagnostics,
    but for degree four and above the sign of the discriminant does not determine
    the number of real roots. Such cells therefore carry root_count=-1 rather
    than promoting one sampled fiber to a cell-wide exact claim.
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
            except _RECOVERABLE_ERRORS:
                pattern = ()
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
    complete for linear, quadratic, and cubic parameter families (including
    coefficient-induced degree drops). Quartics use compact classical invariant
    formulas where possible and a subresultant/Sturm sign stratification on
    unresolved multiple-root loci. Higher-degree families use the same exact
    subresultant/Sturm parameter-CAD construction, with conservative unknown
    fallback only if exact stratification cannot be constructed.
    """

    expr = polynomial.as_expr() if isinstance(polynomial, sp.Poly) else sp.sympify(polynomial)
    var = resolve_symbol(variable, context=(expr,))
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
    if degree == 3:
        return _cubic_family(poly, var, params)
    if degree == 4:
        return _quartic_family(poly, var, params)
    sturm = _sturm_parameter_stratification(poly, var, params)
    return sturm if sturm is not None else _sampled_parameter_cells(poly, var, params)


__all__ = ["RootClassificationCell", "RootClassificationResult", "classify_real_roots"]
