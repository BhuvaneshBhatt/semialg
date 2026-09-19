"""Exact convexity and monotonicity classification for semialgebraic functions."""

from __future__ import annotations

from collections.abc import Sequence

import sympy as sp

from .conditional import ParameterStratifiedResult, conditional_result
from .decision import is_satisfiable
from .decomposition import cad
from .derived_geometry import is_connected
from .domain_solve import function_domain
from .formula import parse_formula
from .function_analysis import (
    FunctionPropertyPartitionResult,
    _algebraic_domain_formula,
    _domain_graph,
    _function_graph,
)
from .function_graph import (
    UnsupportedFunctionGraph,
)
from .normalization import normalize_formula, normalize_variables
from .qe import qe_by_complete_cad
from .reasoning_signs import function_sign
from .reconstruct.cylindrical import path_condition
from .sampling import sample_point, sign_at
from .structural_keys import symbol_identity_key


def _function_relation_region(
    expression: sp.Expr,
    domain: sp.Expr,
    retained_symbols: sp.Symbol | Sequence[sp.Symbol],
    relation: str,
) -> sp.Expr | None:
    """Project an exact sign relation for a supported semialgebraic function."""

    try:
        domain_formula, domain_aux = _domain_graph(domain)
        graph_formula, value, graph_aux = _function_graph(
            expression, f"semialg_partition_{relation}_value"
        )
    except UnsupportedFunctionGraph:
        return None
    comparison = {
        "positive": value > 0,
        "negative": value < 0,
        "zero": sp.Eq(value, 0),
    }[relation]
    retained = (
        (retained_symbols,) if isinstance(retained_symbols, sp.Symbol) else tuple(retained_symbols)
    )
    auxiliaries = tuple(dict.fromkeys((value, *domain_aux, *graph_aux)))
    if not auxiliaries:
        return normalize_formula(sp.And(domain_formula, graph_formula, comparison))
    result = qe_by_complete_cad(
        (*retained, *auxiliaries),
        tuple(("exists", auxiliary) for auxiliary in auxiliaries),
        parse_formula(sp.And(domain_formula, graph_formula, comparison)),
        free_variables=retained,
        return_result=True,
    )
    return normalize_formula(result.formula)


def _ordered_univariate_regions(
    regions: Sequence[sp.Expr], variable: sp.Symbol
) -> tuple[sp.Expr, ...]:
    """Return deterministic left-to-right ordering using exact sample points."""

    keyed = []
    for region in regions:
        try:
            point = sample_point(region, (variable,))
            value = point[variable]
        except (TypeError, ValueError, NotImplementedError, sp.PolynomialError, KeyError):
            value = sp.S.NaN
        keyed.append((sp.default_sort_key(value), region))
    return tuple(region for _, region in sorted(keyed, key=lambda item: item[0]))


def _merge_property_partition(
    pieces: Sequence[tuple[str, sp.Expr]],
    variable: sp.Symbol,
    classifier,
) -> tuple[tuple[str, sp.Expr], ...]:
    """Merge adjacent CAD pieces whenever the public property query certifies the union."""

    current = list(pieces)
    changed = True
    while changed and len(current) > 1:
        changed = False
        merged: list[tuple[str, sp.Expr]] = []
        index = 0
        while index < len(current):
            if index + 1 >= len(current):
                merged.append(current[index])
                break
            left_class, left_region = current[index]
            right_class, right_region = current[index + 1]
            union = normalize_formula(sp.Or(left_region, right_region))
            try:
                connected = is_connected(union, (variable,))
            except (TypeError, ValueError, NotImplementedError, sp.PolynomialError):
                connected = False
            if connected:
                try:
                    union_class = classifier(union)
                except (TypeError, ValueError, NotImplementedError, sp.PolynomialError):
                    union_class = "unknown"
                preferred = {left_class, right_class} - {"constant", "affine", "unknown"}
                if union_class == left_class == right_class or (
                    len(preferred) == 1 and union_class in preferred
                ):
                    merged.append((union_class, union))
                    index += 2
                    changed = True
                    continue
            merged.append(current[index])
            index += 1
        current = merged
    return tuple(current)


def _partition_from_derivative_sign(
    expression: sp.Expr,
    variable: sp.Symbol,
    explicit_domain: sp.Expr,
    *,
    derivative_order: int,
    property_name: str,
) -> FunctionPropertyPartitionResult:
    """Partition a univariate domain into cells on which the requested derivative sign is invariant."""
    natural_domain = function_domain(expression, (variable,))
    effective_domain = _algebraic_domain_formula(
        normalize_formula(sp.And(explicit_domain, natural_domain)), (variable,)
    )
    if not is_satisfiable(effective_domain, (variable,)):
        return FunctionPropertyPartitionResult(
            expression,
            variable,
            effective_domain,
            natural_domain,
            property_name,
            (),
            "empty_domain",
            {"derivative_order": derivative_order},
        )

    derivative = sp.diff(expression, variable, derivative_order)
    try:
        sign_cover = sp.Or(derivative < 0, sp.Eq(derivative, 0), derivative > 0)
        cad_result = cad(
            sp.And(effective_domain, sign_cover),
            (variable,),
            output="cells",
            return_result=True,
        )
    except (TypeError, ValueError, NotImplementedError, sp.PolynomialError):
        cad_result = None

    fallback_cells: list[tuple[int, object, object]] = []
    if cad_result is None or not cad_result.cells:
        for relation, sign_value in (("negative", -1), ("zero", 0), ("positive", 1)):
            region = _function_relation_region(derivative, effective_domain, variable, relation)
            if region is None or region is sp.false or region == sp.false:
                continue
            try:
                relation_cad = cad(region, (variable,), output="cells", return_result=True)
            except (TypeError, ValueError, NotImplementedError, sp.PolynomialError):
                continue
            for cell in relation_cad.cells:
                fallback_cells.append((sign_value, cell, relation_cad.cad.cells_by_level))
        if not fallback_cells:

            def classify_global(region):
                from .function_analysis import function_convexity, function_monotonicity

                if property_name == "monotonicity":
                    return function_monotonicity(expression, variable, domain=region)
                return function_convexity(expression, (variable,), domain=region)

            classification = classify_global(effective_domain)
            return FunctionPropertyPartitionResult(
                expression,
                variable,
                effective_domain,
                natural_domain,
                property_name,
                ((classification, effective_domain),),
                "global_property_fallback",
                {"derivative": derivative},
            )

    cell_pieces: list[tuple[str, sp.Expr, bool]] = []
    source_cells = (
        [(None, cell, cad_result.cad.cells_by_level) for cell in cad_result.cells]
        if cad_result is not None and cad_result.cells
        else fallback_cells
    )
    source_cells.sort(key=lambda item: item[1].index)
    for known_sign, cell, cells_by_level in source_cells:
        region = path_condition(cell, (variable,), cells_by_level, closed=False)
        sample_value = cell.sample_exprs[0]
        derivative_sign = known_sign
        if derivative_sign is None:
            try:
                derivative_sign = sign_at(
                    derivative, {variable: sample_value}, variables=(variable,)
                )
            except (TypeError, ValueError, NotImplementedError, sp.PolynomialError):
                value = sp.simplify(derivative.subs(variable, sample_value))
                derivative_sign = (
                    1
                    if value.is_positive
                    else -1
                    if value.is_negative
                    else 0
                    if value.is_zero
                    else None
                )
        if property_name == "monotonicity":
            classification = {
                1: "strictly_increasing",
                -1: "strictly_decreasing",
                0: "constant",
            }.get(derivative_sign, "unknown")
        else:
            classification = {1: "convex", -1: "concave", 0: "affine"}.get(
                derivative_sign, "unknown"
            )
        cell_pieces.append((classification, region, cell.is_section))

    # Merge sign-compatible neighboring sectors and absorb isolated zero
    # sections.  Isolated derivative zeros do not destroy strict monotonicity;
    # isolated second-derivative zeros do not split a convex/concave region.
    merged: list[tuple[str, sp.Expr, bool]] = []
    index = 0
    neutral = "constant" if property_name == "monotonicity" else "affine"
    while index < len(cell_pieces):
        if (
            index + 2 < len(cell_pieces)
            and cell_pieces[index + 1][2]
            and cell_pieces[index + 1][0] == neutral
            and cell_pieces[index][0] == cell_pieces[index + 2][0]
            and cell_pieces[index][0] not in {neutral, "unknown"}
        ):
            merged.append(
                (
                    cell_pieces[index][0],
                    normalize_formula(
                        sp.Or(
                            cell_pieces[index][1],
                            cell_pieces[index + 1][1],
                            cell_pieces[index + 2][1],
                        )
                    ),
                    False,
                )
            )
            index += 3
            continue
        current = cell_pieces[index]
        if (
            current[2]
            and current[0] == neutral
            and merged
            and merged[-1][0] not in {neutral, "unknown"}
        ):
            previous = merged.pop()
            merged.append((previous[0], normalize_formula(sp.Or(previous[1], current[1])), False))
            index += 1
            continue
        merged.append(current)
        index += 1

    pieces = tuple((classification, region) for classification, region, _ in merged)
    return FunctionPropertyPartitionResult(
        expression,
        variable,
        effective_domain,
        natural_domain,
        property_name,
        pieces,
        "univariate_derivative_sign_cad",
        {"derivative": derivative, "derivative_order": derivative_order},
    )


def _parameterized_derivative_partition(
    expression: sp.Expr,
    variable: sp.Symbol,
    parameters: tuple[sp.Symbol, ...],
    explicit_domain: sp.Expr,
    *,
    derivative_order: int,
    property_name: str,
) -> ParameterStratifiedResult:
    """Return a cylindrical parameter-stratified univariate property partition."""

    all_variables = (*parameters, variable)
    natural_domain = function_domain(expression, all_variables)
    effective_domain = _algebraic_domain_formula(
        normalize_formula(sp.And(explicit_domain, natural_domain)), all_variables
    )
    derivative = sp.diff(expression, variable, derivative_order)
    sign_cover = sp.Or(derivative < 0, sp.Eq(derivative, 0), derivative > 0)
    try:
        decomposition = cad(
            sp.And(effective_domain, sign_cover),
            all_variables,
            output="cells",
            return_result=True,
        )
    except (TypeError, ValueError, NotImplementedError, sp.PolynomialError):
        return conditional_result(
            parameters,
            ((sp.true, "unknown_partition"),),
            method=f"function_{property_name}_partition_parameter_unknown",
        )

    parameter_level = len(parameters)
    grouped: dict[tuple[int, ...], list[tuple[str, sp.Expr]]] = {}
    parameter_conditions: dict[tuple[int, ...], sp.Expr] = {}
    for cell in sorted(decomposition.cells, key=lambda item: item.index):
        prefix = cell.index[:parameter_level]
        if prefix not in parameter_conditions:
            parameter_cell = next(
                item
                for item in decomposition.cad.cells_by_level[parameter_level]
                if item.index == prefix
            )
            parameter_conditions[prefix] = path_condition(
                parameter_cell,
                parameters,
                decomposition.cad.cells_by_level,
                closed=False,
            )
        full_region = path_condition(
            cell, all_variables, decomposition.cad.cells_by_level, closed=False
        )
        point = dict(zip(all_variables, cell.sample_exprs, strict=True))
        try:
            derivative_sign = sign_at(derivative, point, variables=all_variables)
        except (TypeError, ValueError, NotImplementedError, sp.PolynomialError):
            value = sp.simplify(derivative.subs(point))
            derivative_sign = (
                1
                if value.is_positive
                else -1
                if value.is_negative
                else 0
                if value.is_zero
                else None
            )
        if property_name == "monotonicity":
            classification = {
                1: "strictly_increasing",
                -1: "strictly_decreasing",
                0: "constant",
            }.get(derivative_sign, "unknown")
        else:
            classification = {1: "convex", -1: "concave", 0: "affine"}.get(
                derivative_sign, "unknown"
            )
        grouped.setdefault(prefix, []).append((classification, full_region))

    branches = tuple(
        (parameter_conditions[prefix], tuple(grouped[prefix])) for prefix in sorted(grouped)
    )
    if not branches:
        branches = ((sp.true, ()),)
    return conditional_result(
        parameters,
        branches,
        method=f"function_{property_name}_partition_parameter_cad",
        diagnostics={"derivative": derivative, "derivative_order": derivative_order},
    )


def _function_sign_partition_without_parameters(
    expression: sp.Expr,
    variable: sp.Symbol,
    explicit_domain: sp.Expr,
) -> FunctionPropertyPartitionResult:
    """Partition a univariate function domain into connected exact sign regions."""

    natural_domain = function_domain(expression, (variable,))
    effective_domain = _algebraic_domain_formula(
        normalize_formula(sp.And(explicit_domain, natural_domain)), (variable,)
    )
    if not is_satisfiable(effective_domain, (variable,)):
        return FunctionPropertyPartitionResult(
            expression,
            variable,
            effective_domain,
            natural_domain,
            "sign",
            (),
            "empty_domain",
            {},
        )

    # Polynomial/rational functions need only one sign-invariant CAD. Keeping
    # the sign cover unevaluated preserves numerator, denominator, and root
    # boundaries in the same decomposition.
    try:
        direct_rational = bool(expression.is_rational_function(variable))
    except (TypeError, ValueError):
        direct_rational = False
    if direct_rational:
        try:
            numerator, denominator = sp.fraction(sp.cancel(expression))
            numerator_cover = sp.Or(
                numerator < 0, sp.Eq(numerator, 0), numerator > 0, evaluate=False
            )
            denominator_cover = sp.Or(
                denominator < 0, sp.Eq(denominator, 0), denominator > 0, evaluate=False
            )
            direct_cad = cad(
                sp.And(effective_domain, numerator_cover, denominator_cover, evaluate=False),
                (variable,),
                output="cells",
                return_result=True,
            )
            direct_pieces: list[tuple[str, sp.Expr]] = []
            for cell in sorted(direct_cad.cells, key=lambda item: item.index):
                region = path_condition(
                    cell, (variable,), direct_cad.cad.cells_by_level, closed=False
                )
                value = cell.sample_exprs[0]
                sign_value = sign_at(expression, {variable: value}, variables=(variable,))
                classification = {1: "positive", 0: "zero", -1: "negative"}[sign_value]
                direct_pieces.append((classification, region))
            return FunctionPropertyPartitionResult(
                expression,
                variable,
                effective_domain,
                natural_domain,
                "sign",
                tuple(direct_pieces),
                "univariate_rational_sign_cad",
                {"cad_cell_count": len(direct_pieces)},
            )
        except (TypeError, ValueError, NotImplementedError, sp.PolynomialError, KeyError):
            pass

    relation_regions: list[tuple[str, sp.Expr]] = []
    for relation in ("negative", "zero", "positive"):
        region = _function_relation_region(expression, effective_domain, variable, relation)
        if region is None or region == sp.false:
            continue
        try:
            region_cad = cad(region, (variable,), output="cells", return_result=True)
        except (TypeError, ValueError, NotImplementedError, sp.PolynomialError):
            continue
        cells = sorted(region_cad.cells, key=lambda cell: cell.index)
        for cell in cells:
            relation_regions.append(
                (
                    relation,
                    path_condition(cell, (variable,), region_cad.cad.cells_by_level, closed=False),
                )
            )

    if not relation_regions:
        classification = function_sign(expression, (variable,), assumptions=effective_domain)
        return FunctionPropertyPartitionResult(
            expression,
            variable,
            effective_domain,
            natural_domain,
            "sign",
            ((classification, effective_domain),),
            "global_sign_fallback",
            {},
        )

    ordered = list(
        _ordered_univariate_regions(tuple(region for _, region in relation_regions), variable)
    )
    by_region = {sp.srepr(region): classification for classification, region in relation_regions}
    pieces = tuple((by_region[sp.srepr(region)], region) for region in ordered)

    # CADs for a single sign relation can contain artificial adjacent cells.
    # Merge only when their union is connected and the exact sign query proves
    # the same classification on that union.  Zero boundaries therefore remain
    # explicit and positive/negative components are never joined across them.
    pieces = _merge_property_partition(
        pieces,
        variable,
        lambda region: function_sign(expression, (variable,), assumptions=region),
    )
    return FunctionPropertyPartitionResult(
        expression,
        variable,
        effective_domain,
        natural_domain,
        "sign",
        pieces,
        "univariate_function_graph_sign_cad",
        {},
    )


def _function_sign_partition_with_parameters(
    expression: sp.Expr,
    variable: sp.Symbol,
    parameters: tuple[sp.Symbol, ...],
    explicit_domain: sp.Expr,
) -> ParameterStratifiedResult:
    """Return a parameter-first cylindrical exact sign partition."""

    retained = (*parameters, variable)
    natural_domain = function_domain(expression, retained)
    effective_domain = _algebraic_domain_formula(
        normalize_formula(sp.And(explicit_domain, natural_domain)), retained
    )
    relation_formulas = {
        relation: _function_relation_region(expression, effective_domain, retained, relation)
        for relation in ("negative", "zero", "positive")
    }
    if any(formula is None for formula in relation_formulas.values()):
        return conditional_result(
            parameters,
            ((sp.true, "unknown_partition"),),
            method="function_sign_partition_parameter_unknown",
        )

    # Build one CAD containing every polynomial boundary appearing in the three
    # projected sign formulas.  The conjunction with an unevaluated tautological
    # sign cover prevents simplification from discarding those boundaries.
    cover = sp.Or(*(relation_formulas.values()), evaluate=False)
    try:
        decomposition = cad(
            sp.And(effective_domain, cover, evaluate=False),
            retained,
            output="cells",
            return_result=True,
        )
    except (TypeError, ValueError, NotImplementedError, sp.PolynomialError):
        return conditional_result(
            parameters,
            ((sp.true, "unknown_partition"),),
            method="function_sign_partition_parameter_unknown",
        )

    parameter_level = len(parameters)
    grouped: dict[tuple[int, ...], list[tuple[str, sp.Expr]]] = {}
    parameter_conditions: dict[tuple[int, ...], sp.Expr] = {}
    for cell in sorted(decomposition.cells, key=lambda item: item.index):
        prefix = cell.index[:parameter_level]
        if prefix not in parameter_conditions:
            parameter_cell = next(
                item
                for item in decomposition.cad.cells_by_level[parameter_level]
                if item.index == prefix
            )
            parameter_conditions[prefix] = path_condition(
                parameter_cell, parameters, decomposition.cad.cells_by_level, closed=False
            )
        full_region = path_condition(cell, retained, decomposition.cad.cells_by_level, closed=False)
        classification = "unknown"
        for relation in ("negative", "zero", "positive"):
            try:
                if is_satisfiable(sp.And(full_region, relation_formulas[relation]), retained):
                    classification = relation
                    break
            except (TypeError, ValueError, NotImplementedError, sp.PolynomialError):
                pass
        grouped.setdefault(prefix, []).append((classification, full_region))

    branches = tuple(
        (parameter_conditions[prefix], tuple(grouped[prefix])) for prefix in sorted(grouped)
    )
    return conditional_result(
        parameters,
        branches or ((sp.true, ()),),
        method="function_sign_partition_parameter_cad",
        diagnostics={"relations": relation_formulas, "domain": effective_domain},
    )


def function_sign_partition(
    expression,
    variable: sp.Symbol | str | None = None,
    *,
    domain=sp.true,
    parameters: Sequence[sp.Symbol | str] | None = None,
    return_result: bool = False,
):
    """Partition a univariate semialgebraic function into exact sign regions.

    Each connected piece is classified as ``'positive'``, ``'negative'``, or
    ``'zero'``.  The natural real function domain is intersected automatically
    with ``domain``.  Other free symbols are treated as parameters when the
    function variable is explicit, yielding a :class:`ParameterStratifiedResult`.
    """

    expression = sp.sympify(expression)
    domain = normalize_formula(domain)
    combined = sp.Tuple(expression, domain)
    if variable is None:
        explicit_parameters = (
            normalize_variables(parameters, combined, append_context_symbols=False)
            if parameters is not None
            else ()
        )
        candidates = tuple(
            sorted(combined.free_symbols - set(explicit_parameters), key=symbol_identity_key)
        )
        if len(candidates) != 1:
            raise ValueError(
                "variable must be supplied unless exactly one non-parameter symbol is present"
            )
        variable = candidates[0]
    else:
        variable = normalize_variables((variable,), combined, append_context_symbols=False)[0]
    if parameters is None:
        parameters = tuple(sorted(combined.free_symbols - {variable}, key=symbol_identity_key))
    else:
        parameters = normalize_variables(parameters, combined, append_context_symbols=False)
    if parameters:
        return _function_sign_partition_with_parameters(
            expression, variable, tuple(parameters), domain
        )
    result = _function_sign_partition_without_parameters(expression, variable, domain)
    return result if return_result else result.pieces


def function_monotonic_partition(
    expression,
    variable: sp.Symbol | str | None = None,
    *,
    domain=sp.true,
    parameters: Sequence[sp.Symbol | str] | None = None,
    return_result: bool = False,
):
    """Partition a univariate semialgebraic domain into monotonic regions.

    Candidate regions come from the exact sign decomposition of the first
    derivative (including graph algebraization for supported nonsmooth
    derivatives), and every region is re-certified with
    :func:`function_monotonicity`.  Adjacent regions are merged whenever the
    combined region has the same strongest exact monotonicity property.
    """

    expression = sp.sympify(expression)
    domain = normalize_formula(domain)
    combined = sp.Tuple(expression, domain)
    if variable is None:
        explicit_parameters = (
            normalize_variables(parameters, combined, append_context_symbols=False)
            if parameters is not None
            else ()
        )
        candidates = tuple(
            sorted(combined.free_symbols - set(explicit_parameters), key=symbol_identity_key)
        )
        if len(candidates) != 1:
            raise ValueError(
                "variable must be supplied unless exactly one non-parameter symbol is present"
            )
        variable = candidates[0]
    else:
        variable = normalize_variables((variable,), combined, append_context_symbols=False)[0]
    if parameters is None:
        parameters = tuple(sorted(combined.free_symbols - {variable}, key=symbol_identity_key))
    else:
        parameters = normalize_variables(parameters, combined, append_context_symbols=False)
    if parameters:
        return _parameterized_derivative_partition(
            expression,
            variable,
            tuple(parameters),
            domain,
            derivative_order=1,
            property_name="monotonicity",
        )
    result = _partition_from_derivative_sign(
        expression, variable, domain, derivative_order=1, property_name="monotonicity"
    )
    return result if return_result else result.pieces


def function_convex_partition(
    expression,
    variable: sp.Symbol | str | None = None,
    *,
    domain=sp.true,
    parameters: Sequence[sp.Symbol | str] | None = None,
    return_result: bool = False,
):
    """Partition a univariate semialgebraic domain into convexity regions.

    Candidate regions come from the exact sign decomposition of the second
    derivative and are re-certified with :func:`function_convexity`, so
    isolated inflection/zero-curvature points do not force spurious splits.
    """

    expression = sp.sympify(expression)
    domain = normalize_formula(domain)
    combined = sp.Tuple(expression, domain)
    if variable is None:
        explicit_parameters = (
            normalize_variables(parameters, combined, append_context_symbols=False)
            if parameters is not None
            else ()
        )
        candidates = tuple(
            sorted(combined.free_symbols - set(explicit_parameters), key=symbol_identity_key)
        )
        if len(candidates) != 1:
            raise ValueError(
                "variable must be supplied unless exactly one non-parameter symbol is present"
            )
        variable = candidates[0]
    else:
        variable = normalize_variables((variable,), combined, append_context_symbols=False)[0]
    if parameters is None:
        parameters = tuple(sorted(combined.free_symbols - {variable}, key=symbol_identity_key))
    else:
        parameters = normalize_variables(parameters, combined, append_context_symbols=False)
    if parameters:
        return _parameterized_derivative_partition(
            expression,
            variable,
            tuple(parameters),
            domain,
            derivative_order=2,
            property_name="convexity",
        )
    result = _partition_from_derivative_sign(
        expression, variable, domain, derivative_order=2, property_name="convexity"
    )
    return result if return_result else result.pieces
