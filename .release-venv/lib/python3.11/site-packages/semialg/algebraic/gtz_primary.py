"""Full certified GTZ primary decomposition over ``QQ``.

This module implements the recursive positive-dimensional stages of the
Gianni--Trager--Zacharias strategy.  A maximal independent set ``U`` reduces
the current ideal to a zero-dimensional algebra over ``QQ(U)``.  That finite
algebra is decomposed exactly from its multiplication matrices and trace
form.  The localized primary factors are contracted to ``QQ[U,X]`` and the
lower-dimensional residual branch is recovered by an exact saturation split.

The certificate tree records every localization, local Artinian split,
contraction, saturation exponent, recursive residual branch, and the final
irredundant reconstruction.  Replay recomputes all algebraic identities and
does not rerun the GTZ search heuristics.
"""

from __future__ import annotations

from collections.abc import Sequence
from functools import lru_cache
from itertools import product

import sympy as sp

from ._gtz_primary_types import (  # noqa: E402, F401
    FractionFieldZeroDimensionalCertificate,
    FractionFieldZeroDimensionalResult,
    GTZContractedComponentCertificate,
    GTZExecutionPlan,
    GTZNodeCertificate,
    GTZPrimaryDecompositionCertificate,
    GTZPrimaryDecompositionResult,
    LocalPrimaryFactor,
)
from .equality_ideal import EqualityIdealContext
from .groebner_utils import compute_groebner_basis
from .gtz import (
    _colon_once,
    certified_independent_localization,
    contract_localized_ideal,
    saturation_stabilization,
)
from .gtz_zero_dim import (
    zero_dimensional_primary_decomposition,
)
from .hilbert import ideal_degree
from .ideal_ops import canonical_qq_basis_uncached, ideal_intersection_qq, intersection_all_qq


def _prefer_modular_qq(generators: Sequence[sp.Expr], variables: Sequence[sp.Symbol]) -> bool:
    max_degree = 0
    for generator in generators:
        try:
            poly = sp.Poly(generator, *variables, domain=sp.QQ)
        except (sp.PolynomialError, sp.polys.polyerrors.CoercionFailed):
            return False
        max_degree = max(max_degree, int(poly.total_degree()))
    complexity = len(tuple(generators)) * max(1, len(tuple(variables))) * max(1, max_degree)
    return complexity >= 8 or len(tuple(variables)) >= 3


@lru_cache(maxsize=512)
def _canonical_qq_basis_cached(
    generators: tuple[sp.Expr, ...], variables: tuple[sp.Symbol, ...]
) -> tuple[sp.Expr, ...]:
    use_modular = _prefer_modular_qq(generators, variables)
    return canonical_qq_basis_uncached(
        generators, variables, modular=None if use_modular else False
    )


def _canonical_qq_basis(
    generators: Sequence[sp.Expr], variables: Sequence[sp.Symbol]
) -> tuple[sp.Expr, ...]:
    vars_ = tuple(variables)
    expressions = []
    for generator in generators:
        expanded = sp.expand(generator)
        if expanded != 0:
            expressions.append(expanded)
    expressions = tuple(expressions)
    return _canonical_qq_basis_cached(expressions, vars_)


def _qq_ideal_equal(
    left: Sequence[sp.Expr], right: Sequence[sp.Expr], variables: Sequence[sp.Symbol]
) -> bool:
    return _canonical_qq_basis(left, variables) == _canonical_qq_basis(right, variables)


def _ideal_intersection(
    left: Sequence[sp.Expr], right: Sequence[sp.Expr], variables: Sequence[sp.Symbol]
) -> tuple[sp.Expr, ...]:
    return ideal_intersection_qq(left, right, variables, canonicalizer=_canonical_qq_basis)


def _intersection_all(
    ideals: Sequence[Sequence[sp.Expr]], variables: Sequence[sp.Symbol]
) -> tuple[sp.Expr, ...]:
    return intersection_all_qq(ideals, variables, canonicalizer=_canonical_qq_basis)


def _ff_domain(parameters: Sequence[sp.Symbol]):
    params = tuple(parameters)
    return sp.QQ.frac_field(*params) if params else sp.QQ


@lru_cache(maxsize=256)
def _canonical_ff_basis_cached(
    generators: tuple[sp.Expr, ...],
    variables: tuple[sp.Symbol, ...],
    parameters: tuple[sp.Symbol, ...],
) -> tuple[sp.Expr, ...]:
    if not generators:
        return tuple()
    domain = _ff_domain(parameters)
    gb = compute_groebner_basis(generators, variables, order="grevlex", domain=domain, modular=None)
    return tuple(sp.cancel(p.as_expr()) for p in gb.polys)


def _canonical_ff_basis(
    generators: Sequence[sp.Expr], variables: Sequence[sp.Symbol], parameters: Sequence[sp.Symbol]
) -> tuple[sp.Expr, ...]:
    vars_ = tuple(variables)
    params = tuple(parameters)
    expressions = []
    for generator in generators:
        reduced = sp.cancel(generator)
        if reduced != 0:
            expressions.append(reduced)
    expressions = tuple(expressions)
    return _canonical_ff_basis_cached(expressions, vars_, params)


def _ff_ideal_equal(
    left: Sequence[sp.Expr],
    right: Sequence[sp.Expr],
    variables: Sequence[sp.Symbol],
    parameters: Sequence[sp.Symbol],
) -> bool:
    return _canonical_ff_basis(left, variables, parameters) == _canonical_ff_basis(
        right, variables, parameters
    )


def _ff_saturate(
    generators: Sequence[sp.Expr],
    splitter: sp.Expr,
    variables: Sequence[sp.Symbol],
    parameters: Sequence[sp.Symbol],
) -> tuple[sp.Expr, ...]:
    vars_ = tuple(variables)
    domain = _ff_domain(parameters)
    u = sp.Dummy("gtz45_sat")
    equations = tuple(sp.cancel(g) for g in generators) + (1 - u * sp.cancel(splitter),)
    gb = sp.groebner(
        equations,
        u,
        *vars_,
        order="lex",
        domain=domain,
    )
    eliminated = tuple(
        sp.cancel(p.as_expr()) for p in gb.polys if u not in p.as_expr().free_symbols
    )
    return _canonical_ff_basis(eliminated, vars_, parameters)


def _standard_monomials(
    basis: sp.polys.polytools.GroebnerBasis, variables: Sequence[sp.Symbol]
) -> tuple[tuple[int, ...], ...] | None:
    vars_ = tuple(variables)
    if not vars_:
        return (tuple(),)
    leading = [tuple(int(e) for e in p.LM(order=basis.order).exponents) for p in basis.polys]
    bounds: list[int] = []
    for i in range(len(vars_)):
        pure = [exp[i] for exp in leading if exp[i] and sum(exp) == exp[i]]
        if not pure:
            return None
        bounds.append(min(pure))
    out = []
    for exponent in product(*(range(bound) for bound in bounds)):
        if not any(all(e >= g for e, g in zip(exponent, lm, strict=True)) for lm in leading):
            out.append(tuple(exponent))
    return tuple(out)


def _monomial_expr(exponent: Sequence[int], variables: Sequence[sp.Symbol]) -> sp.Expr:
    out = sp.Integer(1)
    for variable, power in zip(variables, exponent, strict=True):
        out *= variable ** int(power)
    return out


def _remainder_vector(
    expression: sp.Expr,
    basis: sp.polys.polytools.GroebnerBasis,
    standard: Sequence[tuple[int, ...]],
    variables: Sequence[sp.Symbol],
    domain,
) -> tuple[sp.Expr, ...]:
    _quotients, remainder = basis.reduce(sp.cancel(expression))
    poly = sp.Poly(sp.cancel(remainder), *variables, domain=domain)
    coeffs = {tuple(int(e) for e in mon): domain.to_sympy(coeff) for mon, coeff in poly.terms()}
    return tuple(sp.cancel(coeffs.get(tuple(mon), 0)) for mon in standard)


def _multiplication_matrix(
    expression: sp.Expr,
    basis: sp.polys.polytools.GroebnerBasis,
    standard: Sequence[tuple[int, ...]],
    variables: Sequence[sp.Symbol],
    domain,
) -> sp.Matrix:
    columns = []
    for monomial in standard:
        columns.append(
            _remainder_vector(
                expression * _monomial_expr(monomial, variables),
                basis,
                standard,
                variables,
                domain,
            )
        )
    return sp.Matrix.hstack(*(sp.Matrix(column) for column in columns))


def _trace_data(
    basis: sp.polys.polytools.GroebnerBasis,
    standard: Sequence[tuple[int, ...]],
    variables: Sequence[sp.Symbol],
    domain,
) -> tuple[sp.Matrix, tuple[sp.Expr, ...], tuple[sp.Matrix, ...]]:
    basis_exprs = tuple(_monomial_expr(mon, variables) for mon in standard)
    multiplication = tuple(
        _multiplication_matrix(expr, basis, standard, variables, domain) for expr in basis_exprs
    )
    size = len(standard)
    trace = sp.zeros(size)
    for i in range(size):
        for j in range(i, size):
            value = sp.cancel(sp.trace(multiplication[i] * multiplication[j]))
            trace[i, j] = value
            trace[j, i] = value
    return trace, basis_exprs, multiplication


def _kernel_generators(trace: sp.Matrix, basis_exprs: Sequence[sp.Expr]) -> tuple[sp.Expr, ...]:
    out = []
    for vector in trace.nullspace():
        expr = sp.cancel(sum(vector[i] * basis_exprs[i] for i in range(len(basis_exprs))))
        if expr != 0:
            out.append(expr)
    return tuple(out)


def _primitive_coefficients(count: int):
    if count == 0:
        yield tuple()
        return
    if count == 1:
        yield (sp.Integer(1),)
        return
    c = 0
    while True:
        yield tuple(sp.Integer(c) ** i for i in range(count))
        c += 1


def _compute_fraction_field_zero_dimensional(
    generators: Sequence[sp.Expr],
    variables: Sequence[sp.Symbol],
    parameters: Sequence[sp.Symbol],
) -> FractionFieldZeroDimensionalResult:
    vars_ = tuple(variables)
    params = tuple(parameters)
    domain = _ff_domain(params)
    source = _canonical_ff_basis(generators, vars_, params)
    gb = compute_groebner_basis(source, vars_, order="grevlex", domain=domain, modular=None)
    if len(gb.polys) == 1 and gb.polys[0].is_one:
        standard: tuple[tuple[int, ...], ...] = tuple()
        certificate = FractionFieldZeroDimensionalCertificate(
            params,
            vars_,
            source,
            standard,
            0,
            tuple(),
            tuple(),
            sp.Integer(0),
            sp.Integer(1),
            tuple(),
            tuple(),
        )
        return FractionFieldZeroDimensionalResult(params, vars_, tuple(), True, certificate)
    standard = _standard_monomials(gb, vars_)
    if standard is None:
        raise ValueError("localized ideal is not zero-dimensional")
    trace, basis_exprs, _basis_multiplication = _trace_data(gb, standard, vars_, domain)
    trace_rank = int(trace.rank())
    nilradical = _kernel_generators(trace, basis_exprs)

    z = sp.Symbol("_gtz45_primitive")
    primitive_coeffs = None
    primitive = None
    charpoly_expr = None
    factored: tuple[tuple[sp.Expr, int], ...] | None = None
    for coeffs in _primitive_coefficients(len(vars_)):
        candidate = sp.expand(sum(c * v for c, v in zip(coeffs, vars_, strict=True)))
        matrix = _multiplication_matrix(candidate, gb, standard, vars_, domain)
        charpoly = sp.Poly(matrix.charpoly(z).as_expr(), z, domain=domain)
        _unit, factors = charpoly.factor_list()
        normalized_factors = tuple(
            (sp.cancel(f.monic().as_expr()), int(multiplicity)) for f, multiplicity in factors
        )
        reduced_degree = sum(sp.Poly(f, z, domain=domain).degree() for f, _m in normalized_factors)
        if reduced_degree == trace_rank:
            primitive_coeffs = tuple(coeffs)
            primitive = candidate
            charpoly_expr = sp.cancel(charpoly.monic().as_expr())
            factored = normalized_factors
            break
    if primitive_coeffs is None or primitive is None or factored is None:
        raise ArithmeticError("failed to find a primitive element for the reduced quotient")

    radical_source = _canonical_ff_basis((*source, *nilradical), vars_, params)
    components: list[LocalPrimaryFactor] = []
    for index, (factor, multiplicity) in enumerate(factored):
        p_of_t = sp.cancel(factor.subs(z, primitive))
        others = [f.subs(z, primitive) for j, (f, _m) in enumerate(factored) if j != index]
        separator = sp.cancel(sp.prod(others)) if others else sp.Integer(1)
        component = source if separator == 1 else _ff_saturate(source, separator, vars_, params)
        radical = _canonical_ff_basis((*radical_source, p_of_t), vars_, params)
        components.append(LocalPrimaryFactor(component, radical, factor, multiplicity, separator))

    certificate = FractionFieldZeroDimensionalCertificate(
        params,
        vars_,
        source,
        standard,
        trace_rank,
        nilradical,
        primitive_coeffs,
        primitive,
        charpoly_expr if charpoly_expr is not None else sp.Integer(1),
        factored,
        tuple(components),
    )
    result = FractionFieldZeroDimensionalResult(params, vars_, tuple(components), True, certificate)
    if not verify_fraction_field_zero_dimensional_certificate(certificate):
        raise ArithmeticError("internal fraction-field primary certificate verification failed")
    return result


def verify_fraction_field_zero_dimensional_certificate(
    certificate: FractionFieldZeroDimensionalCertificate,
) -> bool:
    try:
        vars_ = tuple(certificate.variables)
        params = tuple(certificate.parameters)
        domain = _ff_domain(params)
        source = _canonical_ff_basis(certificate.source_generators, vars_, params)
        if source != tuple(certificate.source_generators):
            return False
        gb = compute_groebner_basis(source, vars_, order="grevlex", domain=domain, modular=None)
        standard = _standard_monomials(gb, vars_)
        if standard is None or tuple(standard) != tuple(certificate.standard_monomials):
            return False
        trace, basis_exprs, _multiplication = _trace_data(gb, standard, vars_, domain)
        if int(trace.rank()) != certificate.trace_rank:
            return False
        nilradical = _kernel_generators(trace, basis_exprs)
        if _canonical_ff_basis(nilradical, vars_, params) != _canonical_ff_basis(
            certificate.nilradical_generators, vars_, params
        ):
            return False
        if len(certificate.primitive_coefficients) != len(vars_):
            return False
        primitive = sp.expand(
            sum(c * v for c, v in zip(certificate.primitive_coefficients, vars_, strict=True))
        )
        if sp.expand(primitive - certificate.primitive_element) != 0:
            return False
        characteristic_symbols = (
            sp.sympify(certificate.characteristic_polynomial).free_symbols
            - set(vars_)
            - set(params)
        )
        z = next(iter(characteristic_symbols), None)
        # Stored factors use a private dummy that survives in the expressions.
        if z is None and certificate.factors:
            z = next(iter(sp.sympify(certificate.factors[0][0]).free_symbols - set(params)), None)
        if z is None:
            z = sp.Dummy("gtz45_replay")
        matrix = _multiplication_matrix(primitive, gb, standard, vars_, domain)
        charpoly = sp.Poly(matrix.charpoly(z).as_expr(), z, domain=domain).monic()
        stored_charpoly = sp.Poly(certificate.characteristic_polynomial, z, domain=domain).monic()
        # Dummy identity can differ after serialization/replacement; compare coefficients.
        if tuple(charpoly.all_coeffs()) != tuple(stored_charpoly.all_coeffs()):
            return False
        _unit, factors = charpoly.factor_list()
        normalized = tuple((sp.cancel(f.monic().as_expr()), int(m)) for f, m in factors)
        stored = tuple(certificate.factors)
        if len(normalized) != len(stored):
            return False
        for (left, lm), (right, rm) in zip(normalized, stored, strict=True):
            if lm != rm:
                return False
            if tuple(sp.Poly(left, z, domain=domain).all_coeffs()) != tuple(
                sp.Poly(right, z, domain=domain).all_coeffs()
            ):
                return False
        reduced_degree = sum(sp.Poly(f, z, domain=domain).degree() for f, _m in normalized)
        if reduced_degree != certificate.trace_rank:
            return False
        radical_source = _canonical_ff_basis(
            (*source, *certificate.nilradical_generators), vars_, params
        )
        if len(certificate.components) != len(normalized):
            return False
        for index, component in enumerate(certificate.components):
            factor, multiplicity = normalized[index]
            if component.multiplicity != multiplicity:
                return False
            if tuple(sp.Poly(component.irreducible_factor, z, domain=domain).all_coeffs()) != tuple(
                sp.Poly(factor, z, domain=domain).all_coeffs()
            ):
                return False
            p_of_t = sp.cancel(factor.subs(z, primitive))
            others = [f.subs(z, primitive) for j, (f, _m) in enumerate(normalized) if j != index]
            separator = sp.cancel(sp.prod(others)) if others else sp.Integer(1)
            if sp.cancel(component.separator - separator) != 0:
                return False
            expected = source if separator == 1 else _ff_saturate(source, separator, vars_, params)
            if not _ff_ideal_equal(component.generators, expected, vars_, params):
                return False
            expected_radical = _canonical_ff_basis((*radical_source, p_of_t), vars_, params)
            if not _ff_ideal_equal(component.radical, expected_radical, vars_, params):
                return False
            # The reduced local quotient must be exactly the irreducible field factor.
            rgb = sp.groebner(expected_radical, *vars_, order="grevlex", domain=domain)
            rstandard = _standard_monomials(rgb, vars_)
            if rstandard is None:
                return False
            if len(rstandard) != sp.Poly(factor, z, domain=domain).degree():
                return False
        return True
    except (
        ArithmeticError,
        ValueError,
        TypeError,
        ZeroDivisionError,
        sp.PolynomialError,
        sp.polys.polyerrors.CoercionFailed,
    ):
        return False


def plan_gtz_primary_decomposition(
    generators: Sequence[sp.Expr], variables: Sequence[sp.Symbol]
) -> GTZExecutionPlan:
    """Return a deterministic cost plan for the recursive GTZ driver.

    The plan only selects acceleration policy.  It never weakens exact
    verification or changes the mathematical result.
    """
    vars_ = tuple(variables)
    source = _canonical_qq_basis(generators, vars_)
    context = EqualityIdealContext.build(source, vars_)
    dimension = int(context.dimension)
    max_degree = 0
    for generator in source:
        try:
            poly = sp.Poly(generator, *vars_, domain=sp.QQ)
            max_degree = max(max_degree, int(poly.total_degree()))
        except (sp.PolynomialError, sp.polys.polyerrors.CoercionFailed):
            pass
    try:
        degree = int(ideal_degree(source, vars_))
    except (
        ArithmeticError,
        ValueError,
        TypeError,
        sp.PolynomialError,
        sp.polys.polyerrors.CoercionFailed,
    ):
        degree = None
    quotient_degree = degree if degree is not None and dimension == 0 else None
    use_modular = _prefer_modular_qq(source, vars_)
    reason = (
        "modular_qq_groebner_preferred"
        if use_modular
        else "small_exact_problem_modular_fallback_available"
    )
    return GTZExecutionPlan(
        vars_,
        len(source),
        dimension,
        degree,
        max_degree,
        quotient_degree,
        use_modular,
        True,
        "trace_rank_then_deterministic_linear_form",
        reason,
    )


def _parameter_annihilator(
    source: Sequence[sp.Expr],
    target_generator: sp.Expr,
    independent: Sequence[sp.Symbol],
    dependent: Sequence[sp.Symbol],
    variables: Sequence[sp.Symbol],
) -> sp.Expr | None:
    colon = _colon_once(source, target_generator, variables)
    if _qq_ideal_equal(colon, (1,), variables):
        return sp.Integer(1)
    if not independent:
        return None
    gb = compute_groebner_basis(
        colon,
        (*dependent, *independent),
        order="lex",
        domain=sp.QQ,
        modular=None,
    )
    independent_set = set(independent)
    candidates = [
        sp.expand(p.as_expr())
        for p in gb.polys
        if p.as_expr() != 0 and p.as_expr().free_symbols.issubset(independent_set)
    ]
    return candidates[0] if candidates else None


def _contraction_splitter(
    source: Sequence[sp.Expr],
    hull: Sequence[sp.Expr],
    independent: Sequence[sp.Symbol],
    dependent: Sequence[sp.Symbol],
    variables: Sequence[sp.Symbol],
) -> sp.Expr | None:
    if _qq_ideal_equal(source, hull, variables):
        return sp.Integer(1)
    factors = []
    for generator in hull:
        h = _parameter_annihilator(source, generator, independent, dependent, variables)
        if h is None or h == 0:
            return None
        factors.append(h)
    product_h = sp.expand(sp.prod(factors))
    if product_h == 0:
        return None
    if independent:
        _content, primitive = sp.Poly(product_h, *independent, domain=sp.QQ).primitive()
        product_h = sp.expand(primitive.as_expr())
    return product_h


def _component_from_zero_dim(component, variables):
    from ..algebraic_decomposition import PrimaryComponent

    return PrimaryComponent(
        _canonical_qq_basis(component.generators, variables),
        _canonical_qq_basis(component.radical, variables),
        0,
        int(component.degree),
        True,
        "gtz_zero_dimensional_primary",
    )


def _recursive_gtz(
    source_generators: Sequence[sp.Expr],
    variables: Sequence[sp.Symbol],
    memo: dict[tuple[sp.Expr, ...], tuple[tuple[object, ...], GTZNodeCertificate]] | None = None,
) -> tuple[tuple[object, ...], GTZNodeCertificate]:
    from ..algebraic_decomposition import PrimaryComponent

    vars_ = tuple(variables)
    source = _canonical_qq_basis(source_generators, vars_)
    if memo is not None and source in memo:
        return memo[source]
    context = EqualityIdealContext.build(source, vars_)
    dimension = int(context.dimension)
    if context.inconsistent:
        result = (tuple(), GTZNodeCertificate(vars_, source, dimension, "unit"))
        if memo is not None:
            memo[source] = result
        return result
    if dimension == 0:
        zero = zero_dimensional_primary_decomposition(source, vars_)
        if not zero.complete or zero.certificate is None:
            raise ArithmeticError("zero-dimensional coefficient-field decomposition was incomplete")
        components = tuple(
            _component_from_zero_dim(component, vars_) for component in zero.components
        )
        result = (
            components,
            GTZNodeCertificate(
                vars_,
                source,
                0,
                "zero_dimensional",
                zero_dimensional_certificate=zero.certificate,
            ),
        )
        if memo is not None:
            memo[source] = result
        return result

    localization = certified_independent_localization(source, vars_)
    local = _compute_fraction_field_zero_dimensional(
        localization.localized_basis,
        localization.dependent_variables,
        localization.independent_variables,
    )
    contracted: list[GTZContractedComponentCertificate] = []
    top_components: list[PrimaryComponent] = []
    for local_component in local.components:
        q_contraction = contract_localized_ideal(
            local_component.generators,
            localization.independent_variables,
            localization.dependent_variables,
        )
        p_contraction = contract_localized_ideal(
            local_component.radical,
            localization.independent_variables,
            localization.dependent_variables,
        )
        equations = _canonical_qq_basis(q_contraction.generators, vars_)
        radical = _canonical_qq_basis(p_contraction.generators, vars_)
        qctx = EqualityIdealContext.build(equations, vars_)
        pctx = EqualityIdealContext.build(radical, vars_)
        if qctx.dimension != dimension or pctx.dimension != dimension:
            raise ArithmeticError("localized component contraction changed GTZ dimension")
        degree = int(ideal_degree(equations, vars_))
        evidence = GTZContractedComponentCertificate(
            equations,
            radical,
            dimension,
            degree,
            q_contraction.certificate,
            p_contraction.certificate,
        )
        contracted.append(evidence)
        top_components.append(
            PrimaryComponent(
                equations,
                radical,
                dimension,
                degree,
                True,
                "gtz_localized_primary_contraction",
            )
        )

    hull = _intersection_all([component.equations for component in top_components], vars_)
    if _qq_ideal_equal(hull, source, vars_):
        node = GTZNodeCertificate(
            vars_,
            source,
            dimension,
            "localized",
            localization_certificate=localization.certificate,
            local_primary_certificate=local.certificate,
            contracted_components=tuple(contracted),
            hull_generators=hull,
        )
        result = (tuple(top_components), node)
        if memo is not None:
            memo[source] = result
        return result

    splitter = _contraction_splitter(
        source,
        hull,
        localization.independent_variables,
        localization.dependent_variables,
        vars_,
    )
    if splitter in (None, 0, 1):
        raise ArithmeticError("could not construct the GTZ contraction splitter")
    split = saturation_stabilization(source, splitter, vars_)
    if not _qq_ideal_equal(split.generators, hull, vars_):
        raise ArithmeticError("GTZ splitter saturation did not recover the localized hull")
    residual_context = EqualityIdealContext.build(split.companion_generators, vars_)
    if not residual_context.inconsistent and residual_context.dimension > dimension:
        raise ArithmeticError("GTZ residual branch increased dimension")
    if _qq_ideal_equal(split.companion_generators, source, vars_):
        raise ArithmeticError("GTZ residual branch made no ideal progress")
    residual_components, residual_node = _recursive_gtz(split.companion_generators, vars_, memo)
    node = GTZNodeCertificate(
        vars_,
        source,
        dimension,
        "localized_split",
        localization_certificate=localization.certificate,
        local_primary_certificate=local.certificate,
        contracted_components=tuple(contracted),
        hull_generators=hull,
        split_certificate=split.certificate,
        residual=residual_node,
    )
    result = (tuple(top_components) + tuple(residual_components), node)
    if memo is not None:
        memo[source] = result
    return result


def _merge_same_radical(components: Sequence[object], variables: Sequence[sp.Symbol]):
    from ..algebraic_decomposition import PrimaryComponent

    vars_ = tuple(variables)
    groups: dict[tuple[sp.Expr, ...], list[object]] = {}
    for component in components:
        radical = _canonical_qq_basis(component.radical, vars_)
        groups.setdefault(radical, []).append(component)
    merged = []
    for radical, group in groups.items():
        if len(group) == 1:
            merged.append(group[0])
            continue
        equations = _intersection_all([component.equations for component in group], vars_)
        context = EqualityIdealContext.build(equations, vars_)
        merged.append(
            PrimaryComponent(
                equations,
                radical,
                int(context.dimension),
                int(ideal_degree(equations, vars_)),
                True,
                "gtz_same_prime_intersection",
            )
        )
    return tuple(merged)


def _remove_redundant(components: Sequence[object], variables: Sequence[sp.Symbol]):
    vars_ = tuple(variables)
    kept = list(components)
    changed = True
    while changed and len(kept) > 1:
        changed = False
        for index, component in enumerate(tuple(kept)):
            others = [c.equations for j, c in enumerate(kept) if j != index]
            intersection = _intersection_all(others, vars_)
            qctx = EqualityIdealContext.build(component.equations, vars_)
            if all(qctx.normal_form(g) == 0 for g in intersection):
                kept.pop(index)
                changed = True
                break
    return tuple(kept)


def _cleanup_components(components: Sequence[object], variables: Sequence[sp.Symbol]):
    return _remove_redundant(_merge_same_radical(components, variables), variables)


def gtz_primary_decomposition(
    generators: Sequence[sp.Expr], variables: Sequence[sp.Symbol]
) -> GTZPrimaryDecompositionResult:
    """Compute a full recursive GTZ primary decomposition over ``QQ``."""
    vars_ = tuple(variables)
    source = _canonical_qq_basis(generators, vars_)
    try:
        raw_components, root = _recursive_gtz(source, vars_, {})
        components = _cleanup_components(raw_components, vars_)
        reconstructed = _intersection_all([c.equations for c in components], vars_)
        if not _qq_ideal_equal(reconstructed, source, vars_):
            raise ArithmeticError("GTZ components failed final exact reconstruction")
        radicals = [tuple(_canonical_qq_basis(c.radical, vars_)) for c in components]
        irredundant = len(set(radicals)) == len(radicals)
        certificate = GTZPrimaryDecompositionCertificate(
            vars_, source, root, tuple(components), irredundant
        )
        plan = plan_gtz_primary_decomposition(source, vars_)
        result = GTZPrimaryDecompositionResult(
            vars_, tuple(components), True, irredundant, "gtz_recursive", certificate, plan
        )
        if not verify_gtz_primary_decomposition_certificate(certificate):
            raise ArithmeticError("internal recursive GTZ certificate verification failed")
        return result
    except (
        ArithmeticError,
        ValueError,
        TypeError,
        ZeroDivisionError,
        sp.PolynomialError,
        sp.polys.polyerrors.CoercionFailed,
    ):
        return GTZPrimaryDecompositionResult(
            vars_,
            tuple(),
            False,
            False,
            "gtz_incomplete",
            None,
            plan_gtz_primary_decomposition(source, vars_),
        )


def clear_gtz_caches() -> None:
    """Clear process-local GTZ canonical-basis caches."""
    _canonical_qq_basis_cached.cache_clear()
    _canonical_ff_basis_cached.cache_clear()


def gtz_cache_info():
    """Return canonical-basis cache statistics for production profiling."""
    return {
        "qq": _canonical_qq_basis_cached.cache_info(),
        "fraction_field": _canonical_ff_basis_cached.cache_info(),
    }


__all__ = [
    "FractionFieldZeroDimensionalCertificate",
    "FractionFieldZeroDimensionalResult",
    "GTZContractedComponentCertificate",
    "GTZNodeCertificate",
    "GTZPrimaryDecompositionCertificate",
    "GTZExecutionPlan",
    "GTZPrimaryDecompositionResult",
    "plan_gtz_primary_decomposition",
    "LocalPrimaryFactor",
    "clear_gtz_caches",
    "gtz_cache_info",
    "gtz_primary_decomposition",
    "verify_fraction_field_zero_dimensional_certificate",
    "verify_gtz_primary_decomposition_certificate",
]

# Replay/validation is separated from recursive GTZ realization.
from ._gtz_primary_certification import (  # noqa: E402
    verify_gtz_primary_decomposition_certificate,
)
