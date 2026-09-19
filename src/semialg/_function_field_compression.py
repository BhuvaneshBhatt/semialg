"""Primitive-element compression, transport, and replay certification."""

from __future__ import annotations

from itertools import count

import sympy as sp

from .algebraic_function_fields import (
    FunctionFieldError,
    MonogenicElement,
    MonogenicFunctionField,
    PrimitiveElementCompression,
    PrimitiveElementCompressionStep,
    RationalFunctionField,
    evaluate_expression,
    field_element_to_expr,
)


def _tower_levels(field):
    levels = []
    current = field
    while isinstance(current, MonogenicFunctionField):
        levels.append(current)
        current = current.base
    if not isinstance(current, RationalFunctionField):
        raise FunctionFieldError("primitive compression requires a rational-function base")
    levels.reverse()
    return current, tuple(levels)


def _flatten_two_stage_element(value, top: MonogenicFunctionField):
    """Coordinates in the basis alpha^i beta^j over the rational base."""
    lower = top.base
    if not isinstance(lower, MonogenicFunctionField):
        raise FunctionFieldError("expected a two-stage tower")
    value = top.convert(value)
    coords = []
    for beta_coeff in value.coeffs:
        beta_coeff = lower.convert(beta_coeff)
        coords.extend(lower.base.normalize(c) for c in beta_coeff.coeffs)
    return tuple(coords)


def _linear_combination_of_powers(coefficients, theta_powers, top):
    out = top.zero
    for coefficient, power in zip(coefficients, theta_powers, strict=True):
        out = out + top.convert(coefficient) * power
    return out


def _transport_element_to_extended_base(value, source, target):
    """Transport an element while extending only the rational-function base."""
    if isinstance(source, RationalFunctionField):
        return target.normalize(source.normalize(value))
    if not isinstance(target, MonogenicFunctionField):
        raise FunctionFieldError("tower shape changed during base extension")
    value = source.convert(value)
    return MonogenicElement(
        target,
        tuple(
            _transport_element_to_extended_base(c, source.base, target.base) for c in value.coeffs
        ),
    )


def _extend_two_stage_base(lower, top, parameter_symbol):
    """Extend K(alpha)(beta) from K to K(c) exactly."""
    base = lower.base
    if not isinstance(base, RationalFunctionField) or top.base != lower:
        raise FunctionFieldError("expected adjacent two-stage tower over a rational base")
    extended_base = RationalFunctionField(base.parameters + (parameter_symbol,))
    lower_ext = MonogenicFunctionField(
        extended_base,
        lower.generator,
        tuple(extended_base.normalize(c) for c in lower.modulus),
    )
    top_modulus = tuple(
        _transport_element_to_extended_base(c, lower, lower_ext) for c in top.modulus
    )
    top_ext = MonogenicFunctionField(lower_ext, top.generator, top_modulus)
    return extended_base, lower_ext, top_ext


def _primitive_bad_parameter_guard(lower, top):
    """Return the exact determinant guard for ``theta = alpha + c*beta``.

    The power-basis matrix is assembled with the binomial theorem, so the
    transcendental parameter c never enters the algebraic tower arithmetic.
    Its determinant vanishes exactly when theta fails to generate the full
    compositum over the rational-function base.
    """
    if top.base != lower or not isinstance(lower.base, RationalFunctionField):
        raise FunctionFieldError("expected adjacent two-stage tower over a rational base")
    c = sp.Dummy("primitive_c")
    dimension = lower.degree * top.degree
    alpha = MonogenicElement(top, (lower.alpha,) + (lower.zero,) * (top.degree - 1))
    beta = top.alpha
    alpha_powers = [top.one]
    beta_powers = [top.one]
    for _ in range(1, dimension):
        alpha_powers.append(alpha_powers[-1] * alpha)
        beta_powers.append(beta_powers[-1] * beta)

    columns = []
    for j in range(dimension):
        entries = [sp.Integer(0)] * dimension
        for k in range(j + 1):
            term = alpha_powers[j - k] * beta_powers[k]
            coords = _flatten_two_stage_element(term, top)
            scale = sp.binomial(j, k) * c**k
            for row, coeff in enumerate(coords):
                entries[row] += scale * sp.cancel(coeff)
        columns.append([sp.cancel(entry) for entry in entries])
    matrix = sp.Matrix.hstack(*(sp.Matrix(column) for column in columns))
    # Symbolic determinant expansion in c grows rapidly with tower degree.
    # The j-th power column has c-degree at most j, so det(M(c)) has degree
    # at most d(d-1)/2.  For larger composita recover the entire determinant
    # exactly from degree-bound+1 evaluations over the rational-function base.
    degree_bound = dimension * (dimension - 1) // 2
    if dimension <= 3:
        determinant = sp.cancel(matrix.det(method="domain-ge"))
    else:
        points = []
        for value in range(degree_bound + 1):
            specialized = matrix.applyfunc(
                lambda entry, value=value: sp.cancel(entry.subs(c, value))
            )
            points.append((sp.Integer(value), sp.cancel(specialized.det(method="domain-ge"))))
        determinant = sp.cancel(sp.interpolate(points, c))
        check_value = degree_bound + 1
        specialized = matrix.applyfunc(lambda entry: sp.cancel(entry.subs(c, check_value)))
        if sp.cancel(determinant.subs(c, check_value) - specialized.det(method="domain-ge")) != 0:
            raise FunctionFieldError("primitive guard interpolation failed exact verification")
    numerator, _denominator = sp.fraction(sp.together(determinant))
    domain = sp.QQ.frac_field(*lower.base.parameters) if lower.base.parameters else sp.QQ
    guard = sp.Poly(sp.expand(numerator), c, domain=domain)
    if guard.is_zero:
        raise FunctionFieldError("primitive-element guard vanished identically")
    return c, sp.factor(guard.as_expr())


def _first_integer_off_guard(guard, parameter_symbol, base_parameters):
    """Choose a small integer provably outside a nonzero univariate guard."""
    domain = sp.QQ.frac_field(*base_parameters) if base_parameters else sp.QQ
    poly = sp.Poly(sp.expand(guard), parameter_symbol, domain=domain)
    if poly.is_zero:
        raise FunctionFieldError("bad-parameter guard must be nonzero")
    candidates = [0]
    for k in count(1):
        candidates.extend((k, -k))
        while candidates:
            value = candidates.pop(0)
            if sp.cancel(poly.eval(value)) != 0:
                return value


def _merge_two_monogenic_fields(
    lower: MonogenicFunctionField,
    top: MonogenicFunctionField,
    *,
    parameter_override: int | None = None,
    guard_data=None,
):
    """Compress two adjacent simple extensions using an exact bad-parameter guard."""
    if top.base != lower or not isinstance(lower.base, RationalFunctionField):
        raise FunctionFieldError("two-stage compression requires adjacent simple extensions")
    base = lower.base
    dimension = lower.degree * top.degree
    parameter_symbol, guard = (
        _primitive_bad_parameter_guard(lower, top) if guard_data is None else guard_data
    )
    parameter = (
        _first_integer_off_guard(guard, parameter_symbol, base.parameters)
        if parameter_override is None
        else int(parameter_override)
    )
    if sp.cancel(guard.subs(parameter_symbol, parameter)) == 0:
        raise FunctionFieldError("chosen primitive parameter lies on the bad-parameter guard")

    alpha_in_top = MonogenicElement(top, (lower.alpha,) + (lower.zero,) * (top.degree - 1))
    beta_in_top = top.alpha
    theta_old = alpha_in_top + top.convert(parameter) * beta_in_top
    powers = [top.one]
    for _ in range(1, dimension + 1):
        powers.append(powers[-1] * theta_old)
    columns = [
        [sp.cancel(c) for c in _flatten_two_stage_element(power, top)]
        for power in powers[:dimension]
    ]
    matrix = sp.Matrix.hstack(*(sp.Matrix(column) for column in columns))
    determinant = sp.cancel(matrix.det())
    if determinant == 0:
        raise FunctionFieldError("guard specialization did not yield a primitive basis")

    rhs_theta_d = sp.Matrix(_flatten_two_stage_element(powers[dimension], top))
    # Invert the exact change-of-basis matrix once.  Recomputing the inverse
    # for each right-hand side was a major source of repeated symbolic work.
    inverse = matrix.inv()
    relation = inverse * rhs_theta_d
    rhs_alpha = sp.Matrix(_flatten_two_stage_element(alpha_in_top, top))
    rhs_beta = sp.Matrix(_flatten_two_stage_element(beta_in_top, top))
    alpha_coeffs = inverse * rhs_alpha
    beta_coeffs = inverse * rhs_beta

    modulus = tuple(base.normalize(-relation[i]) for i in range(dimension)) + (base.one,)
    new_symbol = sp.Dummy(f"pe_{lower.generator}_{top.generator}")
    new_field = MonogenicFunctionField(base, new_symbol, modulus)
    image_a = sum(
        new_field.convert(base.normalize(alpha_coeffs[i])) * new_field.alpha**i
        for i in range(dimension)
    )
    image_b = sum(
        new_field.convert(base.normalize(beta_coeffs[i])) * new_field.alpha**i
        for i in range(dimension)
    )

    lower_rel = new_field.zero
    for i, c in enumerate(lower.modulus):
        lower_rel = lower_rel + new_field.convert(base.normalize(c)) * image_a**i
    top_rel = new_field.zero
    for i, c in enumerate(top.modulus):
        c_expr = field_element_to_expr(c, lower)
        c_image = evaluate_expression(c_expr, new_field, {lower.generator: image_a})
        top_rel = top_rel + c_image * image_b**i
    if not new_field.is_zero(lower_rel) or not new_field.is_zero(top_rel):
        raise FunctionFieldError("primitive compression failed defining-relation verification")

    theta_image = image_a + new_field.convert(parameter) * image_b
    if not new_field.is_zero(theta_image - new_field.alpha):
        raise FunctionFieldError("primitive compression failed generator verification")
    return new_field, image_a, image_b, parameter_symbol, guard, parameter


def compress_primitive_element(field) -> PrimitiveElementCompression:
    """Compress an arbitrary finite monogenic tower to one simple extension.

    Every adjacent merge constructs the exact bad-parameter guard for
    ``theta = alpha + c*beta`` and chooses a certified integer outside that
    finite exceptional set.  Generator images and the final primitive element
    in the original tower are retained for bidirectional exact transport.
    """
    base, levels = _tower_levels(field)
    if not levels:
        return PrimitiveElementCompression(
            base, tuple(), True, None, "already_rational", source_field=field
        )
    if len(levels) == 1:
        return PrimitiveElementCompression(
            levels[0],
            ((levels[0].generator, levels[0].alpha),),
            True,
            None,
            "already_monogenic",
            primitive_expression=levels[0].generator,
            source_field=field,
        )

    current = levels[0]
    images: dict[sp.Symbol, object] = {levels[0].generator: current.alpha}
    primitive_expression = sp.sympify(levels[0].generator)
    steps = []
    last_parameter = None
    for original_top in levels[1:]:
        original_prev = original_top.base
        transported = []
        for coeff in original_top.modulus:
            expr = field_element_to_expr(coeff, original_prev)
            transported.append(evaluate_expression(expr, current, images))
        rebased_top = MonogenicFunctionField(current, original_top.generator, tuple(transported))
        (
            new_field,
            image_current,
            image_top,
            parameter_symbol,
            guard,
            parameter,
        ) = _merge_two_monogenic_fields(current, rebased_top)
        new_images: dict[sp.Symbol, object] = {}
        for symbol, old_value in images.items():
            expr = field_element_to_expr(old_value, current)
            new_images[symbol] = evaluate_expression(
                expr, new_field, {current.generator: image_current}
            )
        new_images[original_top.generator] = image_top
        primitive_expression = sp.expand(primitive_expression + parameter * original_top.generator)
        steps.append(
            PrimitiveElementCompressionStep(
                current.generator,
                original_top.generator,
                parameter_symbol,
                parameter,
                guard,
                current.degree * rebased_top.degree,
                new_field.degree,
                primitive_expression,
            )
        )
        current = new_field
        images = new_images
        last_parameter = parameter
    return PrimitiveElementCompression(
        current,
        tuple(images.items()),
        True,
        last_parameter,
        "successive_exact_guard_primitive_compression",
        tuple(steps),
        primitive_expression,
        field,
    )


def _same_field_element(a, b, field) -> bool:
    try:
        return field.is_zero(field.sub(a, b))
    except (ArithmeticError, ValueError, TypeError, ZeroDivisionError, AttributeError):
        try:
            return field.is_zero(a - b)
        except (ArithmeticError, ValueError, TypeError, ZeroDivisionError, AttributeError):
            return False


def verify_primitive_element_compression(compression: PrimitiveElementCompression) -> bool:
    """Replay a primitive-element compression independently of the search."""
    if not compression.complete or compression.source_field is None:
        return False
    try:
        base, levels = _tower_levels(compression.source_field)
        if not levels:
            return compression.field == base
        if not isinstance(compression.field, MonogenicFunctionField):
            return False
        expected_degree = 1
        for level in levels:
            expected_degree *= level.degree
        if compression.field.degree != expected_degree or compression.field.base != base:
            return False
        if len(compression.steps) != max(0, len(levels) - 1):
            return False

        # Replay every adjacent merge using the certificate's chosen parameter.
        current = levels[0]
        images: dict[sp.Symbol, object] = {levels[0].generator: current.alpha}
        primitive_expression = sp.sympify(levels[0].generator)
        for original_top, step in zip(levels[1:], compression.steps, strict=True):
            transported = []
            for coeff in original_top.modulus:
                expr = field_element_to_expr(coeff, original_top.base)
                transported.append(evaluate_expression(expr, current, images))
            rebased_top = MonogenicFunctionField(
                current, original_top.generator, tuple(transported)
            )
            recomputed_symbol, recomputed_guard = _primitive_bad_parameter_guard(
                current, rebased_top
            )
            domain = sp.QQ.frac_field(*base.parameters) if base.parameters else sp.QQ
            certified_guard = sp.Poly(
                sp.expand(step.bad_parameter_guard), step.parameter_symbol, domain=domain
            ).monic()
            replay_guard = sp.Poly(
                sp.expand(recomputed_guard.subs(recomputed_symbol, step.parameter_symbol)),
                step.parameter_symbol,
                domain=domain,
            ).monic()
            if certified_guard != replay_guard:
                return False
            if step.degree_before != current.degree * rebased_top.degree:
                return False
            guard_at_parameter = recomputed_guard.subs(recomputed_symbol, step.parameter)
            if sp.cancel(guard_at_parameter) == 0:
                return False
            merged = _merge_two_monogenic_fields(
                current,
                rebased_top,
                parameter_override=step.parameter,
                guard_data=(recomputed_symbol, recomputed_guard),
            )
            new_field, image_current, image_top, _ps, _guard, _parameter = merged
            if new_field.degree != step.degree_after:
                return False
            new_images = {}
            for symbol, old_value in images.items():
                expr = field_element_to_expr(old_value, current)
                new_images[symbol] = evaluate_expression(
                    expr, new_field, {current.generator: image_current}
                )
            new_images[original_top.generator] = image_top
            primitive_expression = sp.expand(
                primitive_expression + step.parameter * original_top.generator
            )
            if sp.expand(primitive_expression - step.primitive_expression) != 0:
                return False
            current, images = new_field, new_images

        # Do not depend on generated Dummy names: compare modulus coefficients.
        if current.degree != compression.field.degree:
            return False
        for a, b in zip(current.modulus, compression.field.modulus, strict=True):
            if not _same_field_element(a, b, base):
                return False
        if sp.expand(primitive_expression - compression.primitive_expression) != 0:
            return False
        for symbol, expected in compression.generator_images:
            if symbol not in images:
                return False
            # Compare after transporting the replayed element to expressions in
            # the corresponding primitive generator.
            lhs = field_element_to_expr(images[symbol], current).subs(
                current.generator, compression.field.generator
            )
            rhs = field_element_to_expr(expected, compression.field)
            if sp.cancel(lhs - rhs) != 0:
                return False

        # Verify all original defining equations after transport.
        image_map = dict(compression.generator_images)
        for level in levels:
            image = image_map[level.generator]
            relation = compression.field.zero
            for i, coeff in enumerate(level.modulus):
                coeff_expr = field_element_to_expr(coeff, level.base)
                coeff_image = evaluate_expression(coeff_expr, compression.field, image_map)
                relation = relation + coeff_image * image**i
            if not compression.field.is_zero(relation):
                return False
        primitive_image = evaluate_expression(
            compression.primitive_expression, compression.field, image_map
        )
        return compression.field.is_zero(primitive_image - compression.field.alpha)
    except (
        ArithmeticError,
        ValueError,
        TypeError,
        ZeroDivisionError,
        KeyError,
        IndexError,
        AttributeError,
        FunctionFieldError,
        sp.PolynomialError,
        sp.polys.polyerrors.CoercionFailed,
    ):
        return False
