import sympy as sp

from semialg.algebraic_function_fields import (
    MonogenicFunctionField,
    RationalFunctionField,
    _trager_shift_is_good,
    field_element_to_expr,
)


def test_monogenic_inverse_is_exact():
    x, y = sp.symbols("x y")
    base = RationalFunctionField((x,))
    field = MonogenicFunctionField(base, y, (-x, 0, 1))
    inverse = field.inv(field.alpha + 1)
    assert field.mul(field.alpha + 1, inverse) == field.one


def test_trager_factorization_splits_over_algebraic_function_field():
    x, y = sp.symbols("x y")
    base = RationalFunctionField((x,))
    field = MonogenicFunctionField(base, y, (-x, 0, 1))
    # z^2-x = z^2-y^2 in QQ(x)(y).
    _unit, factors = field.factor_univariate((-x, 0, 1))
    assert len(factors) == 2
    exprs = [
        sp.expand(
            sum(field_element_to_expr(c, field) * sp.Symbol("z") ** i for i, c in enumerate(factor))
        )
        for factor, multiplicity in factors
        if multiplicity == 1
    ]
    z = sp.Symbol("z")
    residual = sp.expand((exprs[0] * exprs[1]) - (z**2 - x))
    assert sp.rem(sp.Poly(residual, y), sp.Poly(y**2 - x, y)).as_expr() == 0


def test_trager_bad_shift_guard_replaces_bounded_search():
    x, y = sp.symbols("x y")
    base = RationalFunctionField((x,))
    field = MonogenicFunctionField(base, y, (-x, 0, 1))
    coeffs = tuple(field.convert(c) for c in (-x, 0, 1))

    # At s=0 the norm is (z**2-x)**2, so the discriminant vanishes.
    assert _trager_shift_is_good(coeffs, field, 0) is None


def _mul_coeffs(left, right, field):
    out = [field.zero] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            out[i + j] = out[i + j] + a * b
    return tuple(out)


def _pow_coeffs(poly, exponent, field):
    out = (field.one,)
    for _ in range(exponent):
        out = _mul_coeffs(out, poly, field)
    return out


def test_squarefree_factorization_preserves_multiplicities_and_unit():
    x, y = sp.symbols("x y")
    base = RationalFunctionField((x,))
    field = MonogenicFunctionField(base, y, (-x, 0, 1))

    minus = (-field.alpha, field.one)  # z-y
    plus = (field.alpha, field.one)  # z+y
    polynomial = _mul_coeffs(
        _pow_coeffs(minus, 2, field),
        _pow_coeffs(plus, 3, field),
        field,
    )
    polynomial = tuple(field.convert(7) * coefficient for coefficient in polynomial)

    sf_unit, squarefree_parts = field.squarefree_decomposition(polynomial)
    assert sf_unit == field.convert(7)
    assert sorted(multiplicity for _factor, multiplicity in squarefree_parts) == [2, 3]

    unit, factors = field.factor_univariate(polynomial)

    assert unit == field.convert(7)
    assert sorted(multiplicity for _factor, multiplicity in factors) == [2, 3]

    reconstructed = (unit,)
    for factor, multiplicity in factors:
        reconstructed = _mul_coeffs(
            reconstructed,
            _pow_coeffs(factor, multiplicity, field),
            field,
        )
    assert reconstructed == polynomial


def test_squarefree_decomposition_precedes_trager_on_repeated_nonlinear_factor():
    x, y, z = sp.symbols("x y z")
    base = RationalFunctionField((x,))
    field = MonogenicFunctionField(base, y, (-x, 0, 1))

    # (z**2-y)**2 is not squarefree, while z**2-y is irreducible over
    # QQ(x)(y) when y**2=x.  The outer factorizer must extract multiplicity 2
    # before invoking the squarefree-only Trager core.
    quadratic = (-field.alpha, field.zero, field.one)
    polynomial = _pow_coeffs(quadratic, 2, field)
    unit, factors = field.factor_univariate(polynomial)

    assert unit == field.one
    assert len(factors) == 1
    factor, multiplicity = factors[0]
    assert multiplicity == 2
    expr = sum(field_element_to_expr(c, field) * z**i for i, c in enumerate(factor))
    assert sp.expand(expr - (z**2 - y)) == 0


def test_factor_univariate_constant_and_zero_contracts():
    x, y = sp.symbols("x y")
    base = RationalFunctionField((x,))
    field = MonogenicFunctionField(base, y, (-x, 0, 1))

    unit, factors = field.factor_univariate((5,))
    assert unit == field.convert(5)
    assert factors == ()

    import pytest

    with pytest.raises(Exception, match="zero polynomial"):
        field.factor_univariate(())


def test_certified_factorization_reconstructs_exactly():
    from semialg.algebraic_function_fields import certified_factor_univariate

    x, y = sp.symbols("x y")
    base = RationalFunctionField((x,))
    field = MonogenicFunctionField(base, y, (-x, 0, 1))
    result = certified_factor_univariate(field, (-x, 0, 1))
    assert result.reconstruction_verified
    assert len(result.factors) == 2


def test_primitive_element_compresses_two_nonlinear_stages():
    from semialg.algebraic_function_fields import compress_primitive_element

    x, y, z = sp.symbols("x y z")
    base = RationalFunctionField((x,))
    first = MonogenicFunctionField(base, y, (-x, 0, 1))
    second = MonogenicFunctionField(first, z, (-first.alpha, first.zero, first.one))
    compression = compress_primitive_element(second)
    assert compression.complete
    assert isinstance(compression.field, MonogenicFunctionField)
    assert compression.field.base == base
    assert compression.field.degree == 4
    iy = compression.image(y)
    iz = compression.image(z)
    assert compression.field.is_zero(iy**2 - compression.field.convert(x))
    assert compression.field.is_zero(iz**2 - iy)


def test_primitive_guard_characterizes_bad_zero_parameter_and_replays():
    from dataclasses import replace

    from semialg.algebraic_function_fields import (
        compress_primitive_element,
        verify_primitive_element_compression,
    )

    x, y, z = sp.symbols("x y z")
    base = RationalFunctionField((x,))
    first = MonogenicFunctionField(base, y, (-x, 0, 1))
    second = MonogenicFunctionField(first, z, (-first.alpha, first.zero, first.one))
    compression = compress_primitive_element(second)
    assert compression.steps
    step = compression.steps[0]
    assert sp.cancel(step.bad_parameter_guard.subs(step.parameter_symbol, 0)) == 0
    assert sp.cancel(step.bad_parameter_guard.subs(step.parameter_symbol, step.parameter)) != 0
    assert verify_primitive_element_compression(compression)

    forged_step = replace(step, parameter=0)
    forged = replace(compression, steps=(forged_step,))
    assert not verify_primitive_element_compression(forged)

    forged_guard = replace(step, bad_parameter_guard=step.bad_parameter_guard + 1)
    forged = replace(compression, steps=(forged_guard,))
    assert not verify_primitive_element_compression(forged)


def test_primitive_compression_bidirectional_transport():
    from semialg.algebraic_function_fields import compress_primitive_element

    x, y, z = sp.symbols("x y z")
    base = RationalFunctionField((x,))
    first = MonogenicFunctionField(base, y, (-x, 0, 1))
    second = MonogenicFunctionField(first, z, (-first.alpha, first.zero, first.one))
    compression = compress_primitive_element(second)
    value = second.alpha + second.convert(first.alpha) + second.convert(x)
    roundtrip = compression.to_source(compression.to_compressed(value))
    assert second.is_zero(roundtrip - value)


def test_deeper_tower_compression_and_automatic_heuristic():
    from semialg.algebraic_function_fields import (
        compress_primitive_element,
        should_compress_tower,
        tower_degree,
        tower_depth,
        verify_primitive_element_compression,
    )

    x, y, z, w = sp.symbols("x y z w")
    base = RationalFunctionField((x,))
    first = MonogenicFunctionField(base, y, (-x, 0, 1))
    second = MonogenicFunctionField(first, z, (-first.alpha, first.zero, first.one))
    # A third degree-one layer keeps total degree modest while exercising the
    # arbitrary-depth transport/compression machinery.
    third = MonogenicFunctionField(second, w, (-second.alpha, second.one))
    assert tower_depth(third) == 3
    assert tower_degree(third) == 4
    assert should_compress_tower(third)
    compression = compress_primitive_element(third)
    assert compression.field.degree == 4
    assert len(compression.steps) == 2
    assert verify_primitive_element_compression(compression)


def test_factorization_can_use_automatic_primitive_compression():
    x, y, z, w = sp.symbols("x y z w")
    base = RationalFunctionField((x,))
    first = MonogenicFunctionField(base, y, (-x, 0, 1))
    second = MonogenicFunctionField(first, z, (-first.alpha, first.one))
    third = MonogenicFunctionField(second, w, (-second.alpha, second.one))
    # The depth-three tower has total degree two, so compression removes tower
    # overhead without inflating the finite extension.  u^2-x then splits.
    unit, factors = third.factor_univariate((-x, third.zero, third.one))
    reconstructed = [third.convert(unit)]
    from semialg.algebraic_function_fields import _poly_mul

    for factor, multiplicity in factors:
        for _ in range(multiplicity):
            reconstructed = _poly_mul(reconstructed, list(factor), third)
    target = [third.convert(-x), third.zero, third.one]
    assert len(reconstructed) == len(target)
    assert all(third.is_zero(a - b) for a, b in zip(reconstructed, target, strict=True))
    assert len(factors) == 2


def test_primitive_compression_handles_five_layer_tower():
    from semialg.algebraic_function_fields import (
        compress_primitive_element,
        tower_degree,
        tower_depth,
        verify_primitive_element_compression,
    )

    x = sp.symbols("x")
    generators = sp.symbols("a0:5")
    base = RationalFunctionField((x,))
    field = MonogenicFunctionField(base, generators[0], (-x, 0, 1))
    for generator in generators[1:]:
        field = MonogenicFunctionField(field, generator, (-field.alpha, field.one))
    assert tower_depth(field) == 5
    assert tower_degree(field) == 2
    compression = compress_primitive_element(field)
    assert len(compression.steps) == 4
    assert compression.field.degree == 2
    assert verify_primitive_element_compression(compression)
