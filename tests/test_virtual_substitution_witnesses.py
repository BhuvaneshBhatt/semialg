import sympy as sp

from semialg.formula import parse_quant_form_text
from semialg.qe.virtual_substitution import (
    VirtualSubstitutionWitnessResult,
    reconstruct_vs_value,
    try_quadratic_virtual_substitution_witness,
)
from semialg.solve.find_instance import find_instance_text


def _is_true(expr):
    value = expr
    if value == sp.true or value is sp.true or value is True:
        return True
    return bool(value)


def test_reconstructs_value_from_disk_cross_section():
    x, y = sp.symbols("x y", real=True)
    formula = sp.And(x > 0, x**2 + y**2 <= 1, evaluate=False)

    value = reconstruct_vs_value(formula, x, {y: sp.Rational(1, 2)})

    assert value is not None
    assert _is_true(formula.subs({x: value, y: sp.Rational(1, 2)}))


def test_virtual_substitution_witness_for_one_quantified_variable():
    x, y = sp.symbols("x y", real=True)
    parsed = parse_quant_form_text(
        "exists x. (x > 0) & (x^2 + y^2 <= 1) & (y == 1/2)",
        symbols={"x": x, "y": y},
        variable_order=[x, y],
    )

    def base_finder(reduced_formula, variables):
        assert variables == (y,)
        return {y: sp.Rational(1, 2)}

    result = try_quadratic_virtual_substitution_witness(
        parsed.vars, parsed.quantifiers, parsed.matrix_expr, base_finder
    )

    assert isinstance(result, VirtualSubstitutionWitnessResult)
    assert result.instance is not None
    assert x in result.instance and y in result.instance
    assert _is_true(parsed.matrix_expr.subs(result.instance))


def test_find_instance_text_uses_virtual_substitution_witness_backend():
    x, y = sp.symbols("x y", real=True)

    result = find_instance_text(
        "exists x. (x > 0) & (x^2 + y^2 <= 1) & (y == 1/2)",
        symbols={"x": x, "y": y},
        variable_order=[x, y],
        domain="reals",
        strategy="auto",
        use_preprocess=False,
        return_result=True,
    )

    assert result.method == "quadratic_virtual_substitution_instance"
    assert result.found
    inst = result.first()
    assert inst is not None
    formula = sp.And(x > 0, x**2 + y**2 <= 1, sp.Eq(y, sp.Rational(1, 2)), evaluate=False)
    assert _is_true(formula.subs(inst))


def test_virtual_substitution_witness_reconstructs_two_quantified_variables():
    x, y, z = sp.symbols("x y z", real=True)
    parsed = parse_quant_form_text(
        "exists x, y. (x^2 + y^2 <= 1) & (x > 0) & (y > 0) & (z == 0)",
        symbols={"x": x, "y": y, "z": z},
        variable_order=[x, y, z],
    )

    def base_finder(reduced_formula, variables):
        return {z: 0}

    result = try_quadratic_virtual_substitution_witness(
        parsed.vars, parsed.quantifiers, parsed.matrix_expr, base_finder
    )

    assert result is not None
    assert result.instance is not None
    assert result.instance[x] > 0
    assert result.instance[y] > 0
    assert _is_true(parsed.matrix_expr.subs(result.instance))
