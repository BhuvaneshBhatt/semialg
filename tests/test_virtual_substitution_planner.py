import sympy as sp

from semialg.formula import parse_quant_form_text
from semialg.qe.virtual_substitution import (
    VirtualSubstitutionQEResult,
    try_quadratic_virtual_substitution_qe,
)
from semialg.solve import reduce_text


def test_planner_uses_virtual_substitution_for_existential_quadratic_disk_projection():
    x, y = sp.symbols("x y", real=True)

    result = reduce_text(
        "exists x. x^2 + y^2 <= 1",
        symbols={"x": x, "y": y},
        variable_order=[y, x],
        strategy="auto",
        use_preprocess=False,
        return_result=True,
    )

    assert result.method == "quadratic_virtual_substitution"
    assert isinstance(result.metadata["qe_result"], VirtualSubstitutionQEResult)
    assert x not in result.result.free_symbols
    assert sp.simplify(result.result.subs(y, 0)) == sp.true
    assert sp.simplify(result.result.subs(y, 1)) == sp.true
    assert sp.simplify(result.result.subs(y, 2)) == sp.false


def test_virtual_substitution_prepass_falls_back_when_degree_is_too_high():
    x, y = sp.symbols("x y", real=True)

    result = reduce_text(
        "exists x. x^3 + y <= 0",
        symbols={"x": x, "y": y},
        variable_order=[y, x],
        strategy="auto",
        use_preprocess=False,
        return_result=True,
    )

    assert result.method == "cad"
    assert x not in result.result.free_symbols
    assert result.metadata["qe_result"].backend != "quadratic-virtual-substitution-qe"


def test_try_virtual_substitution_qe_removes_vacuous_existential_after_quadratic_step():
    x, y = sp.symbols("x y", real=True)
    parsed = parse_quant_form_text(
        "exists x, y. x^2 <= 1",
        symbols={"x": x, "y": y},
        variable_order=[x, y],
    )

    result = try_quadratic_virtual_substitution_qe(
        parsed.vars, parsed.quantifiers, parsed.matrix_expr
    )

    assert result is not None
    assert result.status == "complete"
    assert result.eliminated_variables == (x, y)
    assert result.formula == sp.true or sp.simplify(result.formula) == sp.true


def test_try_virtual_substitution_qe_returns_none_for_mixed_quantifier_prefix_in_full_mode():
    x, y = sp.symbols("x y", real=True)
    parsed = parse_quant_form_text(
        "forall x. (y^2 + x <= 0)",
        symbols={"x": x, "y": y},
        variable_order=[x, y],
    )

    result = try_quadratic_virtual_substitution_qe(
        parsed.vars, (("forall", x), ("exists", y)), parsed.matrix_expr, full=True
    )

    assert result is None
