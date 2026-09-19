from __future__ import annotations

import sympy as sp

from semialg import Exists, ForAll


def test_quantifier_free_symbols_and_bound_symbols():
    x, y, z = sp.symbols("x y z", real=True)
    formula = Exists((x, y), sp.Eq(x + y + z, 0))

    assert formula.variables == (x, y)
    assert formula.bound_symbols == (x, y)
    assert formula.free_symbols == {z}
    assert formula.formula == sp.Eq(x + y + z, 0)


def test_quantifiers_are_sympy_boolean_nodes():
    x, y = sp.symbols("x y", real=True)
    quantified = Exists(x, sp.Eq(x + y, 0))

    combined = sp.And(quantified, y > 0)
    assert quantified in combined.args
    assert combined.free_symbols == {y}
    assert str(quantified) == "Exists(x, Eq(x + y, 0))"
    assert str(ForAll(x, x > y)) == "ForAll(x, x > y)"


def test_vacuous_and_constant_quantifiers_simplify():
    x, y = sp.symbols("x y", real=True)

    assert Exists(x, y > 0) == (y > 0)
    assert ForAll(x, sp.true) is sp.true
    assert Exists((), sp.Eq(y, 0)) == sp.Eq(y, 0)


def test_substitution_does_not_replace_bound_variable():
    x, y, z = sp.symbols("x y z", real=True)
    quantified = Exists(x, sp.Eq(x + y, 0))

    assert quantified.subs(x, z) == quantified
    assert quantified.subs(y, z) == Exists(x, sp.Eq(x + z, 0))


def test_substitution_alpha_renames_to_prevent_capture():
    x, y = sp.symbols("x y", real=True)
    quantified = Exists(x, sp.Eq(x + y, 0))

    replaced = quantified.subs(y, x)

    assert isinstance(replaced, Exists)
    assert replaced.free_symbols == {x}
    assert replaced.variables[0] != x
    assert isinstance(replaced.variables[0], sp.Dummy)
    assert replaced.formula == sp.Eq(replaced.variables[0] + x, 0)


def test_nested_quantifiers_preserve_lexical_scope():
    x, y, z = sp.symbols("x y z", real=True)
    quantified = ForAll(x, Exists(y, sp.Eq(x + y, z)))

    assert quantified.free_symbols == {z}
    replaced = quantified.subs(z, x)
    assert replaced.free_symbols == {x}
    assert replaced.variables[0] != x
    assert isinstance(replaced.variables[0], sp.Dummy)


def test_apply_and_split_quantifiers_round_trip():
    from semialg import apply_quantifiers
    from semialg.quantifiers import split_quantifiers

    x, y, z = sp.symbols("x y z", real=True)
    matrix = sp.Eq(x + y, z)
    prefix = (("forall", z), ("exists", x), ("exists", y))
    quantified = apply_quantifiers(matrix, prefix)

    got_prefix, got_matrix = split_quantifiers(quantified)
    assert got_prefix == prefix
    assert got_matrix == matrix


def test_parsed_prenex_formula_exposes_quantified_expression():
    from semialg.formula import parse_quant_form_text, parse_quantified_expr

    x, y = sp.symbols("x y", real=True)
    parsed = parse_quant_form_text("forall x. exists y. x + y = 0", symbols={"x": x, "y": y})
    assert parsed.quantified_expr == ForAll(x, Exists(y, sp.Eq(x + y, 0)))

    reparsed = parse_quantified_expr(parsed.quantified_expr)
    assert reparsed.quantifiers == (("forall", x), ("exists", y))
    assert reparsed.matrix_expr == sp.Eq(x + y, 0)


def test_complete_qe_accepts_quantified_expression():
    from semialg.solve import reduce_complete_expr

    x = sp.Symbol("x", real=True)
    result = reduce_complete_expr(ForAll(x, x**2 >= 0))
    assert result is sp.true or result == sp.true
