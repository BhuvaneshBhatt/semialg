import sympy as sp

from semialg.normalization import (
    normalize_problem_variables,
    normalize_sampling_variables,
    normalize_symbol_sequence,
)


def test_problem_variables_append_context_symbols():
    x, y = sp.symbols("x y", real=True)
    assert normalize_problem_variables((x,), x + y) == (x, y)


def test_sampling_variables_respect_explicit_variable_list():
    x, y = sp.symbols("x y", real=True)
    assert normalize_sampling_variables((x,), x + y) == (x,)
    assert normalize_sampling_variables(None, x + y) == (x, y)


def test_symbol_sequence_preserves_order_and_deduplicates():
    x, y = sp.symbols("x y", real=True)
    assert normalize_symbol_sequence((y, x, y)) == (y, x)
