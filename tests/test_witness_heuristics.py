import sympy as sp

from semialg import find_negative_witness_fast
from semialg.witness_heuristics import odd_degree_negative_witness, random_line_negative_witness


def test_odd_degree_witness_is_exactly_verified():
    x, y = sp.symbols("x y")
    result = odd_degree_negative_witness(x**3 + y**2, (x, y))
    assert result.found
    assert result.assignment is not None
    assert sp.simplify((x**3 + y**2).subs(result.assignment)) < 0


def test_random_line_failure_has_no_unsat_semantics():
    x, y = sp.symbols("x y")
    result = random_line_negative_witness(x**2 + y**2, (x, y), attempts=0)
    assert result.found is False
    assert result.certified is False


def test_fast_witness_search_finds_simple_negative_ray():
    x, y = sp.symbols("x y")
    result = find_negative_witness_fast(x**2 + y**2 - 2, (x, y))
    assert result.found
    assert result.assignment is not None
    assert sp.simplify((x**2 + y**2 - 2).subs(result.assignment)) < 0


def test_random_line_search_validates_tuning_parameters():
    x = sp.symbols("x")
    for kwargs in ({"attempts": -1}, {"coefficient_bound": 0}, {"coefficient_bound": -1}):
        try:
            random_line_negative_witness(x**2 - 1, (x,), **kwargs)
        except ValueError:
            pass
        else:
            raise AssertionError(f"expected ValueError for {kwargs}")


def test_fast_witness_search_rejects_negative_random_line_count():
    x = sp.symbols("x")
    try:
        find_negative_witness_fast(x**2 - 1, (x,), random_lines=-1)
    except ValueError:
        pass
    else:
        raise AssertionError("expected ValueError")


def test_univariate_witness_uses_exact_sector_between_algebraic_roots():
    from semialg.witness_heuristics import _univariate_negative_sample

    x = sp.symbols("x")
    # The negative interval is bounded by irrational roots +/-sqrt(2).
    sample = _univariate_negative_sample(x**2 - 2, x)
    assert sample is not None
    assert sample.is_Rational is True
    assert sp.expand((x**2 - 2).subs(x, sample)) < 0
