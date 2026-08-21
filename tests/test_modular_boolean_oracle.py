import itertools

import pytest
import sympy as sp

from semialg.solve.integer.congruence import combine_mod_crt, solve_quant_free_mod_sys


def brute_points(expr, variables, modulus):
    out = []
    for pt in itertools.product(range(modulus), repeat=len(variables)):
        subs = dict(zip(variables, pt, strict=True))
        truth = expr.subs(subs)
        if isinstance(truth, (sp.logic.boolalg.BooleanTrue, sp.logic.boolalg.BooleanFalse)):
            ok = bool(truth)
        else:
            ok = bool(sp.simplify(truth))
        if ok:
            out.append(pt)
    return out


@pytest.mark.parametrize("modulus", [2, 3, 4, 5, 6])
def test_boolean_formulas_match_bruteforce(modulus):
    x, y = sp.symbols("x y", integer=True)
    atoms = [
        sp.Eq(sp.Mod(x, modulus), 0),
        sp.Ne(sp.Mod(x, modulus), 0),
        sp.Eq(sp.Mod(y, modulus), 1 % modulus),
        sp.Ne(sp.Mod(y, modulus), 1 % modulus),
    ]
    formulas = [
        sp.true,
        sp.false,
        *atoms,
        sp.Or(atoms[0], atoms[2]),
        sp.And(atoms[1], atoms[3]),
        sp.Not(atoms[0]),
        sp.Or(sp.Not(atoms[0]), atoms[2]),
        sp.And(sp.Or(atoms[0], atoms[2]), sp.Not(atoms[3])),
    ]
    for formula in formulas:
        got = solve_quant_free_mod_sys(formula, (x, y), modulus, max_points=10000)
        assert got.complete
        assert sorted(got.points) == sorted(brute_points(formula, (x, y), modulus))


def test_composite_inequality_uses_correct_crt_semantics():
    x = sp.symbols("x", integer=True)
    result = solve_quant_free_mod_sys(sp.Ne(x, 0), (x,), 6)
    assert result.points == [(1,), (2,), (3,), (4,), (5,)]


def test_invalid_modulus_is_not_truncated():
    x = sp.symbols("x", integer=True)
    with pytest.raises(TypeError):
        solve_quant_free_mod_sys(sp.Eq(x, 0), (x,), 2.5)
    with pytest.raises(TypeError):
        solve_quant_free_mod_sys(sp.Eq(x, 0), (x,), True)


def test_crt_recombination_enforces_budget():
    with pytest.raises(RuntimeError, match="CRT recombination exceeds"):
        combine_mod_crt(
            [[(0,), (1,)], [(0,), (1,), (2,)]],
            (sp.Symbol("x"),),
            (2, 3),
            max_points=4,
        )
