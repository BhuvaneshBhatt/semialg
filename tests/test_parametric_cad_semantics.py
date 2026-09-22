from __future__ import annotations

import sympy as sp

from semialg import equivalent, parametric_cad, solvability_conditions


def test_parametric_cad_strata_cover_parameter_domain() -> None:
    x, a = sp.symbols("x a", real=True)
    result = parametric_cad(sp.Eq(a * x, 1), (x,), parameters=(a,))
    union = sp.Or(*(stratum.condition for stratum in result.strata))
    assert equivalent(union, result.parameter_domain, (a,))


def test_parametric_cad_generic_and_exceptional_strata_are_disjoint() -> None:
    x, a = sp.symbols("x a", real=True)
    result = parametric_cad(sp.Eq(a * x, 1), (x,), parameters=(a,))
    assert equivalent(
        sp.And(result.generic_condition, result.exceptional_condition), sp.false, (a,)
    )


def test_parametric_cad_representative_fibers_preserve_original_formula() -> None:
    x, a = sp.symbols("x a", real=True)
    formula = sp.Eq(x**2, a)
    result = parametric_cad(formula, (x,), parameters=(a,))
    for stratum in result.strata:
        if not stratum.has_solution:
            continue
        expected = formula.subs(dict(stratum.sample))
        assert equivalent(stratum.specialized_formula, expected, (x,))


def test_parametric_cad_nonempty_strata_agree_with_solvability_conditions() -> None:
    x, a = sp.symbols("x a", real=True)
    formula = sp.Eq(a * x, 1)
    result = parametric_cad(formula, (x,), parameters=(a,))
    expected = solvability_conditions(formula, (x,), (a,))
    assert equivalent(result.parameter_condition, expected, (a,))


def test_parametric_cad_exceptional_analysis_is_provenance() -> None:
    x, a = sp.symbols("x a", real=True)
    result = parametric_cad(sp.Eq(a * x**2 - 1, 0), (x,), parameters=(a,))
    assert result.exceptional_analysis is not None
    assert a in result.exceptional_analysis.branching_polynomials
