from __future__ import annotations

import sympy as sp

from semialg import parametric_cad
from semialg.instances.real_fallbacks import satisfies_formula


def test_parametric_cad_samples_satisfy_their_stratum_conditions() -> None:
    x, a = sp.symbols("x a", real=True)
    result = parametric_cad(sp.And(a > 3, sp.Eq(x**2, a)), (x,), parameters=(a,))
    for stratum in result.strata:
        assert satisfies_formula(stratum.condition, stratum.sample, strict=True)


def test_parametric_cad_feasible_samples_respect_nonzero_parameter_boundary() -> None:
    x, a = sp.symbols("x a", real=True)
    result = parametric_cad(sp.And(a > 4, sp.Eq(x**2, a)), (x,), parameters=(a,))
    feasible = [stratum for stratum in result.strata if stratum.has_solution]
    assert feasible
    assert all(stratum.sample[a] > 4 for stratum in feasible)


def test_parametric_cad_specialized_formula_uses_certified_sample() -> None:
    x, a = sp.symbols("x a", real=True)
    result = parametric_cad(
        sp.And(a > 4, sp.Eq(x**2, a)),
        (x,),
        parameters=(a,),
        specialize_fibers=False,
    )
    stratum = next(stratum for stratum in result.strata if stratum.has_solution)
    assert satisfies_formula(stratum.condition, stratum.sample, strict=True)
    assert stratum.specialized_formula == stratum.solution_formula.subs(dict(stratum.sample))
    assert stratum.solution is None
