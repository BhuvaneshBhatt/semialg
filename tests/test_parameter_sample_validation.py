from __future__ import annotations

import pytest
import sympy as sp

from semialg.instances.real_fallbacks import satisfies_formula
from semialg.parameter_stratification import (
    _make_parameter_stratum,
    _validated_parameter_sample,
)


def test_parameter_sample_replaces_invalid_zero_fallback():
    a = sp.symbols("a", real=True)
    condition = a > 3

    sample = _validated_parameter_sample(condition, (a,), {a: sp.Integer(0)})

    assert sample[a] != 0
    assert satisfies_formula(condition, sample, strict=True)


def test_parameter_sample_fails_conservatively_without_representative():
    a = sp.symbols("a", real=True)

    with pytest.raises(NotImplementedError, match="empty parameter stratum"):
        _validated_parameter_sample(sp.false, (a,))


def test_parameter_stratum_builder_uses_validated_sample_for_specialization():
    x, a = sp.symbols("x a", real=True)
    stratum = _make_parameter_stratum(
        index=0,
        expr=sp.Eq(x**2, a),
        parameters=(a,),
        variables=(x,),
        parameter_cell=None,
        condition=a > 4,
        candidate_sample={a: sp.Integer(0)},
        specialize_fibers=False,
    )

    assert satisfies_formula(stratum.condition, stratum.sample, strict=True)
    assert stratum.sample[a] > 4
    assert stratum.specialized_formula == sp.Eq(x**2, stratum.sample[a])
    assert stratum.solution is None
