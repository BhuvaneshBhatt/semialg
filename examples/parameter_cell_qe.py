"""Quantified parameter conditions, bounded computation and statement replay."""

import sympy as sp

from semialg import cad, quantifier_eliminate, replay_certificate

parameter, x, y = sp.symbols("parameter x y", real=True)
result = quantifier_eliminate(
    sp.Eq(y**2, x**2 + parameter),
    (("forall", x), ("exists", y)),
    strategy="partial",
    return_result=True,
)
assert result.formula == (parameter >= 0)
proof = result.backend_result
assert proof.truth_by_index == {(0,): False, (1,): True, (2,): True}
assert replay_certificate(proof).verified is True

limited = cad(x > 0, (x,), max_cells=2, return_result=True)
assert limited.status == "unknown"
complete = cad(x > 0, (x,), max_cells=3, return_result=True)
assert complete.as_function()(1)
assert not complete.as_function()(0)
assert replay_certificate(complete).verified is True
