"""Exhaust short quantifier prefixes against explicit two-point truth tables."""

from itertools import product

import pytest
import sympy as sp

from semialg import quantifier_eliminate


def _expected(prefix, values=()):
    if not prefix:
        # A dependency chain: the final choice must disagree with the first,
        # and each intermediate bit equals its predecessor.
        return values[-1] != values[0] and all(
            a == b for a, b in zip(values[:-2], values[1:-1], strict=False)
        )
    answers = [_expected(prefix[1:], (*values, value)) for value in (0, 1)]
    return (any if prefix[0] == "exists" else all)(answers)


@pytest.mark.parametrize(
    "prefix",
    [p for n in range(1, 5) for p in product(("exists", "forall"), repeat=n)],
    ids=lambda p: "-".join(p),
)
def test_finite_prefix(prefix):
    variables = sp.symbols(f"x:{len(prefix)}", real=True)
    body = sp.And(
        sp.Ne(variables[-1], variables[0]),
        *(sp.Eq(a, b) for a, b in zip(variables[:-2], variables[1:-1], strict=False)),
    )
    for quantifier, variable in reversed(tuple(zip(prefix, variables, strict=False))):
        guard = sp.Or(sp.Eq(variable, 0), sp.Eq(variable, 1))
        body = sp.And(guard, body) if quantifier == "exists" else sp.Or(sp.Not(guard), body)
    assert bool(
        quantifier_eliminate(body, tuple(zip(prefix, variables, strict=False)))
    ) == _expected(prefix)
