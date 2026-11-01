"""Adaptive and fixed control choices have different quantified specifications."""

import sympy as sp

from semialg import quantifier_eliminate

u, d = sp.symbols("u d", real=True)
disturbance = sp.And(d >= -1, d <= 1)
response = sp.And(u >= -1, u <= 1, sp.Eq(u, d))
specification = sp.Or(sp.Not(disturbance), response)
adaptive = quantifier_eliminate(specification, (("forall", d), ("exists", u)))
fixed = quantifier_eliminate(specification, (("exists", u), ("forall", d)))
assert adaptive is sp.true
assert fixed is sp.false
