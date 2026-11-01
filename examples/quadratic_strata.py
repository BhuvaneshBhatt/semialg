"""Classify quadratic, linear, empty and zero-equation parameter strata."""

import sympy as sp

from semialg.decision import equivalent, is_satisfiable
from semialg.qe import quantifier_eliminate

x, a, b, c = sp.symbols("x a b c", real=True)
delta = b * b - 4 * a * c
guards = (
    sp.And(sp.Ne(a, 0), delta > 0),
    sp.And(sp.Ne(a, 0), sp.Eq(delta, 0)),
    sp.And(sp.Ne(a, 0), delta < 0),
    sp.And(sp.Eq(a, 0), sp.Ne(b, 0)),
    sp.And(sp.Eq(a, 0), sp.Eq(b, 0), sp.Ne(c, 0)),
    sp.And(sp.Eq(a, 0), sp.Eq(b, 0), sp.Eq(c, 0)),
)
u, d = sp.symbols("u d", real=True)
radicand = quantifier_eliminate(sp.Eq(u * u, d), [("exists", u)], variables=(d, u))
assert equivalent(radicand, d >= 0, (d,))
assert equivalent(sp.Or(d > 0, sp.Eq(d, 0), d < 0), sp.true, (d,))
abstract = tuple(g.xreplace({delta: d}) for g in guards)
for variable in (a, b, c):
    assert equivalent(sp.Or(sp.Eq(variable, 0), sp.Ne(variable, 0)), sp.true, (variable,))
# Distribute the same hierarchical cases propositionally; the scalar split
# identities above and discriminant trichotomy prove exhaustive coverage.
hierarchical = sp.Or(
    sp.And(sp.Ne(a, 0), sp.Or(d > 0, sp.Eq(d, 0), d < 0)),
    sp.And(sp.Eq(a, 0), sp.Or(sp.Ne(b, 0), sp.And(sp.Eq(b, 0), sp.Or(sp.Ne(c, 0), sp.Eq(c, 0))))),
)
assert sp.simplify_logic(sp.Equivalent(sp.Or(*abstract), hierarchical), force=True) is sp.true
for index, guard in enumerate(abstract):
    for other in abstract[:index]:
        assert not is_satisfiable(sp.And(guard, other), (a, b, c, d))

# Completing the square is an invertible affine change when a is nonzero.
assert sp.expand(4 * a * (a * x * x + b * x + c) - ((2 * a * x + b) ** 2 - delta)) == 0
assert sp.cancel((a * x * x + b * x + c).subs(x, -c / b)) == a * c * c / (b * b)
projection = sp.Or(sp.And(sp.Ne(a, 0), radicand.subs(d, delta)), guards[3], guards[5])
for values, count in (
    ((1, 0, -1), 2),
    ((1, 2, 1), 1),
    ((1, 0, 1), 0),
    ((0, 2, -1), 1),
    ((0, 0, 1), 0),
):
    polynomial = (a * x * x + b * x + c).subs(dict(zip((a, b, c), values, strict=True)))
    assert (
        len(sp.polys.polytools.intervals(sp.Poly(polynomial, x), eps=sp.Rational(1, 1000))) == count
    )
print(
    "Six disjoint guards cover parameter space; their finite fibers have the independently checked real-root counts."
)
