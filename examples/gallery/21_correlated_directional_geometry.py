"""Exact joint images retain correlations lost by scalar ranges."""

import sympy as sp

from semialg import correlated_map_image

t, u, v = sp.symbols("t u v", real=True)
image = correlated_map_image((t, t**2), sp.And(t >= -1, t <= 1), (t,), image_variables=(u, v))
assert sp.simplify(image.formula.subs({u: sp.Rational(1, 2), v: sp.Rational(1, 4)})) is sp.true
assert sp.simplify(image.formula.subs({u: sp.Rational(1, 2), v: sp.Rational(1, 3)})) is sp.false
