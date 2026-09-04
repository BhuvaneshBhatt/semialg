"""Executable companion to the semialg documentation example."""

import sympy as sp

from semialg import (
    function_range,
    integrate_over_region,
    semialgebraic_measure,
    semialgebraic_minimize,
)

x, a = sp.symbols("x a", real=True)
region = (x >= 0) & (x <= a)

length = semialgebraic_measure(region, [x], parameters=[a], return_stratified=True)
integral = integrate_over_region(x**2, region, [x], parameters=[a], return_stratified=True)
range_result = function_range(
    x,
    region,
    [x],
    parameters=[a],
    return_stratified=True,
    eliminate_quantifiers=True,
)
minimum = semialgebraic_minimize(
    x,
    region,
    [x],
    parameters=[a],
    return_stratified=True,
    eliminate_quantifiers=True,
)

assert length.select({a: -1}) == 0
assert length.select({a: 3}) == 3
assert integral.select({a: 3}) == 9

r3 = range_result.select({a: 3})
m3 = minimum.select({a: 3})
assert r3.quantifier_free and r3.certified
assert m3.quantifier_free and m3.certified
