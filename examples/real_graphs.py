"""Check the domain and the selected principal square-root graph separately."""

import sympy as sp

from semialg.function_graph import semialgebraic_function_graph
from semialg.function_verification import verify_real_composition, verify_real_function

x, y, z = sp.symbols("x y z", real=True)
principal = sp.And(sp.Eq(y * y, x), y >= 0)
check = verify_real_function(principal, input=x, output=y, domain=x >= 0)
assert check.total and check.single_valued and check.certified
assert not verify_real_function(principal, input=x, output=y).total
assert not verify_real_function(sp.Eq(y * y, x), input=x, output=y, domain=x > 0).single_valued
target = semialgebraic_function_graph(x, z)
outer = semialgebraic_function_graph(1 / y, z)
inner = semialgebraic_function_graph(1 / x, y)
assert not verify_real_composition(
    target, outer, inner, input=x, output=z, intermediate=y
).equivalent
assert verify_real_composition(
    target, outer, inner, input=x, output=z, intermediate=y, domain=sp.Ne(x, 0)
).equivalent
print("A principal root is a function on x>=0; reciprocal composition retains its hole at zero.")
