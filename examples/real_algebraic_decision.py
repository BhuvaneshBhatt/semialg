"""Certified algebraic feasibility, decomposition, and polynomial sign examples."""

import sympy as sp

from semialg import (
    find_negative_point,
    find_negative_witness_fast,
    polynomial_nonnegative,
    real_algebraic_feasibility,
)
from semialg.algebraic_decomposition import (
    equidimensional_decomposition,
    verify_decomposition_certificate,
)

x, y = sp.symbols("x y", real=True)

circle = real_algebraic_feasibility((x**2 + y**2 - 1,), (x, y), return_result=True)
assert circle.complete and circle.satisfiable
assert circle.assignment is not None
assert sp.simplify((x**2 + y**2 - 1).subs(circle.assignment)) == 0
print("circle witness:", circle.assignment)

empty = real_algebraic_feasibility((x**2 + y**2 + 1,), (x, y), return_result=True)
assert empty.complete and empty.satisfiable is False
print("empty real variety certified:", empty.complete)

decomposition = equidimensional_decomposition((x * y,), (x, y))
assert decomposition.complete and decomposition.certificate is not None
assert verify_decomposition_certificate(decomposition.certificate)
print("decomposition methods:", decomposition.certificate.methods)

assert polynomial_nonnegative(x**4 + y**4 + 1, (x, y)) is True
negative = find_negative_point(x**2 + y**2 - 1, (x, y))
assert negative is not None
print("negative point:", negative)

fast = find_negative_witness_fast(x**3 + y**2, (x, y))
assert fast.found and fast.assignment is not None
print("fast exact witness:", fast.assignment)
