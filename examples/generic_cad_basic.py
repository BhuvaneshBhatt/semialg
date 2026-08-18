from __future__ import annotations

import sympy as sp

from semialg import generic_cad

a, x = sp.symbols("a x", real=True)
result = generic_cad(sp.Eq(a * x - 1, 0), variables=[x], parameters=[a])

print("generic:", result.generic_formula)
print("exceptional:", result.exceptional_formula)
for case in result.cases:
    print(
        case.param_condition,
        "exceptional=" + str(case.exceptional),
        "solution=",
        case.solution_formula,
    )
