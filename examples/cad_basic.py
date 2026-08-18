from __future__ import annotations

import sympy as sp

from semialg import cad

x, y = sp.symbols("x y", real=True)

quadrant = cad((x >= 0) & (y >= 0), [x, y])
print("Quadrant formula:", quadrant.formula)
print("Cell counts:", quadrant.cell_count_by_level())

fn = cad((x >= 0) & (y >= 0), [x, y], output="function", return_result=False)
print("Contains (1, 1):", fn.contains({x: 1, y: 1}))
print("Contains (-1, 1):", fn.contains({x: -1, y: 1}))
print("Projected condition on x:", fn.project([x]).to_formula())
