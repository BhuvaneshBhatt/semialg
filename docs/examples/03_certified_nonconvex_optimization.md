# Certified global polynomial optimization

Minimize a polynomial over a curved compact region and inspect attainment and certification.

## Problem and interpretation

The objective is squared distance from \((2,0)\). Its unconstrained minimizer
lies outside the unit disk, so the constrained optimum occurs on the curved
boundary at \((1,0)\).

The example illustrates why the structured result is useful: the exact value,
optimizer witness, attainment status, and global certificate are separate
pieces of information. Candidate generation may use KKT/active-set structure,
while certification rules out every better feasible point exactly.

## Executable example

```python
import sympy as sp

from semialg import semialgebraic_minimize

x, y = sp.symbols("x y", real=True)
disk = x**2 + y**2 <= 1
objective = (x - 2) ** 2 + y**2

result = semialgebraic_minimize(objective, disk, [x, y], return_result=True)

print(result)
assert result.value == 1
assert result.points == ({x: 1, y: 0},)
assert result.attained is True
assert result.certified is True
```

The matching executable file is `examples/gallery/03_certified_nonconvex_optimization.py`; the documentation
test suite runs every gallery script.
