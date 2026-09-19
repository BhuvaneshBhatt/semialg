# Tutorial: factorized projection and partial CAD

CAD is a semantic decomposition, not a promise to preserve the syntactic form of every input polynomial.

Consider

\[
f=(x-1)(x+1)(y-1).
\]

The complete Collins projection builder first forms an exact squarefree basis. Internally it may retain the factors `x - 1`, `x + 1`, and `y - 1` rather than the redundant product `f`.

```python
import sympy as sp
from semialg.cad_algorithms.projection.collins import build_collins_proj_set

x, y = sp.symbols("x y")
f = (x - 1) * (x + 1) * (y - 1)
tower = build_collins_proj_set((f,), (x, y))
print([p.as_expr() for p in tower.original_polynomials])
```

That representation is exact. On every cell, the sign of `f` is the product of the factor signs (with multiplicity), while pairwise resultants preserve the cross-root events needed during projection. Tests and downstream code should therefore ask whether the **source formula has the correct sign/truth on a cell**, not whether the literal expanded source polynomial is a sign-table key.

## Lazy/partial CAD

Existential feasibility often does not require constructing every terminal cell. The lazy engine can stop as soon as a true branch provides a certified witness; proving UNSAT still requires exhaustive relevant traversal.

```python
from semialg.formula import parse_formula
from semialg.partial.qe import lazy_find_inst_form

formula = sp.And(x**2 < 1, y**2 < 1)
result = lazy_find_inst_form((x, y), parse_formula(formula))
assert result.found
print(result.stats.evaluated_leaf_cells)
```

This is a scalability optimization with asymmetric logic: early success is valid for `exists`, while a negative result cannot be inferred merely because some branches were skipped.

## Reduced CAD and fallback

Equational-constraint/reduced projection can remove additional projection work, but only when its well-orientedness and related side conditions are certified. If those checks fail, `semialg` falls back to a complete projection/lifting route rather than returning a result under unproved assumptions.

The key distinction throughout is:

```text
source formula semantics     <- public mathematical contract
projection/sign-table layout <- internal exact representation
```
