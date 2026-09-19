# Specialist APIs live in owning namespaces

The package root is reserved for broad mathematical operations. More specialized workflows are imported from the subsystem that owns them, which keeps `semialg.*` smaller without hiding advanced functionality.

```python
import sympy as sp

from semialg.map_degree import parametric_map_degree
from semialg.parameters import root_count_conditions

x, a = sp.symbols("x a", real=True)

counts = root_count_conditions(x**2 - a, x, (a,))
degree = parametric_map_degree((x**2,), (x,))

assert sp.simplify(counts[2] ^ (a > 0)) is sp.false
assert degree.degree == 2
```

Other specialist entry points follow the same rule: roadmap construction lives in `semialg.roadmaps`, Thom encodings in `semialg.algebraic`, certified triangulation/stratification in `semialg.topology.semialgebraic`, and parameter-space covers in `semialg.parametric_geometry`.
