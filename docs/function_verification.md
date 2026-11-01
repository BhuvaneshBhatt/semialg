# Real graph verification

A relation can be equivalent to another relation without defining a total single-valued function. Check these obligations separately and keep visible denominator exclusions.

<!-- examples: executable -->
```python
import sympy as sp
from semialg.function_verification import verify_real_function

x, y = sp.symbols("x y", real=True)
graph = sp.And(sp.Eq(y**2, x), y >= 0)
result = verify_real_function(graph, input=x, output=y, domain=x >= 0)
assert result.total and result.single_valued and result.certified
assert not verify_real_function(graph, input=x, output=y).total
assert not verify_real_function(sp.Eq(y**2, x), input=x, output=y, domain=x > 0).single_valued
```

`verify_real_composition` compares the target graph with an existentially quantified composition. Graph auxiliaries and the middle coordinate receive fresh scopes. Input domains may contain only input coordinates and declared parameters. Function checks are optional because graph equality and totality are distinct questions.

The [real graph script](https://github.com/BhuvaneshBhatt/semialg/blob/main/examples/real_graphs.py) also verifies reciprocal composition: 1/(1/x) retains a hole at zero even though its simplified expression is x. The [application guide](worked_applications.md) covers parameter and linkage workflows. Unsupported exact decisions raise; numerical sampling is not an equivalence proof. Decomposition-specific adapters belong to [alggeom](https://github.com/BhuvaneshBhatt/alggeom).

## Comparing projected graphs

`compare_real_graphs` compares formulas, graph objects, or existing `QuantifierEliminationResult` records. Existing projections retain their provenance and certification status. Function obligations are independent of graph equality. Set `find_counterexample=True` when an explicit exact witness for a difference is needed; its construction can cost more than a closed existential decision.

```python
from semialg.function_verification import compare_real_graphs
from semialg.qe import quantifier_eliminate

projection = quantifier_eliminate(graph, [], variables=(x, y), return_result=True)
comparison = compare_real_graphs(
    projection, graph, input=x, output=y, domain=x >= 0, check_functions=True
)
assert comparison.equivalent and comparison.certified and comparison.function_verified
assert comparison.eliminations[0] is projection
```

The generic comparison engine and result records live in semialg. alggeom constructs decomposition-specific chains, branch selectors and pole guards, then calls this engine. Generic graph functions are imported directly from `semialg.function_verification`; alggeom's real-decomposition interface exposes `verify_chain` and `verify_decomposition`.
