# CAD and quantifier elimination

This page is the practical bridge between the conceptual and reference documentation.

- New to CAD? Read [CAD concepts](concepts/cad.md).
- Want to understand dispatch among presolve, virtual substitution, RUR, and CAD? Read [How semialg chooses an algorithm](concepts/algorithm_selection.md).
- Looking up functions and result types? Use the [Decision and QE reference](reference/decision_and_qe.md) and [CAD reference](reference/cad.md).
- Tuning a slow computation? See the [Performance guide](guides/performance.md).

## First-class quantified expressions

For programmatic formulas, prefer `Exists` and `ForAll` nodes rather than manually encoding quantifier-prefix tuples.

```python
import sympy as sp
from semialg import Exists, ForAll
from semialg.solve import reduce_complete_expr

x, y = sp.symbols("x y", real=True)

formula = ForAll(x, Exists(y, sp.Eq(x + y, 0)))
reduce_complete_expr(formula)
# True
```

Text interfaces remain supported, but programmatic quantifier nodes preserve Symbol identity and compose more safely.

## Virtual substitution before CAD

Quadratic virtual substitution is an exact backend for the supported low-degree fragment. The high-level planner can try it before CAD for QE and witness tasks. If it declines, that is not a mathematical failure; another exact backend can still solve the formula.

## Variable ordering

CAD is highly order-sensitive. Automatic planning preserves quantifier-block semantics while reordering variables only where logically legal. Use `suggest_variable_order` / `suggest_cad_variable_order` for diagnostics rather than guessing an order solely from printed expression size.

## Finite equality varieties

When common polynomial equalities define a finite complex variety, automatic CAD may use the [Gröbner variety CAD](concepts/groebner_variety_cad.md) backend. It follows only compatible algebraic sections and reconstructs the exact finite solution set. Use `return_result=True` to inspect `result.cad.backend` and the `variety_only`, `equality_dimension`, and `quotient_dimension` diagnostics.

```python
import sympy as sp
from semialg import cad

x, y = sp.symbols("x y", real=True)
formula = sp.And(sp.Eq(x - y, 0), sp.Eq(x**2 + y**2, 2), x > 0)

auto = cad(formula, (x, y), return_result=True)
full = cad(formula, (x, y), strategy="collins", return_result=True)

assert auto.cad.backend == "groebner-variety"
assert full.cad.backend == "collins-complete"
assert auto.formula == sp.And(sp.Eq(x, 1), sp.Eq(y, 1))
```

The specialized backend has an exact fallback contract: inability to certify a fiber or algebraic sign causes a complete-CAD retry, not rejection of a candidate section. Pure existential QE may use this specialization; universal and mixed prefixes retain full-space CAD semantics.

## Reduced projection and fallback

Reduced/equational-constraint paths are used only when their logical and invariance requirements are established. The package falls back to the complete Collins-style path rather than treating a merely algebraic resultant as a logically necessary equality.

## Caching

Exact computation contexts and bounded process-local caches reuse projection, root, sign, specialization, comparison, and RUR work. Cache identity uses structural SymPy/polynomial keys rather than string serialization.

## Exactness

CAD/QE correctness depends on exact root ordering and sign/truth invariance. Fixed-precision numerical comparisons are not used as hidden substitutes on certified paths. See [Exactness and certification](concepts/exactness_and_certification.md).
