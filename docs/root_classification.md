# Root and parameter classification

Root classification is a natural use of CAD/QE because root-existence and root-counting questions can be expressed as real first-order formulas.

## Public API

- `classify_real_roots`
- `solvability_conditions`
- `semialg.parameters.root_count_conditions`

## Examples

```python
import sympy as sp
from semialg import solvability_conditions
from semialg.parameters import root_count_conditions

x, a, b = sp.symbols("x a b", real=True)

solvability_conditions(sp.Eq(x**2 + a * x + b, 0), [x], [a, b])
# a**2 - 4*b >= 0

root_count_conditions(x**2 + a * x + b, x, [a, b])
# 2 roots if a**2 - 4*b > 0
# 1 root if a**2 - 4*b == 0
# 0 roots if a**2 - 4*b < 0
```

## Current scope

Parameterized linear, quadratic, and cubic families are classified exactly. The classifier explicitly stratifies coefficient-induced degree drops, so a cubic family whose leading coefficient vanishes is reclassified on the corresponding quadratic, linear, or constant stratum rather than applying cubic discriminant rules outside their domain.

Quartic families have an exact nonmultiple-root classification. For a quartic with nonzero leading coefficient and nonzero discriminant, the discriminant together with the classical quartic `P` and `D` invariants distinguishes zero, two, and four distinct real roots. Leading-coefficient degree drops recurse through the exact cubic classifier. The discriminant-zero quartic locus contains several distinct multiple-root configurations and remains explicitly unknown (`-1`) unless a stronger specialized argument certifies it.

An unknown root-count stratum is not exposed as a certified value by `return_stratified=True`: the result is marked incomplete and its coverage excludes the unresolved locus. Likewise, an existence query never interprets `-1` as zero roots. `solvability_conditions` falls back to complete QE when root-count classification is insufficient.

For degree five and above, discriminant sign cells are diagnostic only. Their sampled fibers are never promoted to cell-wide exact root counts.

A useful invariant for parameterized root-count code is specialization consistency: for every parameter assignment away from unsupported coefficient domains, the selected exact stratum must agree with the distinct real roots of the specialized polynomial, including assignments on discriminant and leading-coefficient boundaries.

## General parameterized Sturm/subresultant stratification

For parameterized families of degree five and above, `classify_real_roots` and `root_count_conditions` use a general exact stratification when it is constructible. The algorithm forms the subresultant polynomial-remainder sequence of $p$ and $p'$, collects its parameter-dependent coefficient data, and builds a sign-invariant CAD of parameter space. On each CAD cell, the specialized Sturm profile is invariant, so an exact root count at one algebraic sample certifies the **entire cell** rather than extrapolating from a numerical sample.

This handles degree drops, multiple-root strata, and arbitrary polynomial degree without requiring radical formulas. Compact low-degree formulas remain preferred: linear/quadratic/cubic classifiers and the quartic invariant classifier are retained because they are substantially cheaper and produce more readable conditions.

If the subresultant/CAD construction itself cannot be certified, the result remains explicitly partial rather than promoting a sampled count to an exact parameter-wide claim.
