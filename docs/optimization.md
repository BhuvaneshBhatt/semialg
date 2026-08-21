# Exact polynomial optimization

`semialg` provides exact global optimization for supported polynomial objectives over real semialgebraic domains.

## Public API

- `semialgebraic_minimize`
- `semialgebraic_maximize`
- `OptimizationResult`
- `ParametricOptimizationResult`
- `OptimizationCertificationPolicy`
- `polynomial_locus_dimension`
- `function_range`
- `FunctionRangeResult`
- `ParametricFunctionRangeResult`

By default the optimization functions return an `OptimizationResult`. Pass `return_result=False` to request only the optimum value.

```python
import sympy as sp
from semialg import semialgebraic_minimize, semialgebraic_maximize

x, y, z = sp.symbols("x y z", real=True)

semialgebraic_minimize(x**2, [x >= 2], [x]).value
# 4

semialgebraic_minimize(x, [x > 0], [x]).attained
# False

semialgebraic_maximize(
    x + y + z,
    [x**2 + y**2 + z**2 <= 1],
    [x, y, z],
).value
# sqrt(3)
```

## Candidate pipeline

For multivariate polynomial problems, the exact pipeline combines:

1. safe elimination of equality-constrained variables when an equality is globally linear with a nonvanishing constant coefficient;
2. interior stationary points;
3. equality-constrained KKT/Lagrange systems;
4. algebraic pruning before multiplier systems are constructed: strict-inequality boundaries are excluded from attained KKT sets, scaled duplicate/equality-implied boundaries are removed, inconsistent active ideals prune all supersets, and globally independent equality gradients cap the number of inequality-active generators needed;
5. active-boundary intersections and vertices;
6. explicit Groebner/leading-ideal dimension detection for KKT and singular loci;
7. exact zero-dimensional solving, including the RUR backend where applicable;
8. recursive optimization over positive-dimensional projected KKT loci when equality reduction makes progress, with exact CAD-range certification as the conservative terminal step when no safe coordinate elimination is available;
9. rank-deficient/singular active loci;
10. exact feasibility filtering against the original strict and non-strict constraints;
11. exact comparison of algebraic objective values;
12. complete-CAD certification that no strictly better feasible point exists;
13. cost-controlled CAD function-range certification when the optimum is finite but not attained or finite candidate generation is insufficient.

The singular-locus step is important because ordinary KKT equations can miss extrema when active-constraint gradients lose rank. Fixed-precision numerical comparison is never used to certify which exact algebraic candidate is better; candidate ordering is exact or the path declines conservatively.

## Certification policy

Range certification is controlled by `OptimizationCertificationPolicy` and the corresponding public keyword arguments:

- `certification="auto"` (default) estimates the symbolic image-CAD cost and runs the complete range fallback only when the estimate is at most `range_cost_limit`;
- `certification="complete"` permits the exact range fallback regardless of the estimate;
- `certification="candidate"` disables the full range fallback while still allowing the cheaper CAD query that proves that no strictly better feasible point exists;
- `recursion_limit` bounds recursive optimization of positive-dimensional critical loci.

For example, a three-variable open ball can now return an unattained exact infimum through the range fallback:

```python
semialgebraic_minimize(
    x,
    [x**2 + y**2 + z**2 < 1],
    [x, y, z],
    certification="auto",
).value
# -1
```

`polynomial_locus_dimension(equations, variables)` exposes the exact affine-dimension detector used by the optimizer. It returns `-1` for the empty algebraic locus, `0` for a finite/zero-dimensional locus, a positive integer for positive-dimensional loci, or `None` if the rational Groebner computation cannot be performed conservatively.

## Result semantics

`OptimizationResult` records:

- `value`: the exact minimum/maximum, infimum/supremum;
- `points`: exact optimizing points when attained and available;
- `attained`: whether the bound is achieved;
- `kind`: `"min"` or `"max"`;
- `method`: the successful pipeline/certificate path;
- `certified`: whether global optimality was established exactly;
- `diagnostics`: candidate and certification metadata.

Open sets are distinguished from closed ones. For example, minimizing `x` over `x > 0` gives value `0` with `attained=False`.


## Parameter-stratified optimization and ranges

`semialgebraic_minimize`, `semialgebraic_maximize`, and `function_range` accept `parameters=[...]` together with `return_stratified=True`. The return value is a `ParameterStratifiedResult`. Each branch is guarded by an exact parameter condition.

For optimization, a branch value is `ParametricOptimizationResult`; for ranges it is `ParametricFunctionRangeResult`. These objects deliberately keep the exact first-order relation and its explicit quantifier prefix instead of automatically performing another potentially enormous CAD elimination just to obtain a quantifier-free display. For a minimum with value symbol `t`, the exact relation encodes both

$$
\forall x\;(C(x,p)\Rightarrow f(x,p)\ge t)
$$

and the tightness condition

$$
\forall s>t\;\exists x\;(C(x,p)\land f(x,p)<s).
$$

Thus unattained infima are represented correctly.  In programmatic code, these
relations can be represented directly with semialg's `ForAll` and `Exists`
nodes rather than hand-built quantifier tuples.  For example, the lower-bound
condition has the shape `ForAll(x, Implies(C, f >= t))`, while tightness nests
`Exists` beneath `ForAll`.  `sample_result` on a parametric optimization branch
is representative convenience data only; the quantified relation is the
stratum-wide answer.

## Boolean domains

Bounded DNF expansion is used for disjunctive domains. Each feasible branch is optimized exactly and the branch optima are compared algebraically. Expansion is intentionally bounded; very large Boolean formulas may be declined rather than expanded exponentially.

## Current limits

The optimizer is strongest for low-dimensional polynomial problems whose KKT/active loci are zero-dimensional, safely reducible by equalities, or affordable for exact CAD range certification. Positive-dimensional critical loci are detected explicitly and recursively reduced when possible rather than being silently treated as failed finite solves. Higher-dimensional open/unbounded problems may still be declined when the estimated complete image-CAD cost exceeds the configured policy, and noncompact optimization at infinity is not yet handled by a dedicated asymptotic-critical-point algorithm. General non-polynomial optimization and SOS/SDP certificate search are outside the current exact pipeline.
