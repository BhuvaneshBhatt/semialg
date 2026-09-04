# How semialg chooses an algorithm

`semialg` has multiple exact backends. Complete CAD is the general correctness engine, but it is deliberately not the first tool used for every problem.

A simplified decision/QE pipeline is:

```text
input formula
    |
normalization
    |
structure-aware presolve
    |-- safe affine equality substitution
    |-- constant-coefficient linear elimination
    |-- variable/block analysis
    |
specialized exact backend?
    |-- quadratic virtual substitution (low-degree QE/witnesses)
    |-- zero-dimensional algebraic solving / RUR
    |
CAD required
    |-- quantifier-aware variable ordering
    |-- reduced / EC-aware path when certifiable
    |-- complete Collins-style fallback when needed
    |
exact reconstruction + simplification
    |
certified result
```

The exact planner differs by public operation, but the design principle is the same: **use structure without weakening the correctness contract**.

## Presolve

Presolve removes structure only when doing so is semantically safe.

For example,

\[
\exists y\;(y=x+1\land y>0)
\]

can substitute \(y=x+1\).

By contrast, blindly solving

\[
a x=1
\]

as \(x=1/a\) is unsafe unless the exceptional parameter case \(a=0\) is handled. The conservative presolver therefore declines transformations that would divide by a coefficient whose nonzeroness is not established.

## Virtual substitution is a backend

Yes. `semialg` contains a **quadratic virtual-substitution backend**.

The high-level decision/solve planner can try virtual substitution for supported low-degree quantified polynomial formulas before invoking CAD. There is also a virtual-substitution witness path. The implementation lives under `semialg.qe.virtual_substitution`.

Virtual substitution is attractive because low-degree quantified formulas can often be eliminated without constructing a full CAD. If the supported fragment does not apply or the pass declines, the planner can continue to another exact backend.

It should therefore be understood as a real backend, not merely a simplification heuristic.

## RUR and zero-dimensional systems

When polynomial equalities define finitely many algebraic solutions, rational univariate representation can encode the solution set through one algebraic parameter and rational coordinate functions.

This is often much more appropriate than decomposing all of real space cylindrically.

## CAD

CAD handles the general semialgebraic case. It projects polynomials to lower-dimensional spaces, decomposes the real line, and recursively lifts stacks so relevant polynomial signs/truth values are invariant on cells.

The implementation includes complete and reduced/EC-aware paths. Reduced methods are used only when their side conditions can be certified; otherwise the solver falls back conservatively to the complete path.

## Variable ordering

CAD cost is highly order-dependent. Automatic ordering is quantifier-aware:

- free variables may be reordered among themselves;
- variables within the same homogeneous quantifier block may be reordered;
- variables never cross an `exists`/`forall` boundary.

Cheap structural scores are preferred for ordinary planning. More expensive projection-based diagnostics can be requested explicitly through the variable-order suggestion API.

## Optimization

A simplified optimization pipeline is:

```text
normalize region/objective
    |
presolve / structural analysis
    |
candidate generation
    |-- stationary/KKT systems
    |-- active boundaries
    |-- lower-dimensional pieces
    |
exact candidate solving
    |
global comparison / certification
    |
parameter stratification when requested
```

Parameterized optimization/range relations try specialized direct reconstruction for supported affine one-dimensional families before requesting a second full QE.

## Integration

Region integration first tries to reduce a semialgebraic region to exact bounds/cells appropriate for symbolic integration. Parameter-dependent univariate fibers can use parameter CAD and algebraic root-function endpoints. More general singular or multidimensional parametric cases may remain unsupported.


## Common structural fast paths

Several high-level geometry operations avoid generic QE when the structure itself gives an exact transformation:

```text
linear existential conjunction
    -> exact Fourier-Motzkin elimination when coefficients have fixed sign

zero-dimensional equality system
    -> RUR / exact algebraic solving

quadratic quantified formula
    -> virtual substitution when the supported fragment applies

nonsingular square affine image
    -> exact inverse substitution

simple affine-polytope boundedness / full dimensionality
    -> linear exact reasoning

general semialgebraic image or projection
    -> existential QE / CAD
```

A fast path is allowed to decline. It must not return a weaker notion of correctness merely because it is cheaper.

## Why the backend can matter

Two mathematically equivalent formulations may have very different costs. A specialized public API often exposes structure that the planner can exploit. See [Which function should I use?](../guides/choosing_an_api.md) and the [Performance guide](../guides/performance.md).
