# Canonical representations and input conventions

`semialg` operates on exact polynomial and semialgebraic objects, usually expressed with SymPy.

## Formulas and regions

A semialgebraic region is normally represented by a Boolean combination of polynomial relations:

```python
import sympy as sp

x, y = sp.symbols("x y", real=True)
region = (x**2 + y**2 <= 1) & (y >= 0)
```

Use `Eq(f, 0)` for polynomial equalities and SymPy inequalities for strict or non-strict sign conditions. Boolean combinations use `And`, `Or`, and `Not` (or their operator equivalents where SymPy supports them).

Standard-region objects such as `IntervalRegion`, `BoxRegion`, and `BallRegion` are convenience representations with explicit geometric invariants. They can be converted or lowered to semialgebraic formulas as required by downstream algorithms.

## Exact coefficients

Prefer exact coefficients:

```python
sp.Rational(1, 3)
sp.sqrt(2)
```

over binary floating-point literals when exact semantics are intended:

```python
0.333333333333
```

A float is already an inexact input. `semialg` does not silently recover the exact rational number the user may have intended.

## Variables and parameters

Variables are the coordinates being solved, projected, optimized, or integrated over. Parameters are free symbolic quantities whose values determine a family of problems.

For example, in

\[
x^2\le a,
\]

`x` may be an integration variable while `a` is a parameter.

When an API accepts explicit `variables=[...]`, symbols not listed there may remain parameters rather than being implicitly eliminated. Use the API's `parameters=` argument when available to make the distinction explicit.

## Symbol identity matters

Two SymPy symbols with the same printed name but different assumptions are distinct mathematical objects:

```python
x_real = sp.Symbol("x", real=True)
x_pos = sp.Symbol("x", positive=True)
```

`semialg` preserves Symbol identity in caches and internal normalization. String keys are resolved contextually only in APIs that explicitly support them. Prefer exact Symbol keys when ambiguity is possible.

## Quantifiers

For programmatic quantified formulas, prefer `Exists` and `ForAll`:

```python
from semialg import Exists, ForAll

formula = ForAll(x, Exists(y, sp.Eq(x + y, 0)))
```

Textual interfaces remain available, but programmatic nodes preserve symbol identity and are easier to compose safely.

## Algebraic numbers and roots

Exact real algebraic values may be represented by SymPy algebraic expressions, `RootOf`/`CRootOf`, or `semialg` root-function/RUR structures. These are not approximations.

When roots depend on parameters, `AlgebraicRootFunction`-style representations preserve which ordered root branch is meant on a certified parameter cell.

## Solution and result representations

A Boolean formula, a point witness, a CAD decomposition, an optimization result, and a parameter-stratified result answer different questions. Avoid flattening a structured result to a numerical approximation if it will feed another exact computation.

See [Understanding result objects](result_objects.md).

## Internal symbols

Algorithms introduce collision-free internal `Dummy` symbols for graph variables, multipliers, deformation parameters, and related auxiliaries. User symbols that happen to have similar printed names therefore do not alias internal variables.

## Representation versus simplification

Equivalent semialgebraic sets can have many syntactic formulas. The simplification layer computes a **stable canonical form for the supported polynomial fragment**: polynomial relation residuals are primitive and deterministically oriented, repeated zero-set multiplicities are removed, compatible scalar bounds are merged, Boolean branches are normalized, and guarded CAD implication checks remove provable semantic redundancy. Reapplying the simplifier is idempotent on this fragment.

This is not presented as a globally minimum Boolean formula: globally minimizing arbitrary semialgebraic formulas would require a cost model and can be computationally much harder than establishing equivalence. Callers that need mathematical equality rather than stable presentation should still use semantic predicates such as `equivalent` or `is_equal`.
