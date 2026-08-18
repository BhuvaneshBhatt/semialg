# Implementation notes

## Conservative exactness policy

`semialg` prefers exact positive results or explicit `NotImplementedError` over unreliable symbolic guesses. Many functions have fast special paths, CAD/QE-backed checks, and then conservative failure.

## Formula representation

The package currently uses SymPy Boolean formulas and relational expressions as the primary user-facing representation. Most public APIs accept either a single formula or a list of constraints.

## CAD/QE as a semantic oracle

Many simplification routines are implemented by asking semantic questions:

- Is this constraint redundant?
- Does this assumption imply this branch condition?
- Is this region empty?
- Are these two formulas equivalent?

This semantic approach is more robust than purely syntactic rewriting, but it can be more expensive.

## Function range as semialgebraic image

`function_range` treats range computation as image projection. For an expression `f` and domain `C`, it constructs a semialgebraic graph relation and eliminates the original variables to obtain a condition on the value symbol.

Graph encodings include:

```text
Abs(u):      (t >= 0) and (t = u or t = -u)
sqrt(u):     (t >= 0) and t**2 = u and u >= 0
Max(u, v):   (t = u and u >= v) or (t = v and v >= u)
Min(u, v):   (t = u and u <= v) or (t = v and v <= u)
Piecewise:   disjunction of branch graph relations
```

## Region integration strategy

Region integration is layered:

```text
recognized shapes
  -> exact formulas or vertical bounds
  -> reduced integral pieces
  -> symbolic / numeric / auto evaluation
  -> intrinsic-measure special cases
```

## Why complete region integration is difficult

CAD can decompose a semialgebraic set, but integration still requires:

- converting cells to usable bounds;
- handling algebraic boundary functions;
- integrating algebraic expressions symbolically;
- computing lower-dimensional metric Jacobians;
- avoiding double counting in mixed-dimensional unions.

The current implementation is therefore a layered exact-first framework, not a complete arbitrary CAD-strata integration engine.

## Testing philosophy

Tests should prefer semantic equivalence and exact mathematical invariants over brittle string comparisons.

Examples:

```python
assert equivalent(result, expected, variables)
assert is_satisfiable(formula, variables)
assert integrate_over_region(...) == expected
```
