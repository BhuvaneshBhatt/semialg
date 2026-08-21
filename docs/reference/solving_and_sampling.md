# Solving and sampling reference

## `solve_semialgebraic`

```text
solve_semialgebraic(
    constraints, variables=None, *, parameters=None, domain="reals",
    count=1, samples=None, sample_mode=None, strategy=None, method="auto",
    variable_order=None, projection_order=None, normalize_domains=True,
    return_formula=False, output=None
)
```

Solves a real semialgebraic system and returns a `SemialgebraicSolution` by default. Depending on `output`/`return_formula`, it can expose formulas, cells, samples, or other structured representations.

`method="auto"` may use specialized exact methods such as rational-univariate representation for finite zero-dimensional equality systems before falling back to broader semialgebraic methods.

## `sample_point` and `sample_points`

Representative sampling is exact by default. Supported workflows include representative/automatic, rational, grid, random, and CAD-cell sampling. Numerical random sampling is explicitly opt-in with `exact=False`.

Returned public samples are checked against the original formula.

## `sign_at` and `sign_vector`

Evaluate polynomial/expression signs at exact points, including algebraic and RUR-backed points. Certified exact paths do not use a hidden floating-point sign guess.

## Instance helpers

`find_instance_formula`, `component_instances`, and related helpers expose witness-oriented workflows over formulas/components.

## Edge cases

- `count=0` requests no samples rather than one implicit sample.
- An empty feasible set is distinct from an unsupported solving strategy.
- String variables follow the shared symbol-resolution rules.

See [Solving](../solving.md) for deeper examples.
