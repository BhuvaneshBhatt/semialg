# Future directions

This page collects implementation directions that would substantially improve `semialg`.

## Near-term

- Improve 2D CAD-to-vertical-bounds conversion.
- Add a richer `boundary_description` return object.
- Add `distance_to_region` and `nearest_point` using Lagrange multiplier and critical-point methods.
- Add real-domain utilities:
  - `function_domain`
  - `is_real_valued`
- Add better parameterized `Piecewise` or conditional answers.
- Improve output formatting for intervals, unions, and semialgebraic pieces.

## Medium-term

- Generalize intrinsic integration via parametrizations and metric Jacobians.
- Add critical-point sampling.
- Add roadmap-style connectivity algorithms.
- Add exact algebraic sample-point infrastructure and stronger sign determination.
- Improve radical conversion and exact ordering for algebraic functions.

## Long-term

- Full exact global polynomial optimization.
- SOS/SDP certificate backends.
- Real radical and Positivstellensatz certificate verification.
- RUR/geometric-resolution support for zero-dimensional systems.
- Topological invariants in restricted cases:
  - connected components;
  - Euler characteristic;
  - Betti numbers.
- CAD-backed plotting and discretization tools inspired by symbolic region workflows.

## Region integration roadmap

The most important region-integration target is a robust general pipeline:

```text
CAD decomposition
  -> disjoint cell/stratum representation
  -> vertical or parametric bounds
  -> exact reduced integral pieces
  -> symbolic/numeric evaluation
```

The current implementation provides useful parts of this pipeline for standard shapes and common low-dimensional forms, but the complete general version remains future work.
