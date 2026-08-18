# semialg documentation

`semialg` provides exact symbolic algorithms for real semialgebraic geometry and CAD-backed symbolic reasoning.

Core themes:

1. Decision procedures over the real numbers
2. Semialgebraic solving and sampling
3. Region operations and integration
4. Symbolic optimization and range computation
5. Assumption-aware simplification.

## Start here

If you want to:

- check whether constraints are feasible, see [`decision_procedures.md`](decision_procedures.md) and [`user_guide/decision_sampling_signs.md`](user_guide/decision_sampling_signs.md);
- solve systems of equations and inequalities, see [`solving.md`](solving.md);
- classify parameterized roots, see [`root_classification.md`](root_classification.md);
- compute ranges or extrema, see [`function_range.md`](function_range.md) and [`optimization.md`](optimization.md);
- manipulate regions, see [`region_operations.md`](region_operations.md);
- integrate over regions, see [`region_integration.md`](region_integration.md);
- compute moments, centroids, or covariance matrices, see [`moments.md`](moments.md);
- simplify formulas or expressions under assumptions, see [`symbolic_simplification.md`](symbolic_simplification.md);
- understand public decision-layer quality contracts, see [`quality/decision_contracts.md`](quality/decision_contracts.md);
- understand internal design choices, see [`implementation_notes.md`](implementation_notes.md);
- see planned work, see [`future_directions.md`](future_directions.md).

## Philosophy

`semialg` follows an exact-first policy. When it can prove a result exactly, it returns an exact symbolic answer. When a case is outside the current supported algorithms, it should fail conservatively rather than produce an unreliable expression.

The package uses CAD/QE as a semantic oracle for many operations, but it also includes specialized faster paths for common low-dimensional and standard-region problems.
