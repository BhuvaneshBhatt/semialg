# Applications reference

`semialg.applications` contains domain-oriented workflows that compose the certified core APIs. They do not duplicate CAD, QE, optimization, range, measure, or integration primitives.

## Robust design and tolerance analysis

### `robust_parameter_analysis(constraints, operating_variables, parameters, *, operating_domain=None)`

Returns a `RobustParameterResult` containing exact parameter conditions for existential feasibility, universal robustness, and existence of a violating operating point.

### `robust_parameter_region(..., quantifier="forall", operating_domain=None)`

Returns only the requested universal or existential parameter condition.

## Certified symbolic-math validation

### `validate_identity(lhs, rhs, variables=None, *, assumptions=None)`

Certifies an expression identity globally or under semialgebraic assumptions.

### `validate_formula_equivalence(original, proposed, variables=None)`

Certifies equality of two real semialgebraic sets/formulas.

### `validate_range(expression, proposed_range, variables, *, constraints=None, value_symbol="t")`

Computes the exact range and certifies equivalence with a proposed range formula.

## Numerical-optimizer benchmark oracle

### `exact_optimization_benchmark(objective, constraints, variables, *, kind="min", certification="auto")`

Returns an `OptimizationBenchmark` retaining the exact `OptimizationResult` and global certification status.

### `validate_numeric_optimization(benchmark, numeric_value, *, atol=1e-8)`

Checks a reported numerical objective value against the exact reference optimum.

## Polynomial control stability

### `polynomial_stability_analysis(polynomial, variable, parameters=None)`

Returns a `PolynomialStabilityResult` with the strict Hurwitz-stability condition, Hurwitz matrix, and principal determinants. The characteristic polynomial must have positive degree and real coefficients.

### `polynomial_stability_region(polynomial, variable, parameters=None)`

Convenience function returning only the exact stability condition.

## Polynomial safety invariants

### `verify_polynomial_invariant(invariant, transition, variables, *, initial_condition=None, unsafe_condition=None, domain=None)`

Returns an `InvariantVerificationResult` containing exact initiation, inductiveness, and safety checks plus available counterexamples. Each state variable must have one polynomial update; undeclared symbolic parameters are rejected.

## Polynomial response surfaces

### `analyze_response_surface(model, variables, *, domain=None, thresholds=(), certification="auto")`

Returns a `ResponseSurfaceResult` containing exact minima/maxima, function range, gradient, stationary condition, and requested threshold superlevel sets. The current workflow expects exact numeric/algebraic model coefficients and rejects undeclared symbolic parameters.

## Core/application boundary

The following remain core APIs and are intentionally not re-exported as applications:

- `function_range`
- `semialgebraic_measure`
- `integrate_over_region`
- `semialgebraic_minimize`, `semialgebraic_maximize`
- CAD and QE functions

Applications depend on these primitives; moving them into `semialg.applications` would create duplicate public homes and weaken the package architecture.


## Lyapunov verification

### `verify_lyapunov_function(function, dynamics, variables, *, domain=None, equilibrium=None, derivative_strict=True)`

Returns a `LyapunovVerificationResult`. The result records equilibrium/domain validity, positive definiteness, the Lie derivative, the derivative-sign certificate, and available counterexamples. `dynamics` may be a mapping keyed by state variables or a sequence aligned with `variables`.

## Barrier certificates

### `verify_barrier_certificate(barrier, dynamics, variables, *, initial_condition, unsafe_condition, domain=None, derivative_strict=False)`

Returns a `BarrierVerificationResult` proving initial inclusion in `B <= 0`, unsafe-set separation by `B > 0`, and a nonpositive (or strict negative) Lie derivative on `B = 0`.

## Sensitivity and monotonicity

### `analyze_polynomial_sensitivity(model, variables, *, domain=None)`

Returns a `SensitivityAnalysisResult`. Each variable maps to a `SensitivityDirectionResult` containing the exact derivative, an exact `FunctionRangeResult` for that derivative, certified sign flags, and a coordinate-wise monotonicity classification.

## Constraint analysis

### `analyze_constraint_redundancy(constraints, variables=None)`

Returns a `ConstraintRedundancyResult` with zero-based redundant and essential constraint indices and nonredundancy witnesses when available.

### `diagnose_feasible_set(constraints, variables=None, *, find_conflict=True)`

Returns a `FeasibleSetDiagnosticResult`. Feasible systems include a witness and redundancy information. Infeasible systems can include an irreducible conflict index set. Set `find_conflict=False` to skip the repeated satisfiability checks used for conflict extraction.

## Polynomial model comparison

### `compare_polynomial_models(first, second, variables, *, domain=None, certification="auto")`

Returns a `PolynomialModelComparisonResult` with:

- the exact signed difference `first - second`;
- exact minimum and maximum difference `OptimizationResult` objects;
- the exact maximum squared discrepancy and its nonnegative square root;
- certified `first_le_second`, `first_ge_second`, and `equivalent_on_domain` flags;
- available counterexamples to failed dominance/equivalence claims.

Both models must be polynomial in the declared predictor variables. Undeclared symbolic coefficients are rejected.

## Parameter regimes

### `analyze_parameter_regimes(constraints, variables, parameters)`

Returns a `ParameterRegimeResult` whose `stratified_result` partitions parameter space into exact `True`/`False` solvability regimes. `select(assignments)` returns the certified solvability value at a fully specified parameter point.

### `analyze_root_count_regimes(polynomial, variable, parameters=None)`

Returns a `ParameterRegimeResult` whose branch values are the exact number of distinct real roots. The polynomial must be univariate in `variable` after treating the declared parameters as coefficients.

`ParameterRegimeResult` exposes `regime_count`, `certified`, and `select()`, while retaining the underlying `ParameterStratifiedResult` for branch conditions and metadata.

## Polynomial probability

### `polynomial_probability(event, variables, *, support=None, density=1, bounds=None)`

Returns a `PolynomialProbabilityResult`. The density must be polynomial and is certified nonnegative over the effective support. The result contains:

- `normalizing_mass` — exact integral of the density over the support;
- `event_mass` — exact integral over support intersected with the event;
- `probability` — exact ratio `event_mass / normalizing_mass`;
- `density_nonnegative` and `certified` flags.

The normalizing mass must be finite and strictly positive.

### `geometric_probability(event, variables, *, support, bounds=None)`

Uniform-density convenience wrapper around `polynomial_probability(..., density=1)`.
