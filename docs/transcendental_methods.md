# Transcendental solving

The `semialg.solve.transcendental` package handles selected real problems that
contain functions outside the first-order theory of real closed fields. Its
results distinguish exact reductions from bounded or heuristic search; symbolic
output alone is never treated as a completeness certificate.

## Problem representation

`TransProblemState` records the formula, free and quantified variables,
parameters, variable domains, quantifier blocks, variable order, and solver
metadata. `build_trans_state` normalizes these fields and also accepts a leading
`Exists`/`ForAll` prefix from semialg's quantified-expression layer.

Quantifier-aware preprocessing records a `QuantifierDispatchPlan` describing the
leading block and the structural method that can be attempted before numerical
search.

## Structural preprocessing

`preprocess_transcendental_problem` performs Piecewise simplification and can
replace supported function families with fresh auxiliary variables plus exact
defining equations. Supported detectors cover trigonometric, hyperbolic,
exponential, inverse, ProductLog/Lambert-W-style, statistical, argument, and
other registered special-function families.

Auxiliary replacement is restricted to equations when that restriction is
required to preserve the intended quantified semantics. Every replacement step
records the source expressions, new variables, and equations needed for later
reconstruction.

## Univariate roots and inequalities

For supported real univariate expressions, the solver can construct certified
sign-change brackets and refine them by bisection. `CertifiedIntervalRoot`
records the isolating interval and certification data. Inequality decomposition
uses the ordered root brackets to determine truth on the intervening intervals.

Periodic functions use SymPy's periodicity information when available.
Representative solutions can be lifted from a fundamental domain into explicit
periodic formulas using semialg quantifiers rather than encoding quantification
implicitly through `Mod` or `ImageSet`.

## Quantifier elimination

Limited genuinely univariate transcendental quantifier blocks can be reduced by
`eliminate_leading_real_quantifier_block`. The result records whether the
reduction is exact and which variables remain for another backend.

Problems that cannot be reduced exactly are not silently promoted to complete
answers. `ResultSemantics` describes how an incomplete formula relates to the
true solution set, including exact, subset, superset, window-scoped, and bounded
approximation semantics.

## Multivariate fallback

Multivariate systems can use structural family rewrites, `nonlinsolve`, and
seeded `nsolve` as bounded witness-search tools. Candidate points are checked by
residual/certification logic before they are reported. Failure to find a point
in a bounded search window does not prove global unsatisfiability.

## Relevant APIs

- `build_trans_state`
- `preprocess_transcendental_problem`
- `replace_function_families`
- `certified_sign_change_isolation`
- `certified_interval_decomposition`
- `reconstruct_periodic_solution_from_representatives`
- `reconstruct_periodic_intervals_from_fundamental_domain`
- `eliminate_leading_real_quantifier_block`
- `orchestrate_transcendental_system_search`

See [Errors and failure modes](guides/errors_and_failure_modes.md) for the
distinction between unsupported exact reductions and bounded numerical search.
