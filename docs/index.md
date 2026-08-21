# semialg documentation

`semialg` provides exact symbolic algorithms for real polynomial and semialgebraic problems. The documentation is organized by how you use the package rather than by implementation history.

## Learn semialg

1. **[Getting started](getting_started.md)** — installation, formulas, variables, decisions, solving, optimization, regions, and integration.
2. **[Exactness and certification](concepts/exactness_and_certification.md)** — what semialg means by exact, certified, candidate, and numerical.
3. **[Progressive demo notebook](../notebooks/semialg_demo.ipynb)** — executable examples that build from decision problems to optimization and integration.

## Task guides

- [Performance](guides/performance.md) — variable order, equational constraints, caching, and certification cost.
- [Errors and failure modes](guides/errors_and_failure_modes.md) — distinguish invalid input, infeasibility, unsupported cases, cost limits, and failed certification.
- [Symbol handling](guides/symbol_handling.md) — SymPy symbol identity, string resolution, parameters, and assumptions.
- [Region invariants](guides/region_invariants.md) — constructor requirements and geometric validation.

## API reference

Start with the **[API overview](api_overview.md)**, then use the family references:

- [Decision and QE](reference/decision_and_qe.md)
- [Solving and sampling](reference/solving_and_sampling.md)
- [CAD and structured geometry](reference/cad.md)
- [Optimization and ranges](reference/optimization_and_range.md)
- [Regions](reference/regions.md)
- [Integration and moments](reference/integration_and_moments.md)
- [Algebraic roots and exact solving](reference/algebraic.md)
- [Parameters and conditional results](reference/parameters.md)

## Explanations and detailed topics

- [CAD and QE](cad_qe.md)
- [Decision procedures](decision_procedures.md)
- [Solving](solving.md)
- [Optimization](optimization.md)
- [Function ranges](function_range.md)
- [Region operations](region_operations.md)
- [Region integration](region_integration.md)
- [Moments](moments.md)
- [Root classification](root_classification.md)
- [Parameter stratification](parameter_stratification.md)
- [Symbolic simplification](symbolic_simplification.md)
- [Quantified expressions](quantified_expressions.md)
- [Implicit geometry](implicit_utilities.md)

## Scope, quality, and internals

- [Feature matrix](feature_matrix.md)
- [Limitations](limitations.md)
- [Architecture](architecture/design.md)
- [Code quality](architecture/code_quality.md)
- [Robustness](quality/robustness.md)
- [Decision contracts](quality/decision_contracts.md)
- [Reference regression suite](quality/reference_regression_suite.md)
- [Future directions](future_directions.md)
- [Changelog](../CHANGELOG.md)

## Applications

The application layer groups thin certified workflows into three practical families:

- **Control and verification applications:** Hurwitz stability, discrete invariants, Lyapunov functions, and barrier certificates.
- **Design and model analysis:** robust parameter analysis, response surfaces, polynomial model comparison, parameter regimes, sensitivity, constraint redundancy, and feasibility diagnostics.
- **Validation, testing, and probability:** symbolic-math validation, exact numerical-optimizer benchmarks, and exact polynomial/geometric probability.


- [Applied workflows](applications.md)
- [Applications reference](reference/applications.md)
