# Future directions

This page lists areas where the current algorithms can be extended or made more
efficient.

## Engineering

- Migrate remaining public wrappers to the shared symbol resolver.
- Narrow broad fallback handlers into explicit strategy, certification,
  arithmetic, and resource-limit failures.
- Reuse compatible sub-CADs and formula-truth caches within computation
  contexts while preserving bounded memory use.
- Improve formatting for intervals, unions, algebraic root functions, and
  conditional results.
- Expand API-level tests for less frequently used result classes and region
  constructors.

## CAD performance

- Calibrate pilot-lifting and arithmetic-growth cost weights against a larger
  benchmark corpus.
- Propagate more logically necessary equational constraints through projection
  chains.
- Add stronger partial-CAD and truth-table-invariant CAD techniques.
- Reuse compatible sub-CADs across related decision and certification calls.
- Improve projection pruning and formula decomposition before full CAD.

## Parameter-dependent algorithms

- Produce parameter-stratified integration results when the integrand or domain
  changes topology across parameter cells.
- Add on-demand elimination of exact first-order optimization and range
  relations into quantifier-free formulas when affordable.
- Improve normalization and merging of equivalent guarded results.

## Geometry and integration

- Support regular CAD manifold strata that require alternative local charts.
- Develop singular-stratum handling beyond isolated or lower-dimensional
  singularities.
- Add richer topology and roadmap algorithms together with selected topological
  invariants.
- Extend exact algebraic-function integration beyond direct SymPy evaluation.

## Optimization and certificates

- Improve high-dimensional candidate generation and recursive boundary
  optimization.
- Add optional SOS/SDP certificate backends with exact certificate
  verification.
- Add real-radical and Positivstellensatz certificate verification where
  useful.
