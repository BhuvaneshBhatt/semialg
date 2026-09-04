# Reference regression suite

The test suite is organized around mathematical contracts rather than only line coverage. The suite exercises the same conclusion through independent algorithms wherever possible and deliberately varies process state, symbol identity, and boundary geometry.

## Deterministic regression tests

Small named regressions cover CAD, QE, algebraic solving, optimization, integration, topology, parameter reconstruction, exact sampling, and known performance and caching failure patterns. A fixed regression should be added whenever a minimized bug example is found.

## Cross-backend differential tests

When two exact backends apply to the same fragment, tests compare their answers semantically rather than requiring identical printed formulas. Current examples include complete CAD versus virtual substitution, specialized parametric reconstruction versus generic QE, and Fourier–Motzkin elimination versus complete QE.

Useful differential pairs include reduced versus complete CAD, RUR versus alternative zero-dimensional solving, and specialized geometry fast paths versus their generic quantified definitions.

## Boundary and discriminant tests

Dedicated cases exercise locations where the topology or algebraic type changes:

- repeated/coalescing roots;
- polynomial degree drops as a parameter vanishes;
- strict versus non-strict shared endpoints;
- zero-width intervals and singleton sections;
- poles and numerator/denominator coincidences;
- singular-rank changes;
- attained versus unattained extrema.

## Property-based tests

A small Hypothesis suite generates low-degree, small-coefficient semialgebraic atoms and checks metamorphic identities such as Boolean identities, contradiction/tautology laws, and variable-renaming invariance. Generated problems are intentionally tiny so failing examples minimize to useful regression cases instead of becoming performance tests.

## Cache and process-state invariance

Exact conclusions must be independent of whether caches are cold, warm, polluted by unrelated computations, or cleared and recomputed. CI can also run the ordinary suite in reproducibly shuffled orders using `SEMIALG_RANDOM_TEST_ORDER` to expose process-global state dependencies.

## Performance architecture tests

Where possible, performance regressions are tested structurally rather than with fragile wall-clock thresholds. Examples include requiring cache-key construction to remain structural, reusing retained CAD connectivity instead of launching pairwise decision solves, and using inverse substitution for nonsingular affine images.

A smaller number of generous wall-clock tests remain marked `slow` for catastrophic regressions.

## Public API contracts

Public result objects are tested for stable fields, formula conversion, witness membership, representative serialization where supported, and explicit distinction among certified emptiness, unsupported fragments, partial results, and errors.

## Distribution checks

The default suite runs on supported Python versions. Additional publishing checks cover the minimum declared SymPy version, reproducible randomized test orders, and the complete slow suite before distributions are built and published.
