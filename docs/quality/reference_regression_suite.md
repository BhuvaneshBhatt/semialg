# Reference regression suite

The test suite is organized around mathematical contracts instead of only line coverage. The suite exercises the same conclusion through independent algorithms wherever possible and varies process state, symbol identity, and boundary geometry.

## Deterministic regression tests

Small named regressions cover CAD, QE, algebraic solving, optimization, integration, topology, parameter reconstruction, exact sampling, and known performance and caching failure patterns. A fixed regression should be added whenever a minimized bug example is found.

## Cross-backend differential tests

When two exact backends apply to the same fragment, tests compare their answers semantically instead of requiring identical printed formulas. Current examples include complete CAD versus virtual substitution, specialized parametric reconstruction versus generic QE, and Fourier–Motzkin elimination versus complete QE.

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

A small Hypothesis suite generates low-degree, small-coefficient semialgebraic atoms and checks metamorphic identities such as Boolean identities, contradiction/tautology laws, and variable-renaming invariance. Generated problems are tiny so failing examples minimize to useful regression cases instead of becoming performance tests.

## Cache and process-state invariance

Exact conclusions must be independent of whether caches are cold, warm, polluted by unrelated computations, or cleared and recomputed. CI can also run the ordinary suite in reproducibly shuffled orders using `SEMIALG_RANDOM_TEST_ORDER` to expose process-global state dependencies.

## Performance architecture tests

Where possible, performance regressions are tested structurally instead of with fragile wall-clock thresholds. Examples include requiring cache-key construction to remain structural, reusing retained CAD connectivity instead of launching pairwise decision solves, and using inverse substitution for nonsingular affine images.

A smaller number of generous wall-clock tests check for catastrophic regressions. They run in the default suite; the `slow` and `performance` markers allow focused profiling.

## Public API contracts

Public result objects are tested for stable fields, formula conversion, witness membership, representative serialization where supported, and explicit distinction among certified emptiness, unsupported fragments, partial results, and errors.

## Distribution checks

The default suite includes the computational and performance regression cases and runs on supported Python versions. Large `stress` workloads require explicit selection with `pytest -m stress tests/stress`. Publication depends on the ordinary test matrix, quality checks and optional SDP checks. The dependency matrix covers minimum and latest effective arithmetic requirements. A separate weekly run profiles the computational subset and installs Singular for the independent algebra corpus; those extra measurements do not gate ordinary publication. Reproducible randomized orders can also be selected for local runs.

## Package lifecycle and differential corpus

The lifecycle regression layer adds deterministic generated invariants in `tests/test_package_lifecycle.py`. It differentially compares certified modular Groebner bases with direct exact Groebner bases, generates families of primary decompositions with known linear associated primes, checks invariance under multiplication by rational units, replays every resulting certificate, and repeatedly exercises cache-heavy workloads before unified teardown.

`scripts/benchmark_package.py` reports JSON records for cold/warm GTZ, direct-versus-modular Groebner computation, CAD feasibility, and cache state before/after teardown. `scripts/differential_singular.py` is an optional external smoke corpus: when Singular is installed it compares dimension, associated primes, primary components, radical ideals, multiplicities and reconstruction against `primdecGTZ`, while semialg independently replays its own exact certificates. The external script exits successfully with an explicit skip message when Singular is unavailable.

## Shared finite quotient arithmetic

The GTZ driver builds finite algebras over QQ(parameters) with `algroots.quotient.QuotientAlgebra`. Standard monomials, multiplication characteristic polynomials, trace rank and nilradicals use that backend. semialg retains localization/contraction, component assembly, proof records and independent certificate replay. The replay reconstructs the quotient from the stored source and checks each obligation; a backend result alone does not establish a GTZ certificate. Stored monomial coordinates follow the backend's total-degree ordering. Zero-algebra certificates are checked separately.
