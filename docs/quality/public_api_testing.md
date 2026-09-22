# Public API testing

The public API contract is checked at several distinct levels. A direct-call inventory
prevents exported functions from becoming completely untested, but a call alone is not
treated as proof of adequate behavior.

| Level | What it establishes | What it does not establish |
| --- | --- | --- |
| Import contract | Every exported name exists, is callable, documented, and has an inspectable signature | Mathematical correctness |
| Direct-call inventory | Every exported function is invoked somewhere in the recursive test tree | Boundary cases or completeness |
| Behavioral contract | Representative results, failures, witnesses, and result-object fields have assertions | Exhaustive input coverage |
| Property regression | Results agree under specialization, symbol/order changes, cache history, or equivalent formulations | All possible formulas |
| Documentation contract | README, guide examples, links, anchors, signatures, and defaults remain executable/current | Backend completeness |
| Distribution contract | A clean environment can install the wheel and use the public API without the source tree | Every supported Python/platform pair |

Public classes that users normally receive from factories are tested through those
production paths. This includes `CADRegion`, `AffineBoxClip`, and standard-region base
classes. Exception base classes are checked through their inheritance and representative
raising sites rather than by merely instantiating them.

## Adequacy policy

An API is considered adequately covered for release when it has:

1. an import/signature/documentation contract;
2. a direct behavioral assertion or an explicit factory/exception contract;
3. success-path coverage appropriate to its return type;
4. invalid-input or unsupported-path coverage where such a path is part of its contract;
5. cross-checks or mathematical invariants for high-risk exact algorithms.

The direct-call check is therefore a floor, not a coverage score. New APIs should be
accompanied by focused contracts and, for exact algorithms, at least one independent
invariant or differential test. The branch-coverage threshold is a repository-wide
backstop and must not replace API-specific assertions.

## Coverage inventory

`tests/public_api_coverage.toml` is the reviewable source of test ownership and evidence
for every root export. Each entry records the export kind, subsystem owner, risk level,
implemented contract categories, and at least one concrete test file. Factory-returned
and base types record the production symbol through which users encounter them.

Run `python scripts/generate_public_api_coverage.py` after changing the
root API or its test evidence. The contract suite regenerates the file in memory and
requires an exact match, rejects missing or stale exports, validates all test paths, and
checks that function evidence contains a direct call. Known gaps are explicit data: for
example, `UnsupportedFragmentError` is currently tested as a public constructible type
but has no production raising site.

Every public function is covered beyond a single-file smoke contract. Function entries
record `evidence_file_count`, retain two concrete evidence paths for review, and use
`coverage_status = "expanded"`. Contract tests recompute the direct-call inventory and
reject any function whose evidence falls below two distinct test files. Focused risk
suites add outcome pairs, exact witness checks, cross-API invariants, boundary semantics,
and validation of structured results.

## API-family case matrices

`tests/test_public_api_family_contracts.py` exercises related operations against shared,
small exact cases. Each parameter combination is collected separately with an identifiable
case name, so failures locate the mathematical boundary being tested.

| Family | Cases and assertions |
| --- | --- |
| Decision | Empty/universal formulas, strict conflicts, irrational solutions, and punctures; Boolean and structured returns, witness validity, emptiness, and negated tautology |
| Topology | Empty sets, points, open/closed/half-open intervals, and punctured lines; exact closure, interior, boundary, dimension, and component counts |
| Optimization | All four endpoint conventions crossed with increasing, decreasing, and constant objectives; exact bounds, attainment, certified results, feasible optimizers, and optimizer sets |
| Root classification | Zero/nonzero constants, simple/repeated/irrational/nonreal roots, and linear coefficient degree drops; counts, sorted multiplicity patterns, and specialization consistency |
| Connectivity | Excluded points and lines versus included connecting points; decomposition, samples, and univariate roadmaps |
| Integration | Four monomials over open/closed intervals, disconnected unions, empty sets, and ambient-measure singletons; scalar and structured results against elementary antiderivatives |
| Function predicates | Identity, affine, square, constant, and absolute-value maps; injectivity, surjectivity, bijectivity, continuity, and smoothness |
| Matrix predicates | All sign combinations for two diagonal entries; positive definiteness, semidefiniteness, and rank |

Expected sets and values are explicit mathematical oracles; agreement between two package
APIs is an additional check, not the sole correctness criterion. These bounded cases run
in the non-slow contract tier. They do not imply exhaustive API coverage: unsupported-input
contracts, broader multivariate families, and property-based tests remain separate concerns.

## Generated invariants and differential checks

The non-slow property tier includes two complementary public-API suites:

- `tests/properties/test_public_api_metamorphic.py` checks Boolean set algebra,
  invertible affine images and preimages, negative-scale endpoint reversal, integral
  additivity and change of variables, and extrema under objective scaling.
- `tests/properties/test_public_api_differential.py` checks root counts and multiplicities
  against generated factors and the quadratic discriminant, solver/CAD formulas against
  exact intervals, affine feasibility against complete quantifier elimination, and
  decision results across cold, warm, and unrelated cache histories.

These tests generate small integer coefficients, bounded degrees, and at most two
variables. Each property has an explicit example budget; important degenerate cases also
use `@example`, so they run regardless of the random seed. Hypothesis can shrink generated
failures to smaller counterexamples. Timing deadlines are disabled because symbolic
evaluation time varies, while the input bounds limit the scope of each run.

The oracles include explicit SymPy sets, elementary antiderivatives, factor multiplicity
counts, and closed-form extrema. Related public APIs can share implementation paths, so
agreement alone is insufficient. The affine feasibility differential check additionally
disables presolve and the variety shortcut on its complete-CAD path, and rejects an
unknown truth value rather than treating it as false.

To reproduce the generated API checks with a fixed seed:

```sh
python -m pytest -q tests/properties/test_public_api_metamorphic.py tests/properties/test_public_api_differential.py --hypothesis-seed=20260912 --no-cov
```

Use another seed to explore additional inputs. A passing seeded run establishes the
tested invariants for its examples; it does not establish exhaustive input or branch
coverage, and the ordinary test count counts properties rather than generated examples.

## Failure-path contracts

`tests/test_public_api_failure_contracts.py` checks invalid inputs and unsupported modes
separately from valid negative mathematical results. It covers decision domains, integral
bounds and modes, optimization branch limits, matrix shape and symmetry, sampling
controls, nonpolynomial root inputs, and malformed image/preimage coordinates. Duplicate
target coordinates are rejected after symbol resolution, including string names.

Paired direct and structured calls check that failures do not turn into false decisions,
zero integrals, or empty geometry. Recovery cases ensure that a rejected request does
not prevent a subsequent valid call. CAD tests cover enforced preprocessing limits and
the explicit rejection of unimplemented time/cell budgets, including warm-cache calls.

An explicit internal sampler hook checks exception boundaries: an expected unsupported
sampling case returns no witness, but unexpected assertion and runtime errors propagate.
Candidate witnesses are checked for completeness and membership. Tests require specific exception classes and meaningful
message fragments rather than accepting arbitrary exceptions. The public
`UnsupportedFragmentError` production-raising-site gap remains recorded in the inventory;
these tests do not create an artificial raising site merely to satisfy coverage.

## Factory-returned types and result objects

Public classes should be tested the way application code receives them.  The permanent
production-path contracts are:

- `CADRegion` through `as_cad_region`;
- `AffineBoxClip` through `clip_affine_subspace_to_box`;
- `Geometry` and `StandardRegion` through a concrete `Interval`;
- `ResourceLimitError` and its public base classes through a resource-limited `cad` call.

`tests/test_public_type_contracts.py` asserts useful fields and behavior after those
operations; it does not merely instantiate the classes. `UnsupportedFragmentError`
remains the explicit exception: the package currently exposes the type but has no
production raising site, so its construction/inheritance contract and the inventory gap
remain visible rather than fabricating a failure path.

Structured result and certificate classes follow the same rule.
`tests/test_factory_result_contracts.py` obtains representative CAD, parameter, root,
optimization, function-range, integration, and convexity objects from their public
operations.  The tests then validate mathematical fields, certification/attainment
metadata, serialization where promised, and frozen-result behavior.  Direct dataclass
construction can still be useful for isolated low-level tests, but it is not considered
sufficient integration coverage for a public operation/result pair.

## Measured coverage policy

Coverage is a backstop for behavioral tests, not a reason to add execution-only cases.
The repository keeps the current combined statement/branch gate at **74.5%** until a
fresh complete non-slow run justifies a higher threshold.  Higher checkpoints such as
76%, 78%, and 80% are adopted only after meaningful decision,
parameter/root, topology/optimization/integration, or remaining thin-family contracts
raise the measured coverage naturally.

`python scripts/check_coverage.py build/coverage/coverage.json` reports statement,
branch, and combined percentages independently.  This distinction matters because
coverage.py's repository `fail_under` value is a combined statement/branch percentage
when branch measurement is enabled; it should not be described as a branch-only number.
The checker also requires every `src/semialg/*.py` source file to appear in the report
and assigns every file to exactly one subsystem.  Missing files, overlapping ownership,
and unassigned nested subsystems are errors.

See [Coverage measurements](coverage_measurements.md) for the subsystem policy and the
rules for increasing floors.

## Semantic coherence contracts

`tests/test_decision_semantic_coherence.py` keeps a small, non-slow set of identities between the everyday decision APIs. These are deliberately redundant at the **mathematical** level but not at the implementation-contract level:

- `implies(A, B)` agrees with unsatisfiability of `A & ~B`;
- `equivalent(A, B)` agrees with implication in both directions;
- variable renaming preserves satisfiability;
- multiplication of a polynomial inequality by a known positive constant preserves its feasible set.

These checks complement explicit expected-value cases. They are especially useful after planner, normalization, or fast-path refactors because they exercise the same mathematical statement through different public entry points. They do not replace independent mathematical oracles: if two APIs share an incorrect backend, agreement alone cannot detect the error.

## Reader-facing example style

Documentation examples show the expression a reader would evaluate followed by its expected result. For example:

```python
is_satisfiable(x**2 <= 1, [x])
# True
```

Tests use assertions; README and guide examples normally do not repeat the same expectation as both an assertion and an output comment. `tests/test_documentation_presentation_contracts.py` protects this convention and also checks that the computation-flow diagrams remain present and linked.


## Root-level adequacy gate

File-level call coverage is necessary but not sufficient. The root API also has a named-contract adequacy gate: every exported function must be exercised by at least two independently named tests, while critical decision/solver APIs require at least three. The current function-by-function audit is recorded in [Root API adequacy audit](root-api-adequacy.md). The audit deliberately labels functions with exactly two contracts as *adequate-minimum* rather than pretending that a test count proves mathematical completeness.
