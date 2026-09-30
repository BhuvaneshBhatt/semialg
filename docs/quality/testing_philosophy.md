# Testing philosophy

Semialg tests mathematical contracts instead of treating line coverage as evidence of correctness. The permanent suite combines direct public-API behavior, differential and metamorphic identities, exact reference regressions, failure/resource semantics, documentation examples, and installed-distribution checks.

## Public API contract matrix

The live root API, `docs/reference/primary_api_manifest.toml`, and `tests/public_api_coverage.toml` are required to describe the same exported surface. Every public function has independent direct-call evidence; structured types identify their production path; exceptions identify construction or raising sites. The manifests are generated from the live package and maintained documentation, then checked for reproducibility.

Documentation adequacy is stricter than name coverage. Risk-adjusted checks require purpose and return semantics for every root function, exactness and limitations at the owning family level, parameter semantics for high-risk operations, executable examples for high-risk APIs, and algorithm/backend discussion for critical APIs.

Generated details remain available in [root API documentation adequacy](root-api-documentation-adequacy.md), [root API test adequacy](root-api-test-adequacy.md), and the [multidimensional root API audit](root-api-multidimensional-testing.md), but those reports are evidence artifacts instead of separate testing philosophies.

## Reference and regression corpus

Reference regressions encode exact expected mathematical behavior and permanent bug reproductions. See [Reference regression suite](reference_regression_suite.md). New regressions should assert semantic outcomes, certification status, or explicit incompleteness instead of implementation-specific intermediate expressions unless the representation itself is part of the contract.

## Property and metamorphic testing

Generated tests exercise invariants with independent mathematical expectations: Boolean set algebra, affine image/preimage equivalence, change of variables, optimization invariants, and local/approach-geometry transformations. For local geometry in particular, coordinate renaming and invertible affine changes must preserve closure semantics, while correlated images must reject Cartesian-product over-approximations that violate the exact joint relation.

## Certification and soundness

A symbolic-looking result is not automatically certified. Tests distinguish exact/certified conclusions from candidates, diagnostics, numerical presentation, and explicit unknown/incomplete outcomes. Certificate replay tests verify proof objects independently of the search that produced them. See [Decision contracts](decision_contracts.md) and [Certified decisions](certified_decisions.md).

## Documentation as executable behavior

The worked-example gallery has a one-to-one Markdown/Python pairing. CI parses every script, requires every script to have a documented companion, and checks that the fenced Python block is byte-for-byte synchronized with its executable source. Core gallery paths are executed as smoke contracts. Documentation registries and generated adequacy reports must also be reproducible from a source distribution.

## Downstream compatibility

Semialg does not depend on downstream packages. `tests/test_asymptotic_downstream_contract.py` is an optional compatibility contract: release validation can install asymptotic and select the `downstream` marker to exercise the exact Semialg surface used for multivariate limit support. Without asymptotic installed the test skips, preserving the dependency direction.

A downstream compatibility failure should be resolved by deciding which package owns the broken contract; Semialg must not add downstream-specific workarounds that weaken its own mathematical semantics.
