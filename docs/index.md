# semialg documentation

`semialg` is an exact-first Python package for computation with polynomial and semialgebraic sets over the real numbers. It combines complete real quantifier elimination with specialized exact backends for finite algebraic solving, low-degree virtual substitution, optimization, integration, topology, and geometry.

A typical input is a SymPy Boolean formula such as

```python
import sympy as sp
from semialg import (
    function_range,
    is_satisfiable,
    semialgebraic_measure,
    semialgebraic_minimize,
)

x, y = sp.symbols("x y", real=True)
disk = x**2 + y**2 <= 1

is_satisfiable(disk & (x > 0) & (y > 0), [x, y])
# True

function_range(x + y, disk, [x, y])
# exact range

semialgebraic_minimize(x + y, disk, [x, y])
# exact global optimization result

semialgebraic_measure(disk, [x, y])
# pi
```

The central contract is **exactness rather than hidden numerical approximation**. A specialized backend may solve a problem before CAD; if it cannot, the package can conservatively fall back to a complete exact method where supported. See [Exactness and certification](concepts/exactness_and_certification.md).

## What semialg is not

`semialg` is an exact real-algebraic/semialgebraic toolkit, not a universal symbolic or numerical solver. In particular, it is **not** a general transcendental quantifier-elimination system, a floating-point nonlinear optimizer, an SOS/SDP search package, a complete computational-homology/roadmap package, or a universal symbolic integration engine. When an operation falls outside a certified fragment, the API is designed to return an explicit unknown/unsupported result or raise a documented exception rather than silently convert a heuristic numerical answer into an exact claim.

See [Limitations and scope](limitations.md) for the supported boundary of each subsystem.

## Start here

New users should read:

1. [Getting started](getting_started.md) — formulas, variables, regions, solving, optimization, and integration.
2. [Which function should I use?](guides/choosing_an_api.md) — choose between QE, projection, ranges, optimization, images, measure, and related APIs.
3. [Algorithm-selection map](concepts/algorithm_selection.md) — feasibility, QE, optimization, decomposition, positivity, integration, and conservative fallback.
4. [Certificate model](concepts/certificate_model.md) — candidate generation, exact reconstruction, verification, and replay.
5. [Understanding result objects](concepts/result_objects.md) — witnesses, attainment, certification, guarded branches, and `.select(...)`.
6. [Region representations](guides/region_representations.md) — choose between formulas, standard regions, `SemialgebraicRegion`, `CADRegion`, and structured cells.
7. [Exact versus numerical region results](guides/exact_vs_numerical_regions.md) — understand which layers are exact certificates and which are numerical geometry.
8. [Parameterized computation tutorial](guides/parameterized_computation.md) — ranges, optimization, measure, and algebraic-endpoint integration as parameters vary.

## Tutorials

- [Certified primary decomposition](tutorials/primary_decomposition.md)
- [Sparse SOS search and exact certification](tutorials/sos_certification.md)
- [Factorized projection and partial CAD](tutorials/cad_factorized_partial.md)

## Worked examples

- [Worked example gallery](examples/index.md) — 17 executable end-to-end examples covering QE, optimization, topology, moments, algebraic geometry, transforms, and parameterized integration.

## Guides: accomplish a task

- [Choosing an API](guides/choosing_an_api.md)
- [Parameterized computation](guides/parameterized_computation.md)
- [Performance](guides/performance.md)
- [Errors and failure modes](guides/errors_and_failure_modes.md)
- [Symbol handling](guides/symbol_handling.md)
- [Region invariants](guides/region_invariants.md)
- [Region representations](guides/region_representations.md)
- [Reuse one CAD](guides/cad_reuse.md)
- [Exact versus numerical regions](guides/exact_vs_numerical_regions.md)

Existing task-oriented guides remain available for [function ranges](function_range.md), [optimization](optimization.md), [region integration](region_integration.md), [geometry queries](geometry_queries.md), [region operations](region_operations.md), [solving](solving.md), and [decision procedures](decision_procedures.md).

## Learn semialgebraic geometry

If the mathematics is new to you, start with the [Introduction to semialgebraic geometry](semialgebraic_geometry/introduction.md), then continue to [Concepts and theory](semialgebraic_geometry/theory.md), [Algorithms and techniques](semialgebraic_geometry/algorithms.md), and [Applications](semialgebraic_geometry/applications.md). These pages explain the theory behind the APIs rather than assuming prior knowledge of real algebraic geometry.

## Concepts: understand the mathematics and guarantees

- [Exactness and certification](concepts/exactness_and_certification.md)
- [Certificate model](concepts/certificate_model.md)
- [How semialg chooses an algorithm](concepts/algorithm_selection.md)
- [Understanding result objects](concepts/result_objects.md)
- [Canonical representations and input conventions](concepts/representations.md)
- [CAD concepts](concepts/cad.md)
- [Parameter stratification](parameter_stratification.md)
- [Quantified expressions](quantified_expressions.md)

## Reference: look up an API family

- [API overview](api_overview.md)
- [Public API index](reference/public_api.md) — all exported names, checked against `semialg.__all__`.
- [Exception hierarchy](reference/exceptions.md)
- [Decision and QE](reference/decision_and_qe.md)
- [CAD](reference/cad.md)
- [Solving and sampling](reference/solving_and_sampling.md)
- [Optimization and ranges](reference/optimization_and_range.md)
- [Regions and geometry](reference/regions.md)
- [Integration and moments](reference/integration_and_moments.md)
- [Algebraic computation](reference/algebraic.md)
- [Parameters and conditional results](reference/parameters.md)
- [Applications](reference/applications.md)

## Architecture and quality

These pages are intended primarily for contributors:

- [Architecture design](architecture/design.md)
- [Code quality](architecture/code_quality.md)
- [Decision contracts](quality/decision_contracts.md)
- [Reference regression suite](quality/reference_regression_suite.md)
- [Robustness](quality/robustness.md)
- [Region testing strategy](quality/region_testing.md)
- [Region capability matrix](quality/region_capability_matrix.md)
- [Certificate matrix](quality/certificate_matrix.md)
- [Public API testing](quality/public_api_testing.md)
- [Coverage measurements](quality/coverage_measurements.md)
- [Changelog](changelog.md)
- [Region-computation architecture](architecture/region_computation.md)

## What “exact” does and does not mean

Exactness does **not** imply that every problem is cheap. Complete CAD/QE has severe worst-case complexity. `semialg` therefore uses structure-aware presolve and specialized exact backends where possible. In particular, **quadratic virtual substitution is a real QE/witness backend** and is tried for supported low-degree quantified problems before falling back to CAD. Zero-dimensional algebraic systems can use RUR machinery.

For current boundaries of the implementation, see [Limitations](limitations.md) and the [Feature matrix](feature_matrix.md).


## Algorithm provenance

See [Algorithm references](references.md) for the literature behind CAD, ARS-style real solving, regular chains, Zeng-family polynomial decisions, and SOS/Gram certificates.

- [Toric and lattice algebra](toric_algebra.md)
