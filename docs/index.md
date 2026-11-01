# semialg documentation

Start with [getting started](getting_started.md), then choose a [task-oriented API](guides/choosing_an_api.md). [Worked applications](worked_applications.md) state independent expected answers. [Package ownership](package_ownership.md) explains how exact algebra, general geometry and real decisions fit together.

| Reader question | Entry point |
| --- | --- |
| What does a certified, unknown or unsupported result mean? | [Exactness and certification](concepts/exactness_and_certification.md) |
| Which representation or algorithm should I choose? | [Concepts](concepts/index.md), [task guide](guides/choosing_an_api.md) |
| Which outputs occur, and which values are only bounds? | [Function range](function_range.md), [exact images and bounds](range-images.md) |
| How do I handle exceptional parameters? | [Parameterized computation](guides/parameterized_computation.md) |
| How do I preserve a branch or an input hole? | [Graph verification](function_verification.md) |
| What are the signatures and return contracts? | [Reference](reference/index.md), [expert API](reference/expert_api.md) |
| What cannot be decided by this interface? | [Feature matrix](feature_matrix.md), [limitations](limitations.md) |
| How are examples and certificates checked? | [Testing policy](quality/testing_philosophy.md) |

The concepts section explains the mathematics, guides develop workflows, and the reference owns signatures and result attributes. Overlapping overview pages link to these entry points instead of maintaining separate feature promises. The full topic tree remains available in the navigation.

- [Certificate-field audits](certificate_audits.md)

- [Geometry for optimization and parameter families](optimization_geometry.md)


## Learn the mathematics behind semialg

The [mathematical overview of semialgebraic geometry](semialgebraic_geometry/index.md) provides a guided introduction to ideals and varieties, real closed fields, quantifier elimination, Sturm sequences, a complete worked CAD, finite topology, roadmaps, positivity, optimization and complexity. It complements the existing algorithm-focused documentation with proofs, examples, and exercises, and distinguishes general mathematical theorems from supported package capabilities.
