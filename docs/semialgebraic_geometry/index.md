# Mathematical overview of semialgebraic geometry

This guide starts from polynomial equations, introduces real algebraic geometry and quantifier elimination, and develops the ideas behind exact roots, cylindrical algebraic decomposition (CAD), connectivity, optimization and positivity certificates. It is designed for mathematicians, scientists, and Python users who want to understand both **why the methods work** and **what `semialg` can actually compute**.

The chapters extend—not replace—the existing [introduction](introduction.md), [concepts and theory](theory.md), [algorithm survey](algorithms.md), and [applications](applications.md). The mathematical results are general; implementation claims are qualified and linked to the actual API, tests, or [limitations](../limitations.md).

## Recommended learning path

1. [Introduction](introduction.md) — what semialgebraic sets are and why inequalities matter.
2. [Algebraic foundations](foundations.md) — ideals, varieties, radicals, Gröbner bases and dimension.
3. [Real closed fields and quantifier elimination](real_geometry.md) — Tarski–Seidenberg, logical projection and tame geometry.
4. [Exact algebraic computation](algebraic_computation.md) — Sturm sequences, resultants, algebraic roots and signs.
5. [A complete worked CAD](cad_worked_example.md) — projection, lifting, 21 cells and a quantified projection.
6. [Topology, triangulation and roadmaps](topology.md) — finiteness, local geometry, Euler characteristic and certified scope.
7. [Optimization and positivity](optimization_positivity.md) — attained extrema, ellipsoid shortcuts, SOS and SDP certificates.
8. [Complexity and scope](complexity.md) — decidability, doubly exponential CAD costs and practical strategy.
9. [Exercises and selected answers](exercises.md) — test understanding with mathematical examples.

For a shorter route, read [concepts and theory](theory.md) followed by [algorithms and techniques](algorithms.md). For package-level entry points, begin with [getting started](../getting_started.md) and [choosing an API](../guides/choosing_an_api.md).

## How to read the examples

Some examples use SymPy for reproducible algebraic calculations; others are pencil-and-paper proofs. Neither should be mistaken for a claim that the same operation is directly exposed by a particular `semialg` API. Package-specific workflows are linked to existing guides, the [feature matrix](../feature_matrix.md), and the [native CAD proof contract](../native_cad_proofs.md).

## Mathematical references

The mathematical treatment follows standard sources on algebraic geometry, real closed fields, and computational real algebraic geometry. The [bibliography](../references.md#mathematical-background-and-further-reading) identifies the relevant books and papers. Package-specific claims are linked to implementation guides and tests, with the [limitations](../limitations.md) documenting where general mathematical theorems exceed the implemented algorithms.
