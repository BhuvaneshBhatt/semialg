# Release notes

## 1.0.0 - 2026-09-12

`semialg` provides exact symbolic computation for real semialgebraic problems for Python 3.11–3.14.

### Capabilities

- Cylindrical algebraic decomposition, real quantifier elimination, satisfiability, implication, and equivalence.
- Exact semialgebraic solving with affine, univariate, virtual-substitution, Gröbner/RUR, and CAD/QE methods.
- Exact algebraic roots, finite polynomial systems, ideal decomposition, toric/lattice algebra, and "replayable" certificates.
- Region operations, topology, projection and images, distances, singular loci, tangent constructions, and standard geometric regions.
- Polynomial optimization, exact function ranges, integration, ambient/intrinsic measure, moments, centroids, and covariance.
- Parameter-stratified results and exact conditional answers for supported parameter-dependent problems.
- Process-local caches and algorithm selection that prefer cheap exact algebraic certificates before CAD/QE.

### Exactness contract

Successful certified paths use exact symbolic or algebraic reasoning. Numerical proposals may accelerate a computation but are not accepted as exact results without independent certification. Unsupported or undecidable cases return a documented unknown result or raise a documented exception rather than silently converting a numerical approximation into a proof.

### Known limitations

- CAD and complete QE have severe worst-case complexity.
- General transcendental QE is outside the package scope.
- Several topology, integration, optimization, and parameterized algorithms cover exact but restricted fragments.
- Affine-region transformation entry points do not uniformly reject complex-valued matrices at the input boundary.
- Some specialized recognizers depend on syntactic normalization and may decline algebraically equivalent presentations; generic certified methods remain available.

See the [capability matrix](https://github.com/BhuvaneshBhatt/semialg/blob/main/docs/feature_matrix.md) and [limitations](https://github.com/BhuvaneshBhatt/semialg/blob/main/docs/limitations.md) for detailed scope.
