

## Version 1.3.0

- Added certified local + approach geometry. Analyze whether semialgebraic sets approach a point or boundary, construct local germs and branches, study local sign behavior, and obtain proof oriented diagnostics.
- Added local geometry for asymptotic analysis. Supportfor certified approach regions, local bounds, weighted/directional charts, and related geometry needed for multivariate asymptotic analysis.
- Expanded parametrized semialgebraic analysis. Improved parameter conditions and stratification make it possible to partition parameter space into regions with stable geometric or logical behavior and identify exceptional/transition loci.
- Expanded CAD and TTICAD capabilities. Improve equational-constraint handling, projection and lifting, nullification handling, and truth-table-invariant CAD support broaden the class of formulas that can be handled effectively.
- Improved exact real-algebraic solving. Root isolation, sign determination, rational univariate representation (RUR), algebraic-number ordering, and exact solution reconstruction are more capable and robust, particularly for larger and more difficult systems.
- Added cheap exact Boolean-structure certificates and vacuous-quantifier pruning before algebraic decision procedures.
- Added local topology tools: local connected components, branches, approachability, germs, sign strata, etc
- CAD, TTICAD, parameter analysis, RUR ordering, local geometry, satisfiability, and several other exact computations have been optimized, including fixes for bad slowdowns on some large Boolean and parameterized problems.
- Improved robustness on degenerate and exceptional cases. Adds substantially broader handling and regression coverage for zero polynomials, nullification, exceptional parameter values, singular/local behavior, exact boundary cases, and difficult algebraic inputs.
- Expanded documentation and examples, with substantially more guidance on CAD, parameterized problems, certification, approach/local geometry, and using local semialgebraic geometry in asymptotic analysis.

## Version 1.2.0

- Reuse native selected CAD cells directly for `SemialgebraicRegion.components()` instead of reconstructing a structured cylindrical solution, memoize the resulting component graph/formulas.
- Reduce CAD closure-incidence work by computing each closed source path and target sample assignment once, then reuse closure membership across cell-pair intersections.
- Centralize dependency-aware cylindrical sample closure in `CylindricalSampleContext`, with cached prefix/full assignments reused by nested-root incidence.
- Cache closed CAD paths and exact sample assignments on `NativeRegionAnalysis`; the solid-torus invariant regression remains exact and now completes in about 17.5 seconds in the validation environment.
- Strengthen exact rational-power graph contracts, including negative rational exponents and explicit odd real-root behavior.
- Route standard-region topology through certified finite triangulations before general CAD/roadmap machinery.

## Version 1.1.0

`semialg` provides exact symbolic computation for real semialgebraic and algebraic problems on Python 3.11–3.14.

### Decision procedures and solving

- Cylindrical algebraic decomposition, certified real quantifier elimination, satisfiability, implication, and equivalence.
- Exact semialgebraic solving with affine, univariate, virtual-substitution, Gröbner/RUR, and CAD/QE methods.
- Independent quantified variable blocks can be eliminated separately when their formulas are structurally independent.
- Equality-ideal Gröbner consequences are used only as equivalence-preserving preprocessing; Zariski projection is not substituted for real existential projection.
- Reusable parameterized CAD supports exact specialization and formula recovery.

### Algebra and certification

- Exact algebraic roots, finite polynomial systems, ideal decomposition, toric and lattice algebra, and "replayable" certificates.
- Polynomial-map implicitization, Zariski closure, certified radicalization, irreducible components, and reduced-component singular loci.
- Constraint implication, redundancy analysis, nonnegative-combination certificates, and component-aware constraint descriptions.
- Shared certified zero, equality, and sign decisions preserve uncertainty instead of converting unresolved symbolic questions to Python truth values.

### Geometry and topology

- Canonical geometric types such as `Point`, `Box`, `Simplex`, `Polygon`, `Polytope`, `Ball`, `Sphere`, `Zonotope`, and `TetrahedralComplex`.
- Polygonal and polyhedral boundary models support disconnected components, holes, shells, and cavities.
- Exact canonicalization recognizes common polygonal and polyhedral structures without convexifying nonconvex inputs.
- Exact convex hull, triangulation, mixed-cell tetrahedralization, random geometry generation, triangular subdivision, and geodesic refinement.
- Exact connected-component reconstruction handles open, closed, touching, equality-defined, and mixed-dimensional semialgebraic sets.
- Region operations include projection and images, distances, singular loci, tangent constructions, relative topology, and local dimension.

### Optimization, integration, and parameters

- Polynomial optimization, exact function ranges, integration, ambient and intrinsic measure, moments, centroids, and covariance.
- Parameter-stratified results and exact conditional answers for supported parameter-dependent problems.
- Parameterized ellipsoid validation retains unresolved symmetry and positive-definiteness conditions rather than rejecting undecidable symbolic cases.
- Process-local caches and algorithm selection prefer inexpensive exact certificates before complete CAD/QE.
