## Version 1.5

Focused on `algroots` integration: native Gröbner defaults, coefficient-field towers, guarded parameter preprocessing, and backend-aware RUR and primary decomposition.

### Native Gröbner integration

- Use `algroots`’s verified native Gröbner backend by default for Gröbner-backed algebraic operations, including RUR construction, primary decomposition and parameter equality preprocessing.
- Add `groebner_backend` and `groebner_limits` options to the relevant entry points, preserving the selected policy through recursive calculations.
- Retain explicit SymPy construction through `groebner_backend="sympy"`. Native requests report unavailable kernels, exhausted budgets and verification failures without automatically switching backends.

### Coefficient fields and algebraic towers

- Add `as_domain()` and `to_coefficient_domain()` bridges from `semialg`’s rational-function fields and monogenic towers to `algroots`’s exact coefficient domains.
- Support finite primary decomposition over declared number fields, rational parameter fields and parameter-dependent algebraic towers. Returned certificates verify the source ideal, primary components, associated primes and component degrees.
- Extend RUR construction with an explicit `domain` argument, including supported parameter and tower fields. Generic RURs describe the declared coefficient field, while real-root extraction requires suitable specialization and real embeddings.
- Cache converted coefficient domains and detach mutable constructor sequences so later input changes cannot invalidate cached field representations.

### Parameter-dependent solvability

- Add verified Gröbner preprocessing of outer conjunctive equalities in `solvability_conditions`, including algebraic coefficients and parameter-dependent systems.
- Use specialization guards to distinguish generic calculations from exceptional parameter fibers. Retain the original equations and inequalities on exceptional fibers before complete real elimination.
- Preserve visible denominator restrictions when clearing rational equalities, preventing cancelation from introducing solutions at excluded poles.
- Keep Gröbner preprocessing budgets separate from subsequent CAD work. Algebraic elimination does not replace real quantifier elimination.

### Correctness and robustness

- Fix algebraic fiber-tower reduction when coefficient denominators are hidden inside expanded sums. Normalize rational remainders before passing them to polynomial-only reduction stages.
- Fix collisions between bound symbols inside algebraic-root representations against polynomial coordinates, preserving the selected coefficient embedding.
- Include coefficient domains, backend selections and operation limits in RUR and GTZ cache keys.

### Documentation, testing and packaging

- Update documentation and examples for backend selection, coefficient-field towers, exceptional fibers, certificate verification and mathematical limitations
- Expand regression tests for native/SymPy agreement, repeated roots, parameter poles, field embeddings, cache separation, budget failures and corrupted certificates.
- Add broader CI and distribution checks.

## Version 1.4

### CAD, quantifier elimination, and certification

- Add full Lazard projection and native valuation-aware lifting for multivariate polynomials with rational or exact real-algebraic coefficients. The implementation handles degree drops and positive-dimensional nullification, and complete Collins can restart with Lazard when ordinary lifting cannot be certified.
- Strengthen reduced CAD with recursively justified equational constraint (EC) implications and degree/order hypotheses. When these conditions cannot be certified, the solver falls back to an ambient decomposition.
- Improve quantifier elimination with bounded exact elimination of linear Boolean formulas, quantifier-guided partial/lazy lifting and reconstruction of simpler formulas from certified CAD sign data when possible.

### Regions, parameter families, and topology

- Add `bounding_cover` for inexpensive outer coordinate bounds and `optimization_cover` for structural parametrizations that retain boundaries, singular points and lower-dimensional pieces.
- Add `region_family` to partition parameter space by fiber behavior, including empty and dimension-changing fibers. Scalar fibers expose component counts and graph/band structure, supported ball and sphere families also have structural component and Betti-number data. General higher-dimensional fibers expose CAD cells and dimensions, not complete topological invariants.
- Add shortcuts for closure, interior and boundary of supported quadratic and Cartesian-product semialgebraic sets, avoiding CAD incidence computations when exact formulas suffice.

### Exact optimization

- Add certified extrema for univariate polynomial objectives on semialgebraic domains using complete real-root isolation, stationary points, and endpoint or infinite limits.
- Add support-function optimization of affine objectives over supported positive-definite ellipsoids, distinguishing attained optima from unattained bounds on open domains and bypass general CAD when their hypotheses apply.

### Algebra, performance and robustness

- Reuse prepared polynomial data, existing Gröbner bases and complete sign-table cache entries. Also reduce expression round trips in FLINT conversions and resultant computations.
- Improve function-field computations with fraction-free norm elimination, and strengthen exact algebraic-zero validation for rational univariate representation (RUR) points.

### Function graphs

- Add `verify_real_function`, `compare_real_graphs`, and `verify_real_composition` for checks of graph totality, single-valuedness, equivalence and composition on declared real domains. These use decision and quantifier-elimination paths, not numerical sampling.
- Correct Boolean graph negation and inactive `Heaviside` branch domains. Improve real coverage of backend curve parametrizations, retaining isolated real points and necessary domain restrictions. The new `alggeom` package handles algebraic decomposition search and reconstruction.

### Misc

- Update dependencies to `sympy>=1.14,<2`, `python-flint>=0.9,<1`, and `algroots>=0.3.0` (with `exprtest>=0.1.0` unchanged).
- Expand documentation for native CAD proofs, optimization geometry, parameter families, and mathematical limitations.

## Version 1.3.2

- Extended coarea support with vector-map Gram certification, parameter-event candidates and exceptional-fiber handling, added budget/dispatch regressions and measured certificate reuse.
- Added cached chart, cover, null-set and linear metric certificates, including some positive algebraic and canonical polar charts.

## Version 1.3.1

- Fixed affine chart detection under equivalent linear constraint bases.
- Fixed closed-interval wrapper around half-open Sturm counts, including left-root and zero-width endpoint cases.

## Version 1.3.0

- Added certified local + approach geometry. Analyze whether semialgebraic sets approach a point or boundary, construct local germs and branches, study local sign behavior, and obtain proof oriented diagnostics.
- Added local geometry for asymptotic analysis. Supportfor certified approach regions, local bounds, weighted/directional charts, and related geometry needed for multivariate asymptotic analysis.
- Expanded parametrized semialgebraic analysis. Improved parameter conditions and stratification make it possible to partition parameter space into regions with stable geometric or logical behavior and identify exceptional/transition loci.
- Expanded CAD and TTICAD capabilities. Improved equational-constraint handling, projection and lifting, nullification handling, and truth-table-invariant CAD support broaden the class of formulas that can be handled effectively.
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
