# Feature matrix

This matrix separates three questions that are easy to conflate in symbolic software:

- **Exactness** — whether successful results are represented and certified without fixed-precision guesses.
- **Completeness** — whether the implemented algorithm decides every input in the stated mathematical fragment, assuming sufficient resources.
- **Practical scope** — the fragment for which the current implementation is intended to be useful.

“Partial” therefore does not mean “approximate.” A partial feature can still be exact: it may simply decline inputs outside its implemented fragment.

| Capability | Exactness | Completeness | Practical scope / notes |
|---|---|---|---|
| CAD construction | Exact | Complete for supported polynomial input when the complete backend finishes | Real polynomial sign-invariant CAD; reduced paths are used only with certified side conditions. |
| Complete polynomial QE | Exact | Complete over real-closed-field polynomial formulas when the complete CAD backend finishes | Subject to CAD's severe worst-case complexity. |
| Satisfiability / implication / equivalence | Exact | Complete through CAD fallback for supported semialgebraic formulas | Specialized exact backends may answer first. |
| Semialgebraic solving | Exact | Partial as a representation of the full solution locus | Exact feasibility and representative samples; richer cell/component output is available where constructed. |
| Sampling and sign evaluation | Exact in exact mode | Partial | Rational/algebraic samples are certified; `exact=False` is explicitly numerical. |
| Root classification with parameters | Exact | Partial-to-broad | Specialized exact linear/quadratic/cubic/quartic classifiers are retained; degree >=5 uses general subresultant/Sturm parameter stratification when the parameter CAD succeeds. |
| Parameter solvability conditions | Exact | Partial | Specialized fast paths plus complete CAD/QE where supported. |
| Parameter-stratified results | Exact | Partial by producer | Exact guards, branch certification, specialization, coverage/disjointness metadata. |
| Quadratic virtual substitution | Exact | Complete only for the implemented low-degree fragment | Used as a real QE/witness backend before CAD when applicable. |
| Zero-dimensional RUR solving | Exact | Partial | Finite polynomial systems over supported exact coefficient fields. |
| Function range | Exact | Partial | Polynomial/rational and selected `Abs`, `sqrt`, `Min`, `Max`, and `Piecewise` forms. |
| Function convexity | Exact | Partial | Polynomial and graph-supported semialgebraic functions on certified convex domains; Hessian/sign fast paths, affine-relative reduction, epi/hypograph reasoning, Jensen/QE fallback, and exact parameter strata. Unsupported transcendental graphs return `unknown`. |
| Polynomial optimization | Exact | Partial | Active-set/KKT, singular and positive-dimensional loci, RUR candidates, exact comparison, CAD global certification. |
| Region Boolean operations | Exact | Partial | Unified formula/`SemialgebraicRegion` operations plus reusable shared-CAD Boolean combination through `CADRegion`. |
| Region predicates | Exact | Partial | Subset/equality/disjointness and supported boundedness/closedness/compactness paths; convexity has exact affine, 1-D, quadratic, polynomial-Hessian, domain-relative, topology/witness, and complete-QE stages. |
| Topology | Exact | Partial | Connected components, exact CAD cell complexes with codimension-one incidence, local dimension, compact-support Euler characteristic, certified zero/one-dimensional and compact-convex BPR roadmaps, b0 generally, b1 for compact sets of dimension at most one, and trivial positive Betti numbers for certified compact convex sets; oriented higher homology and higher-dimensional critical-point roadmaps remain unimplemented. |
| Explicit path construction | Exact where returned | Partial | Explicit in one dimension; higher-dimensional output is a certified cell/connector chain rather than a full roadmap parameterization. |
| Region integration | Exact geometry; evaluation may be symbolic or explicitly numeric | Partial | Standard shapes and typed CAD-cell iterated integrals; exact symbolic antiderivatives may remain unevaluated. |
| Intrinsic measure | Exact | Partial | Certified regular CAD graph strata with induced Hausdorff metric; singular/uncertified strata are explicitly declined. |
| Moments / centroid / covariance | Exact when underlying integral is exact | Partial | Built on region integration. |
| Plotting / discretization | Numerical presentation over exact CAD provenance | Not applicable | Higher-dimensional CAD-cell simplicial meshing, shared-vertex conforming assembly, and delineable curve/surface evaluation; never an exact topology certificate. |
| Affine transforms | Exact | Complete for nonsingular square affine maps; otherwise QE-limited | Invertible maps use inverse substitution; singular/general images use existential QE. |
| Projection / image / preimage / fiber | Exact | Partial | First-class `SemialgebraicRegion` methods; projection/image use exact existential QE and preimage/fiber use substitution where applicable. |
| Bounding boxes and distances | Exact | Partial | Built on exact global optimization. |
| Support functions / width / diameter | Exact | Partial | Exact optimization with structural fast paths where available. |
| Singular locus / tangent space | Exact | Partial | Jacobian criterion and Zariski tangent-space machinery. |
| Nonsmooth locus / active boundary strata | Exact | Partial | Algebraic singularities, transverse inequality-boundary corners/ridges, and pairwise-disjoint active-inequality strata. |
| Bounded parametric covers | Exact | Structural | Bounded standard regions use natural charts; arbitrary formula regions can be represented exactly inside explicit finite clipping bounds without CAD. |
| Structural dimension from charts | Exact | Fallback-complete | Certified parameter-domain dimension and Jacobian rank are reused before exact CAD fallback. |
| Rich boundary strata | Exact | Complete CAD metadata | Boundary cells record dimension, inclusion/exclusion status, active inequality residuals, and retain the reusable boundary CAD. |
| Smooth-variety nearest point | Exact | Specialized | Pure smooth algebraic varieties use an exact distance-critical system before general optimization fallback. |
| Affine-map analysis | Exact/generic with parameters | Structural | Rank, determinant, inverse certification, and scaled-isometry recognition without QE. |
| Parametric map degree | Exact generic algebraic | Rational maps | Computes generic complex fiber degree; real-domain restrictions remain separate. |
| Low-dimensional affine box clipping | Exact | Specialized | One- and two-dimensional affine subspaces are clipped by exact facet intersections without CAD. |
| Tangent cone | Exact | Complete for the implemented ideal-theoretic construction when Gröbner elimination finishes | Saturated m-adic deformation; practical limit is elimination cost. |
| Critical values / positive-dimensional optimizer loci | Exact | Partial | Zero-dimensional KKT/singular candidates plus attained extrema; positive-dimensional projected KKT and singular active-boundary loci can fall through exact function-range/image-CAD analysis, with exact endpoint witness recovery when attained. `critical_values` itself still does not expose every symbolic image component. |
| Boolean and piecewise simplification | Exact | Partial | Relational atoms are treated propositionally by generic Boolean simplification; semantic pruning uses exact semialgebraic checks. |
| Assumption-based simplification | Exact | Partial | Selected `Abs`, square-root, `Min`, `Max`, and `Piecewise` rewrites. |
| Robust design / tolerance analysis | Exact | Partial | Semialgebraic operating domains. |
| Symbolic-math validation | Exact | Partial | Identity/formula/range validation with counterexamples where available. |
| Numerical-optimizer benchmark oracle | Exact reference value plus numerical comparison | Partial | Exact optimum serves as the oracle; tolerance checking of reported numeric values is explicit. |
| Polynomial control stability | Exact | Partial | Strict continuous-time Hurwitz stability regions for real polynomial characteristic equations. |
| Polynomial safety invariants | Exact | Partial | Verification of supplied invariants for discrete polynomial systems. |
| Polynomial probability | Exact when integration succeeds | Partial | Certified nonnegative polynomial densities over supported regions. |
| Lyapunov / barrier verification | Exact | Partial | Verification of supplied polynomial certificates. |
| Polynomial sensitivity | Exact | Partial | Derivative ranges and coordinate-wise sign/monotonicity classification. |
| Constraint redundancy / feasibility diagnostics | Exact | Partial | Implication-based redundancy, witnesses, and supported infeasible-core diagnostics. |
| Exact sign on regions | Exact | Broad semialgebraic fragment | `function_sign` plus positive/nonnegative/negative/nonpositive/zero/nonzero proof APIs; parameter-dependent proofs return exact guarded strata. |
| Exact sign partitions | Exact | Univariate semialgebraic functions | `function_sign_partition` returns connected positive/zero/negative regions with exact algebraic boundaries and parameter-first conditional partitions. |
| Function-property shared analysis / fast paths | Exact | Semialgebraic function queries | Reuses domain/derivative/Hessian/partition work; univariate rational sign partitions and curvature use dedicated exact paths before general CAD/QE. |
| Symmetric matrix definiteness | Exact | Polynomial/semialgebraic domain | Constant LDL/inertia, principal-minor PSD/NSD tests, Sylvester PD/ND tests, parameter conditions, and counterexample points/vectors where recoverable. |
| Matrix rank stratification | Exact | Determinantal | Constant rank on a region and parameter-space rank strata via exact minors. |
| Relative strict feasibility | Exact | Affine conjunctive systems | Explicit affine equalities define the hull; genuinely inactive affine inequalities are strictified and exact witnesses are returned when feasible. |
| Parametric affine reduction | Exact | Affine equalities with parameter-only pivots | Splits zero/nonzero pivot loci and returns reversible substitutions without generic-nonzero assumptions. |
| Algebraic root functions | Exact | Partial | Ordered-root identity, specialization, comparison, derivatives, regularity, guarded radical presentation. |
| CAD performance planning | Heuristic cost model; exact solver remains authoritative | Not applicable | Structural scores, projection-aware estimates, bounded pilot lifting, and process/solve-scoped caches. |
| Gröbner finite-variety CAD | Exact | Zero-dimensional common polynomial equality ideals | FGLM triangular projection, compatible-section lifting, exact residual-constraint pruning; full Collins fallback otherwise. |
| General transcendental QE | — | Not implemented | Outside the real-closed-field scope. Specialized transcendental solvers are separate and partial. |
| SOS certificate verification / optional search | Exact verification; optional heuristic search | Rational/algebraic Gram certificates | `semialg` verifies Gram identities and PSD exactly; optional `symbopt` search is accepted only after exact reconstruction and independent verification, with Zeng/ARS/CAD fallback otherwise. |

## Presolve and backend selection

Complete QE applies conservative affine substitution and exact Fourier–Motzkin elimination when safe. Automatic CAD ordering is quantifier-aware. Low-degree virtual substitution and zero-dimensional RUR solving can avoid CAD entirely. These are performance choices only: they do not weaken the exactness contract.

## Parameter-dependent computation

First-class stratified results are available for several parameter-dependent operations. Some results retain exact quantified relations rather than launching a second expensive QE solely to produce a compact `Piecewise` expression.

For operational details, see [Exactness and certification](concepts/exactness_and_certification.md), [How semialg chooses an algorithm](concepts/algorithm_selection.md), and [Limitations](limitations.md).


## Additional capabilities

| Capability | Exactness | Coverage | Notes |
|---|---|---|---|
| Arbitrary-degree parameter root counts | Exact when stratification succeeds | Partial-to-broad | Uses subresultant/Sturm coefficient sign invariance plus parameter-space CAD; compact low-degree classifiers remain preferred. |
| Semialgebraic graph conditions | Exact | Broad recursive fragment | Supported graph functions may occur inside relational constraints and `Piecewise` branch conditions, with auxiliary variables eliminated by exact CAD/QE. |
| Function monotonicity | Exact for supported univariate semialgebraic functions; derivative-sign fast path, strictness via derivative-zero dimension, pairwise graph/QE fallback, automatic parameter conditions |
| Exact smoothness loci | Exact for polynomial/rational functions and univariate supported semialgebraic joins | `function_smoothness` reports continuity, C^k order, smoothness, and exceptional loci. |
| Mapping properties | Exact for polynomial/rational and graph-supported semialgebraic maps | `function_mapping_properties` certifies injectivity/surjectivity/bijectivity with image and witnesses. |
