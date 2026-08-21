# Feature matrix

This matrix summarizes the current public API. “Supported” denotes functionality with a public API for the stated scope; “Partial” denotes intentionally conservative functionality that declines unsupported cases rather than approximating silently.

| Area | Status | Notes |
|---|---:|---|
| CAD construction | Supported | Public CAD and text renderers are available. |
| Complete QE | Supported | CAD-backed real quantifier elimination in supported cases. |
| Satisfiability / implication / equivalence | Supported | Exposed through `is_satisfiable`, `is_tautology`, `implies`, and `equivalent`. |
| Semialgebraic solving | Supported | `solve_semialgebraic` returns feasibility and representative samples. |
| Sampling and sign evaluation | Supported | Exact/rational/algebraic cases where supported. |
| Root classification with parameters | Partial | Strongest for univariate linear/quadratic families. |
| Parameter solvability conditions | Partial | Useful fast paths plus CAD/QE-backed workflows where supported. |
| Parameter-stratified results | Supported | First-class guarded values with branch certification, disjointness/coverage verification, specialization, selection, and Piecewise conversion. |
| Robust design / tolerance analysis | Supported | Exact existential feasibility, universal robustness, and violation regions over semialgebraic operating domains. |
| Symbolic-math validation | Supported | Exact identity, formula-equivalence, and proposed-range validation with counterexamples where available. |
| Numerical-optimizer benchmark oracle | Supported | Certified exact reference optima and tolerance checks for reported numerical objective values. |
| Polynomial control stability | Supported | Strict continuous-time Hurwitz stability regions for real polynomial characteristic equations. |
| Polynomial safety invariants | Supported | Exact initiation, inductiveness, and unsafe-state exclusion for supplied invariants of discrete polynomial systems. |
| Polynomial response surfaces | Supported | Exact extrema, ranges, gradients, stationary conditions, and threshold sets for polynomial surrogate models. |
| Polynomial model comparison | Supported | Exact signed discrepancy extrema, maximum absolute error, dominance, equivalence, and counterexamples for polynomial models. |
| Parameter regime analysis | Supported | Exact parameter strata for semialgebraic solvability and univariate real-root counts. |
| Polynomial probability | Supported | Exact normalized probabilities for certified nonnegative polynomial densities over supported semialgebraic regions. |
| Lyapunov verification | Supported | Exact verification of supplied polynomial Lyapunov functions for polynomial vector fields. |
| Barrier certificates | Supported | Exact verification of supplied polynomial continuous-time barrier certificates. |
| Polynomial sensitivity | Supported | Exact derivative ranges and certified coordinate-wise sign/monotonicity classification. |
| Constraint redundancy | Supported | Exact implication-based redundancy classification with witnesses. |
| Feasible-set diagnostics | Supported | Exact witnesses or irreducible infeasible constraint subsets. |
| Algebraic root functions | Supported | Certified ordered-root identity, exact specialization/evaluation, comparisons, derivatives, regularity, and guarded radical presentation. |
| CAD performance planning | Supported | Per-solve computation contexts plus bounded process caches, projection-aware order scoring with lifting/root, algebraic-degree and coefficient-height estimates, bounded pilot lifting, and conservative multi-level EC pruning. |
| Function range | Supported | Uses semialgebraic image projection for polynomial, rational, `Abs`, `sqrt`, `Min`, `Max`, and simple `Piecewise` expressions. |
| Symbolic optimization | Supported | Exact polynomial active-set/KKT pipeline with strict/redundant/inconsistent active-set pruning, singular/positive-dimensional locus handling, RUR candidates, exact comparison, CAD global certification, and first-class parameter-stratified optimum relations. |
| Region operations | Partial | Boolean operations, boundary/closure/interior/dimension/components for supported cases. |
| Region predicates | Supported | Subset/equality/disjointness/boundedness/closedness/compactness wrappers. |
| Region integration | Partial | Standard shapes plus typed arbitrary-dimensional CAD-cell iterated integrals; symbolic/numeric/auto evaluation remains dependent on integrability. |
| Intrinsic measure | Supported | Certified regular CAD graph strata use induced Hausdorff metric; singular/uncertified strata are explicitly stratified and declined for verified integration. |
| Moments / centroid / covariance | Supported | Built on the region integration engine. |
| Boolean and piecewise simplification | Supported | Semantic redundancy and branch pruning. |
| Assumption-based simplification | Supported | `Abs`, `sqrt(square)`, `Min`, `Max`, and `Piecewise` in supported cases. |
| Plotting / discretization | Partial | Public 1D/2D plotting/discretization helpers exist; higher-dimensional views remain unsupported. |
| Roadmaps / topology | Planned | Connected components are preliminary; full roadmaps and Betti numbers are future work. |
| SOS / SDP certificates | Planned | Verification and certificate search are future work. |

