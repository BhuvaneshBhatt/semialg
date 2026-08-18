# Feature matrix

This matrix summarizes the current public API. “Initial” means the public API exists and works for supported exact cases, but the implementation is intentionally conservative outside that scope.

| Area | Status | Notes |
|---|---:|---|
| CAD construction | Initial | Public CAD and text renderers are available. |
| Complete QE | Initial | CAD-backed real quantifier elimination in supported cases. |
| Satisfiability / implication / equivalence | Initial | Exposed through `is_satisfiable`, `is_tautology`, `implies`, and `equivalent`. |
| Semialgebraic solving | Initial | `solve_semialgebraic` returns feasibility and representative samples. |
| Sampling and sign evaluation | Initial | Exact/rational/algebraic cases where supported. |
| Root classification with parameters | Initial | Strongest for univariate linear/quadratic families. |
| Parameter solvability conditions | Initial | Useful fast paths plus CAD/QE-backed workflows where supported. |
| Function range | Initial+ | Uses semialgebraic image projection for polynomial, rational, `Abs`, `sqrt`, `Min`, `Max`, and simple `Piecewise` expressions. |
| Symbolic optimization | Initial | Exact univariate and selected low-dimensional cases. Full global polynomial optimization remains future work. |
| Region operations | Initial | Boolean operations, boundary/closure/interior/dimension/components for supported cases. |
| Region predicates | Initial | Subset/equality/disjointness/boundedness/closedness/compactness wrappers. |
| Region integration | Initial+ | Standard shapes, selected vertical slices, reduction pieces, symbolic/numeric/auto modes. General CAD-cell integration remains future work. |
| Intrinsic measure | Initial | Finite point sets, graph curves, and circles in supported cases. General strata remain future work. |
| Moments / centroid / covariance | Initial | Built on the region integration engine. |
| Boolean and piecewise simplification | Initial | Semantic redundancy and branch pruning. |
| Assumption-based simplification | Initial | `Abs`, `sqrt(square)`, `Min`, `Max`, and `Piecewise` in supported cases. |
| Plotting / discretization | Planned | CAD-backed plotting and meshing are future work. |
| Roadmaps / topology | Planned | Connected components are preliminary; full roadmaps and Betti numbers are future work. |
| SOS / SDP certificates | Planned | Verification and certificate search are future work. |
