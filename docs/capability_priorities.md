# Capability priorities

The [parameter-cell QE tutorial](tutorials/parameter_cell_qe.md) explains shared linear dispatch, partial lifting and resource controls. The [coefficient-field tutorial](tutorials/algebraic_coefficient_fields.md) covers exact fields and embeddings. [Branch continuation](tutorials/failure_aware_closure.md) distinguishes finite limits, infinity and unavailable cases.

| Capability | Supported scope | Boundary |
| --- | --- | --- |
| Shared linear dispatch | Bounded affine Boolean formulas and alternating prefixes before virtual substitution; exact equality presolve and Fourier–Motzkin | Parameter-dependent leading coefficients and excessive Boolean expansion use general exact backends |
| Critical-locus number fields | Exact common fields for validation, active sets, dimension and multiplier elimination; quadratic/compositum and embedded-root coefficients | Symbolic coefficient fields need explicit exceptional-stratum semantics; large composita can be expensive |
| Public budgets | One cooperative budget through preprocessing, projection, lifting, failed attempts, output construction and cache returns | One arithmetic operation can overrun a deadline; hard limits require process isolation |
| Unified replay | Native ambient CAD; public polynomial cell-selection wrappers; partial free-cell QE with a full geometric/truth reference | Topological, preprocessed and finite-variety wrappers need additional proof payloads; range replay checks endpoint consistency only |
| Degree-loss continuation | One base parameter, exact Sturm germs and oriented reciprocal limits, including persistent multiplicities | Multibase degree loss and path-dependent limits require broader continuation proofs |
| Free-parameter partial CAD | Complete free-cell truth coverage with alternating quantified suffixes and full Lazard restart | Projection is still complete and can dominate cost |
| Roadmaps | Dimension at most one and supported compact convex sets | General higher-dimensional roadmaps require recursive critical-point constructions |
| Homology | Certified supported polyhedral and one-dimensional triangulations | General curved triangulation and boundary operators are absent |
| Intrinsic integration | Supported regular graph/chart cases | General singular strata and overlapping nonlinear chart covers require additional machinery |

Variable order, equation-aware projection, exact image reuse and specialized fragments improve practical cost but cannot remove CAD's worst-case complexity. General transcendental analysis belongs in funcprops, optimizer strategies in symbopt, positive-dimensional algebraic geometry in alggeom, and finite quotient arithmetic in algroots.
