# Certificate and replay matrix

Exact computation and independent proof replay are separate capabilities. `replay_certificate` dispatches by the defining class; a positive result means only the checks stated below. Its `verified=None` means an unregistered type or insufficient replay payload, while `False` means a check failed. Exceptions and unknown solver results are not proofs of a false mathematical statement.

| Object / operation | Replay entry point | What is actually checked | Boundary |
| --- | --- | --- | --- |
| Native `cad_algorithms.decomposition.CompleteCAD` | `cad_algorithms.verify_cad` | Projection closure, local delineating families, root selectors, native valuations and stack coverage | Full ambient and supported recursive EC formats; finite-variety and TTI formats are excluded |
| Public `decomposition.cylindrical.CADResult` | `replay_certificate` | Input polynomial family, native ambient geometry, selected-cell truth, canonical output and function tree | Ordinary polynomial cell selection and components; preprocessed, finite-variety and topological wrappers return `None` |
| Structural `model.CADResult` | `replay_certificate` | Structural invariants via `validate_cad_result` | Does not re-prove geometric delineability |
| `model.QEResult` | `replay_certificate` | Underlying CAD invariants | Does not independently prove quantified truth or reconstruction |
| `partial.qe.LazyQEResult` | `replay_certificate` | Partial trace reconstruction, full geometric lifting and independent quantified truth propagation | Replay may cost a full CAD; timing and work counters are diagnostic |
| `partial.qe.LazyResolveResult` | No standalone proof replay | Exact decision and work statistics are returned | No independently replayable lazy proof object is retained |
| `optimization_results.OptimizationResult` | `replay_certificate` | Reruns the optimizer, compares value and attainment, checks retained optimizer feasibility | Computational recomputation, not a standalone proof trace; missing constraints give `None` |
| `optimization_results.FunctionRangeResult` | `replay_certificate` | Retained attained endpoint membership | Partial consistency check; does not prove image completeness or exclude missing interior values |
| Polynomial root interval | `algebraic.verify_root_interval_certificate` | Exact distinct-root count and endpoint semantics | Supply expected polynomial, variable, interval and endpoint flags to bind the proof to a caller's statement |
| Branch correspondence | `topology.branch_correspondence.replay_branch_correspondence` | Source spectrum, target multiplicities, barriers and continuation ranks | Checks its stored statement; callers must also bind that statement to the intended problem |
| Recursive EC projection | `cad_algorithms.projection.ec_hypotheses.replay_ec_projection` | Explicit and resultant implications, degree/order families, operator provenance | Formula-bound projection replay; complete geometric checking uses `verify_cad` |
| SOS identity | `verify_sos_certificate` | Exact Gram identity and positive semidefiniteness | Search heuristics are separate from exact verification |
| Modular Gröbner, resultant and subresultant operations | Family-specific `algebraic.verify_*` functions | Exact reconstruction identities and retained modular evidence | Restricted to each verifier's coefficient domains |
| GTZ primary decomposition | `algebraic.gtz_primary.verify_gtz_primary_decomposition_certificate` | Localization, contraction, saturation and recursive finite components | Nested components are replayed through the enclosing proof |
| H-representation redundancy | `polyhedral.verify_h_redundancy_certificate` | Redundancy against the supplied representation | A different representation must be checked separately |
| CAD cell complexes, Euler characteristic and measure | Recompute from source CAD and exact pieces | Incidence, cell dimensions, alternating sums or integrals | No general unified standalone replay payload |
| Numerical mesh | Structural/provenance checks | Cell provenance and conformity | No exact geometry or isotopy certificate |

Use `result_diagnostics` for method choice and work counters. Diagnostics summarize construction; they cannot replace replay. The [mutation audit](../certificate_audits.md) reports active proof fields, inactive branches and mathematical versus constructor/schema rejections separately. Valid proofs are also checked after cache resets and serialization. The [root debugging tutorial](../tutorials/root_certificate_debugging.md) demonstrates statement binding and endpoint counts.
