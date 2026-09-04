# Certificate and replay matrix

`semialg` separates an exact result from the evidence retained for independently replaying that result. `replay_certificate(result)` returns `verified=True`, `False`, or `None` when the result family does not retain enough replay payload.

| Result / operation | Exact result | Retained evidence | Replay support |
|---|---:|---|---|
| `CADResult` | Yes | CAD cells, projection/lifting data, diagnostics | Yes: CAD invariant validation |
| `QEResult` | Yes | Underlying CAD/certification metadata | Yes where CAD evidence is retained |
| `OptimizationResult` | Yes when certified | Method, optimum, points, constraint formula | Yes for retained modern results |
| Convexity certificates | Yes | Certificate-specific proof data | Yes |
| Function range | Yes when certified | Range metadata and attained endpoints | Partial replay |
| `CADPointLocation` | Yes | Cell index + exact tower signs | Directly auditable by substitution |
| `CADCellComplex` incidence | Yes | Source CAD closure relation | Reconstructible from source CAD |
| Euler characteristic | Yes | Cell dimensions / `f_vector` | Recompute from complex |
| Region measure | Yes when exact evaluation succeeds | CAD integration pieces where requested | Recompute from exact pieces |
| Numerical mesh | No | Source-cell provenance, conformity flag | Structural checks only; not an exact geometric certificate |

## Interpreting `verified=None`

`None` means that the public object does not contain enough information for the unified replay interface to reproduce the proof. It is not the same as a failed certificate.

## Diagnostics

Use `result_diagnostics(result)` for a stable summary of method choice, counters, cache statistics, and certification metadata. This is useful in regression tests because deterministic work counters are more portable than tight wall-clock thresholds.
