# Performance research and benchmark policy

Semialg accepts a specialized algorithm only after a baseline/candidate benchmark on the same corpus. Performance work records wall time, peak Python allocation, peak process RSS where available, route selection, and deterministic work counters. Fresh processes are used for cached algorithms.

## Rational-univariate construction

The permanent corpus is `semialg.benchmarks.rur_corpus`; `scripts/benchmark_rur.py` is the runner. The corpus includes coupled, grid, and three-variable zero-dimensional systems.

Columnwise variable multiplication matrices replace the all-pairs quotient multiplication tensor. For `n` variables and quotient dimension `D`, the construction performs `n D` variable-times-basis reductions instead of `D^2` basis-pair reductions. Downstream RUR semantics are unchanged.

A Krylov/incremental minimal-polynomial kernel is available for measurement and future linear-algebra work, but is not used by the RUR constructor. Straight SymPy rank/solve operations were slower than SymPy's characteristic-polynomial implementation on the benchmark corpus and can also expose a non-cyclic trial vector. A production incremental route needs a packed echelon implementation and separating-form verification before it can replace the current characteristic-polynomial path.

Matrix density is measured explicitly. Sparse storage is considered only below a configurable density threshold; the current RUR corpus is moderately dense, so dense matrices remain the default execution representation.

## CAD and QE

Variable-order benchmarking uses the Wilson CAD bank and projection-set size/sum-of-total-degree as deterministic work metrics in addition to wall time. `chordal` ordering uses a minimum-fill chordal completion seed and converts its elimination order to CAD lifting order. It is an explicit strategy, not the automatic default: on some dense examples it produces a larger projection set than Brown ordering.

Quadratic virtual substitution revisits variables skipped earlier in a same-kind existential block after another elimination changes the formula. This can expose a newly quadratic variable without invoking CAD.

Reduced McCallum/Lazard projection with equational constraints and proof-carrying fallback, plus family-aware TTICAD, are already first-class Semialg algorithms. Performance work should benchmark and improve those implementations instead of adding parallel versions.

Local-projection CAD and cylindrical algebraic coverings remain research candidates. They should be considered only after the RUR, virtual-substitution, ordering, equational-constraint, and TTICAD benchmark suites establish where the remaining time is spent.

## Recorded benchmark findings

Fresh-process measurements on the reference environment show that columnwise RUR construction reduces the `grid_3x3` median from 0.3090 s to 0.1706 s and peak Python allocation from 1,042,293 to 669,957 bytes. The coupled three-variable case moves from 0.07315 s to 0.06592 s and from 414,587 to 338,137 bytes. Raw samples are stored in `docs/benchmark-results/rur-columnwise.json`.

The straightforward SymPy Krylov prototype is retained as a measured kernel, not selected by the RUR constructor. It was slower and more allocation-heavy than `Matrix.charpoly` on every current corpus case; raw results are in `docs/benchmark-results/rur-minimal-polynomial.json`.


## EC/TTICAD and ordering follow-up

Reduced EC and TTICAD projection now avoid repeated algebraic expansion while matching already-expanded equational constraints. The benchmark in `benchmark-results/ec-selection.json` measures the shared EC-selection hot path.

Automatic Brown-versus-chordal ordering uses graph fill first and an abstract multidegree projection-growth estimate when structural fill ties. It does not construct resultants, discriminants, or trial projection towers. The release benchmark over the frozen 38-case Wilson subset is recorded in `benchmark-results/ordering-selector-degree-growth.json`.


## Packed RUR elimination experiment

`packed_krylov_minimal_polynomial` maintains normalized pivot rows and the matching Krylov-combination rows incrementally. Each new Krylov vector is reduced once, so the loop contains no repeated rank computation or generic linear solve. It is substantially faster and smaller than the earlier high-level Krylov prototype, but the current SymPy `charpoly` implementation remains faster on the cyclic companion benchmark. Production RUR construction therefore continues to use `charpoly`; the packed kernel remains an experimental exact kernel. See `benchmark-results/packed-rur-kernel.json`.

Local-projection CAD and cylindrical algebraic coverings are not part of this work.


## Large-domain packed RUR benchmark

The packed Krylov experiment now has a coefficient-domain kernel and a permanent
large-dimension runner (`scripts/benchmark_rur_large.py`).  The benchmark compares
the same cyclic matrix family over `QQ` and `GF(65521)` at dimensions 16 through
128 and records time, Python allocation peak, and process RSS.

No sustained runtime crossover was found.  Over `QQ`, domain `charpoly` is already
faster beyond the smallest case and its advantage grows with dimension.  Over
`GF(65521)`, packed elimination wins at dimensions 16 and 24, but `charpoly`
becomes faster at dimension 32 and remains faster through 128.  Packed elimination
does retain a memory advantage: at dimensions 32, 64, and 128 its measured Python
allocation peak is lower in both domains.  Production RUR therefore continues to
use `charpoly`; the packed kernel remains useful as a low-memory research kernel.
Raw measurements are in `benchmark-results/packed-rur-large-domain.json`.

## Structural Brown/chordal selector

The automatic Brown/chordal selector no longer performs a trial projection level
when graph fill ties.  It propagates only multidegree vectors using standard
degree bounds for resultants and discriminants.  This estimates downstream degree
growth without constructing projection polynomials.  On the frozen 38-case
Wilson subset, the structural selector chooses the better complete Brown/chordal
projection score in all 38 cases; Brown and chordal differ on 13 of those cases.
The earlier selector selected the better candidate in 33 of 38.  The benchmark is
recorded in `benchmark-results/ordering-selector-degree-growth.json`.
