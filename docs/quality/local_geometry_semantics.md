# Local geometry quality status

The local-asymptotic geometry APIs use exact semialgebraic machinery for closure,
sign, image, witness, optimization, and local-topology questions.

## Current guarantees

- Complete-CAD sector samples are checked against their current exact algebraic
  bounds before a cached separator is reused.
- Real `find_instance` calls lower sign assumptions on public SymPy symbols to
  explicit constraints over assumption-neutral internal variables.
- `local_components` localizes a punctured representative before connectivity
  analysis so branches that reconnect away from the base point are not merged
  solely because of that distant connection.
- `local_preimage` obtains an exact fiber witness and certifies that no distinct
  source point exists before returning a single source germ.
- `local_sign_strata` prunes unrealizable sign prefixes incrementally.
- `parameter_strata` partitions feasibility by default and accepts ordered
  caller-defined semialgebraic invariant conditions.

## Validation focus

Regression tests cover the universal-QE formula that exposed stale CAD sector
sampling, sign-bearing Symbol assumptions, distant branch reconnection, unique
and nonunique preimage fibers, sign-stratum equivalence with exhaustive small
systems, and disjoint custom parameter strata. Slow topology contracts include
nodes, cusps, tacnodes, translations, and other singular local models.

## Remaining limits

Local connectivity can require expensive CAD and algebraic sign work on singular
multivariate sets; these cases remain in the slow test corpus. Curve selection
can return a theorem-backed existence certificate when a simple explicit
polynomial witness is not found. Rational pullback operations such as
`vanishing_order` require an explicit curve witness.

`parameter_strata` makes overlapping caller cases disjoint in insertion order.
Callers should therefore put more specific invariant cases before more general
ones.
