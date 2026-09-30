# Local geometry for asymptotic analysis

Semialg exposes exact geometric information that limit and asymptotic algorithms
can consume without moving expansion logic into this package.

`local_geometry` combines local semialgebraic dimension with the Bouligand
tangent cone. `local_sign_strata` separates locally realizable sign patterns,
while `local_germ` and `local_components` represent the behavior of a set near a
closure point. `curve_selection` returns a certified polynomial approach curve
when one is found cheaply and otherwise returns a theorem-backed certificate
that a semialgebraic approach curve exists.

`local_range` certifies extrema on an explicit local ball. `local_bound` can
also solve for a positive radius and independently verify the resulting concrete
bound before returning a certificate.
`parameter_strata` gives exact parameter conditions for feasibility. Curve
witnesses support rational-pullback `vanishing_order` and `contact_order`.
`blowup_charts` constructs ordinary or weighted polynomial charts.
`local_image` and `local_preimage` preserve local map geometry through exact
semialgebraic image/preimage operations. `path_independent` certifies equality
of an expression throughout a germ.

## Scope and limitations

Curve selection uses a bounded explicit polynomial search first. When that
search misses, the result certifies existence from the exact closure hypothesis
and the semialgebraic curve-selection theorem; the fallback sets `explicit=False`
so downstream code cannot use it as if a parametrization had been constructed.

`local_range` uses a caller-supplied radius. `local_bound` accepts a fixed radius
or searches for a positive one. Parametric QE proposes a witness and a separate
existential counterexample check validates the concrete neighborhood before it
is certified.

`parameter_strata` stratifies feasibility by default and accepts
`conditions_by_value` for caller-defined semialgebraic invariants. Overlapping
cases are made disjoint in insertion order, so specific cases should precede
more general ones.

`vanishing_order` and `contact_order` require a rational pullback in the curve
parameter. Weighted blow-up charts use positive integer weights.

`local_preimage` requires the target germ point to have one exact source point.
This avoids pretending that a multi-branch preimage is a single germ.
