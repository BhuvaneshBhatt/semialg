# Approach and local geometry

Multivariate local analysis needs more than independent coordinate ranges. A semialgebraic set can constrain coordinates jointly, can approach a point without containing it, and can split into several algebraic branches at a singularity. Semialg exposes exact operations for those questions so downstream symbolic algorithms can reason about *admissible approaches* instead of rectangular over-approximations.

## Family contract

**Mathematical return.** These APIs decide exact local/closure properties or return exact correlated images, local strata, and normalized proof diagnostics.

**Exactness and certification.** Closure, image, and local-stratum conclusions use exact semialgebraic/algebraic machinery. Diagnostic normalization reports the proof route of an existing result and does not itself upgrade heuristic metadata into a certificate.

**Algorithm.** Closure uses semantic topology; correlated images use exact graph projection/QE; angular images add an exact sphere constraint; local strata combine local branch geometry with certified global stratification.

**Complexity and limitations.** Graph projection, CAD/QE, and algebraic decomposition can be intrinsically expensive. Resource guards or unsupported exact fragments are reported explicitly instead of replaced by numerical guesses.

For a map \(F:S\to\mathbb R^m\), coordinate-wise ranges describe only projections of \(F(S)\). They do not recover the joint image. For example, the map \(t\mapsto(t,t^2)\) on \([-1,1]\) has coordinate ranges \([-1,1]\) and \([0,1]\), but the attainable pairs satisfy \(v=u^2\). `correlated_map_image` retains that relation exactly by projecting the graph of the map.

## point_in_closure

`point_in_closure(region, point, variables=None, *, strategy=None, return_result=False)` decides whether an exact point belongs to the Euclidean closure of a semialgebraic region. This is the appropriate predicate for asking whether admissible points occur arbitrarily near a candidate limit point. The point need not itself belong to the region.

```python
import sympy as sp
from semialg import point_in_closure

x, y = sp.symbols("x y", real=True)
punctured = sp.And(x**2 + y**2 < 1, sp.Ne(x**2 + y**2, 0))
assert point_in_closure(punctured, (0, 0), (x, y))
```

The default path uses CAD-semantic closure. `strategy="syntactic"` explicitly requests the atom-wise fallback. `return_result=True` returns `PointInClosureResult`, including the normalized region, point, closure formula, variables, and method.

## CorrelatedMapImageResult

`CorrelatedMapImageResult` is the immutable structured result returned by correlated-image operations. Its `formula` is the exact joint image in `image_variables`; `mapping`, `variables`, and `domain` record the source problem. `diagnostics` can include source/output dimensions, generic rank, image dimension, and critical-value geometry. `proof_trace` records the graph/projection route used to certify the image.

## correlated_map_image

`correlated_map_image(mapping, domain=True, variables=None, *, image_variables=None)` computes the exact *joint* image of a polynomial or rational map. It preserves algebraic correlation between output coordinates instead of replacing the image by a Cartesian product of scalar ranges.

```python
from semialg import correlated_map_image

t, u, v = sp.symbols("t u v", real=True)
image = correlated_map_image((t, t**2), sp.And(t >= -1, t <= 1), (t,), image_variables=(u, v))
assert sp.simplify(image.formula.subs({u: sp.Rational(1, 2), v: sp.Rational(1, 4)})) is sp.true
assert sp.simplify(image.formula.subs({u: sp.Rational(1, 2), v: sp.Rational(1, 3)})) is sp.false
```

The computation is exact graph projection and therefore can have the cost of quantifier elimination. Unsupported rational/algebraic fragments are declined by the underlying exact image machinery instead of approximated numerically.

## angular_map_image

`angular_map_image(mapping, variables, domain=True, *, image_variables=None, sphere=True)` specializes correlated image computation to directional variables. With `sphere=True`, Semialg adds the exact unit-sphere equation before projecting the map. This makes it suitable for direction-dependent leading terms, tangent-direction invariants, and other local analyses where independently varying direction coordinates would be unsound.

```python
from semialg import angular_map_image

u, v = sp.symbols("u v", real=True)
angular = angular_map_image((x, y), (x, y), image_variables=(u, v))
assert sp.simplify(angular.formula.subs({u: 1, v: 0})) is sp.true
assert sp.simplify(angular.formula.subs({u: 1, v: 1})) is sp.false
```

Set `sphere=False` only when the supplied domain already encodes the intended normalization or when the full unconstrained image is desired.

## LocalAlgebraicStrata

`LocalAlgebraicStrata` combines certified local branch geometry with global singular and local-dimension strata incident at the requested point. It prevents callers from assuming that a branch list alone captures component intersections or dimension jumps.

## local_algebraic_strata

`local_algebraic_strata(equations, point, variables=None, *, max_pieces=None)` returns exact local branches together with every certified singular/dimension stratum incident at `point`. For the crossing \(xy=0\) at the origin, the local branch geometry contains the two coordinate-axis branches and records their intersection.

```python
from semialg import local_algebraic_strata

strata = local_algebraic_strata((x * y,), {x: 0, y: 0}, (x, y))
assert strata.complete
assert strata.branch_geometry.singular
assert len(strata.branch_geometry.branches) == 2
```

`max_pieces` is a resource guard for decomposition. Local algebraic geometry can become expensive as degree, dimension, or the number of components grows.

## ProofDiagnostics

`ProofDiagnostics` is the stable normalized diagnostic view produced by `structured_proof_diagnostics`. It records a method name, de-duplicated proof steps, whether complete QE/projection was used, whether critical-value geometry contributed, and result-specific metadata.

## structured_proof_diagnostics

`structured_proof_diagnostics(result)` **Return**s a `ProofDiagnostics` value that normalizes heterogeneous result metadata into one proof trace. It prefers an object's explicit `proof_trace`; otherwise it derives conservative steps from method/diagnostic metadata. It is diagnostic evidence, not an independent certificate verifier.

```python
from semialg import structured_proof_diagnostics

diagnostics = structured_proof_diagnostics(angular)
assert diagnostics.used_complete_qe
assert "sphere_constraint" in diagnostics.steps
```

## Invariants useful for downstream local analysis

These APIs are designed around semantic invariants instead of expression shape. Equivalent formulas should describe the same closure and image; invertible affine coordinate changes should transport closure membership and correlated images; renaming symbols must not change the mathematical result; and a certified image must reject points that satisfy independent coordinate bounds but violate the joint relation. The test suite contains explicit metamorphic contracts for these properties.

## Local asymptotic geometry API

The local-asymptotic geometry family keeps exact approach-set reasoning in
Semialg while leaving expansions to downstream packages. `LocalGerm` and
`local_germ` normalize a set at a closure point; `LocalGeometry` and
`local_geometry` combine local dimension with the Bouligand tangent cone.
`LocalSignStratum` and `local_sign_strata` enumerate locally realizable sign
vectors, while `local_components` separates punctured branches.

`CurveSelectionWitness` and `curve_selection` first seek an explicit certified
polynomial approach curve and otherwise return a theorem-backed existence
certificate without pretending that an explicit parametrization was found.
`LocalRangeResult`, `LocalBoundCertificate`, `local_range`, and `local_bound`
provide exact fixed-ball ranges and can find a positive radius for a requested
or automatically chosen finite local bound; the concrete bound is independently
checked before certification. `ParameterStratum` and `parameter_strata` expose
exact parameter feasibility regions.

For relative rates, `vanishing_order` and `contact_order` operate on explicit
curve witnesses. `BlowupChart` and `blowup_charts` construct ordinary or
weighted polynomial blow-up charts. `local_image` and `local_preimage` transport
germs through exact semialgebraic maps, and `path_independent` certifies that a
value holds on some sufficiently small punctured neighborhood of the germ.


## Local germs, signs, and components

`local_germ(region, point, variables=None)` normalizes a semialgebraic set at a
closure point. `local_geometry` adds the exact local dimension and Bouligand
tangent cone. `local_sign_strata` incrementally enumerates sign vectors whose
strata remain incident at the point.

`local_components` computes connected components after restricting the
punctured set to a sufficiently small neighborhood chosen from exact
critical-distance information. This prevents branches from being merged only
because they reconnect farther from the base point. Singular multivariate
connectivity can still require expensive CAD and is covered by slow regression
tests.

## Curve selection and relative rates

`curve_selection` first searches for a simple explicit polynomial approach
curve. If that search does not construct one, it can return a theorem-backed
existence certificate marked as non-explicit. Operations that need an actual
pullback, including `vanishing_order` and `contact_order`, require an explicit
witness.

## Local bounds

`local_range` certifies extrema on an explicit local ball. `local_bound` can
also search for a positive neighborhood radius and returns a certificate whose
concrete radius/bound pair is checked for a violating point before it is
reported as certified.

## Parameter strata

`parameter_strata(formula, variables, parameters)` returns the exact feasible
and infeasible parameter regions. Passing `conditions_by_value` instead
classifies caller-defined semialgebraic invariants. Cases are made disjoint in
insertion order.

## Local maps and path independence

`local_image` transports a germ through an exact semialgebraic image.
`local_preimage` returns a single source germ only after exact existence and
uniqueness checks on the fiber. `path_independent` certifies equality throughout
some sufficiently small punctured neighborhood.

`blowup_charts` supplies ordinary and weighted polynomial charts for downstream
algorithms that need to expose relative approach scales.
