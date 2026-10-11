# Chart geometry certificates

The specialist `semialg.chart_certificates` module proves geometry without making
claims about integrability. `Chart` describes real coordinates, forward and inverse
formulas, inner-first ranges, optional null exclusions, a source branch domain,
and a rational/algebraic/canonical-polar kind. Coordinates must be plain real
symbols; sign or arithmetic flags on coordinates cannot stand in for proved
domain restrictions. Such candidates decline, while parameter facts remain
valid in their own scope.

`certify_chart(chart, level_map, variables, assumptions=..., budget=...)` checks
inverse and scalar/vector level identities, full rank, the coarea Gram metric identity and regular
source domains. `certify_cover(region, certificates, ...)` replays certificates,
checks source coverage modulo certified null sets and disjoint chart interiors.
`certify_null_set` only accepts proved subsets of polynomial zero sets whose
polynomials are nonzero on every allowed parameter fiber. Equality fibers cannot
be discarded merely because a coefficient is generically nonzero.

`certify_linear_metric` proves full column rank by a positive Gram determinant;
it serves both square affine maps and intrinsic affine embeddings. Affine region
transformations and intrinsic affine integration use this shared certificate.
All records are frozen; successful proofs use bounded caches keyed by candidates,
explicit assumptions and `GeometryBudget`. Temporary global SymPy assumptions
are ignored by certification; parameter-cell facts must be supplied explicitly. `clear_chart_cache` and `chart_cache_info` expose
cache control. Failed proofs raise `GeometryProofError`; malformed requests raise
`ValueError`. Certificate records are geometry evidence, not mathematical-result
wrappers and not trusted caller flags.

```python
import sympy as sp
from semialg.chart_certificates import Chart, certify_chart, certify_cover, polar_chart
x, y, t, u = sp.symbols("x y t u", real=True)
chart = Chart((t, u), (t*u, u), (x/y, y), ((u, 0, 1), (t, 0, 1)),
              chart_exclusions=sp.Eq(u, 0), source_exclusions=sp.Eq(y, 0))
certificate = certify_chart(chart, x/y, (x, y))
region = sp.And(x >= 0, x <= y, y >= 0, y <= 1)
coverage = certify_cover(region, (certificate,))
assert coverage.source_null.zero_polynomials == (y,)

polar = polar_chart((x, y), (t, u), (0, 1), squared=True)
polar_proof = certify_chart(polar, x*x+y*y, (x, y))
assert certify_cover(x*x+y*y <= 1, (polar_proof,)).charts == (polar_proof,)
```

Algebraic charts support selected positive rational-power branches with checked
inverse identities. Canonical planar polar charts have an explicit center/seam
null certificate and injectivity theorem for the open angular interval. Arbitrary
transcendental charts and automatic nonlinear critical-set stratification are
outside the implemented scope. Null sets may still carry nonintegrable singular
behavior in a density; the integration caller must prove convergence separately.

Regression: `tests/test_chart_certificates.py`.

## Vector maps and parameter events

For a tuple of `k` level maps, `1 <= k < m`, the first `k` chart coordinates
are levels. The map Jacobian squared is `det(Dg*Dg.T)`; the fiber Gram
determinant uses the remaining `m-k` forward-chart columns. Inverse and all
level identities, full rank and the volume/metric identity are replayed just as
for scalar charts. Canonical polar charts remain scalar.

`chart_parameter_events(chart, level_map, variables, parameters, budget=...)`
returns an immutable tuple of bounded parameter-only rank, pole and interval-width
factors. They are candidate transitions, not certificates that a generic formula
holds on equality fibers. Integration callers must validate partitions and
specialize those fibers. The event discovery uses the same degree and expansion
guards before determinant work.

`parameter_geometry_events(expressions, parameters, budget=...)` supplies the
same bounded factoring for shell endpoints, affine rates and coincident level
corners that do not need a supplied chart. Determinant proofs estimate symbolic
product growth before expansion; triangular determinants keep their cheap exact
product rule. These guards may decline a mathematically valid chart without
changing its claimed geometry.
