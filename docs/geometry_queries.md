# Geometric queries

The high-level geometry layer composes CAD/QE, exact optimization, and local polynomial algebra rather than introducing separate numerical approximations.

## Projection, image, preimage, and fibers

```python
import sympy as sp
from semialg import semialgebraic_projection, semialgebraic_image, semialgebraic_preimage, fiber

x, y, u, a = sp.symbols("x y u a", real=True)

semialgebraic_projection((x >= 0) & (x <= y) & (y <= 2), [x], [x, y])
# (y >= 0) & (y <= 2)

semialgebraic_image(x, (x >= -1) & (x <= 2), [x], image_variables=[u])

# Parameters remain free when source variables are explicit:
a = sp.symbols("a", real=True)
semialgebraic_image(x, (x >= 0) & (x <= a), [x], image_variables=[u], parameters=[a])
# (u >= -1) & (u <= 2)

semialgebraic_preimage(x**2 + y**2, u <= 1, [x, y], target_variables=[u])
# x**2 + y**2 <= 1

fiber((x >= 0) & (x <= a), {a: 2})
# (x >= 0) & (x <= 2)
```

Projection and image are genuine existential QE operations. Preimage and fiber are exact symbolic substitutions.

## Bounding boxes, distances, and critical values

`bounding_box` performs coordinate-wise exact optimization. Distance queries minimize squared Euclidean distance and preserve attainment through result objects. `critical_values` reports exact objective values from isolated KKT/singular candidates, certifies constant values on positive-dimensional connected KKT components when possible, and includes attained extrema.

## Convexity and topology

`is_convex` first recognizes affine/polyhedral intersections and basic polynomial intersections whose sublevel/superlevel Hessians are globally PSD/NSD. Polynomial-matrix semidefiniteness is certified by all principal minors, with complete CAD used only for unresolved minor signs in the original ambient variables. If this proof does not apply, `is_convex` decides the full quantified segment definition of convexity. `is_path_connected` uses the semantic CAD adjacency graph. `path_between` distinguishes explicit 1-D paths from higher-dimensional CAD cell/connector chains.

`euler_characteristic(..., compact_support=True)` uses the selected CAD cells and the additive cell formula. This is the compactly-supported semialgebraic Euler characteristic; it equals ordinary Euler characteristic on compact sets.

## Singular loci and tangent geometry

For algebraic varieties, `singular_locus` uses the Jacobian rank criterion and `tangent_space` is the Jacobian nullspace at a point. `tangent_cone` computes the exact ideal-theoretic cone, not merely the lowest-degree parts of the supplied generators. After translating the point to the origin it forms the m-adic deformation, saturates by the deformation parameter, and takes the special fibre at zero. Thus cancellations among generators are included exactly; for example `⟨x**2+y**3, x**2-y**3⟩` has tangent-cone ideal `⟨d_x**2, d_y**3⟩`.

### Symbol identity and points

String variable names and mapping keys are resolved against the symbols already present in the formula. If two SymPy symbols have the same printed name but different assumptions, the API raises an ambiguity error instead of guessing. Geometry APIs share one point-normalization path for this behavior.


## Derived set relations and geometric properties

The convenience predicates are exact compositions of the existing decision, topology, and boundedness layers:

```python
from semialg import (
    is_empty, is_bounded, is_compact, is_open, is_closed,
    is_subset, is_equal, is_disjoint, intersects, contains_point,
    is_connected, is_full_dimensional, has_empty_interior,
)
```

`is_connected` and `is_path_connected` agree for semialgebraic sets. `connected_components` uses the semantic CAD adjacency graph and returns one formula per certified component rather than merely returning top-level Boolean disjuncts.

## Extrema and level sets

`argmin_set` and `argmax_set` return the *entire* global optimizer set, not just isolated witness points. This matters for positive-dimensional extrema such as `argmin_set(x**2, box, [x, y])`, whose answer contains the complete slice `x = 0`. `extrema_set` is the union of the global minimum and maximum sets.

`level_set`, `sublevel_set`, and `superlevel_set` are exact constructors that can be composed with measure, topology, projection, and optimization operations. Strict sub/superlevel variants are selected with `strict=True`.

## Metric and convex-geometric conveniences

`coordinate_range` delegates to `function_range`. `nearest_point` and `closest_points` expose the optimizing witnesses from the exact distance routines. `diameter` maximizes squared pairwise distance. `support_function(region, direction)` computes `sup(direction·x)` and `width` computes the corresponding max-minus-min directional width. These operations inherit the exact optimization backend's attainment and unboundedness semantics.

## Moments and rigid/affine operations

`centroid` and `covariance_matrix` are concise aliases over the existing exact region-moment machinery. `moment_matrix` returns the normalized raw second-moment matrix, while `inertia_tensor` integrates `||x||^2 I - x x^T` at unit density about the origin.

`translate` and nonzero scalar `scale` are direct substitutions. `linear_image`, `affine_transform`, and `minkowski_sum` are exact semialgebraic images; the latter constructs the existential image of `(a, b) -> a + b`. General instances may therefore require full QE. `squared_distance_range` gives the exact polynomial pairwise-distance-squared range, while `distance_set` converts it to the nonnegative Euclidean-distance set via `d**2`.

## Algebraic-geometry conveniences

`is_singular(equations, point)` evaluates the Jacobian singular-locus criterion at a point, `is_smooth(equations)` decides whether the real singular locus is empty, and `tangent_dimension` returns the Zariski tangent-space dimension. They are thin convenience APIs over `singular_locus` and `tangent_space`.


### Irreducible equations and real connectivity

Polynomial irreducibility is an algebraic property and does not imply that the real zero set is connected. For example, the irreducible hyperbola `x*y - 1 = 0` has two real connected components, and singular irreducible curves can contain isolated real components. Connectivity decisions therefore remain CAD/topology questions even when a defining polynomial has only one algebraic factor.
