# Algebraic roots and exact finite solving reference
## Family contract

**Mathematical return.** Algebraic APIs represent, isolate, compare, classify, and solve exact real algebraic roots and finite polynomial systems.

**Exactness and certification.** Root identity/order and algebraic signs use exact arithmetic, isolating information, or number-field/RUR machinery rather than fixed-precision comparison.

**Algorithm.** Facilities include root isolation/classification, rational univariate representation, subresultant and border-basis machinery, and cached exact algebraic comparisons.

**Complexity and limitations.** Degree growth and coefficient growth can dominate. Positive-dimensional algebraic sets generally require region/CAD or algebraic-geometry APIs instead of finite RUR solving.




## Primary API overview

This table is the substantive coverage target for the primary APIs assigned to this reference page. Each entry states the API's primary role; the family contract and detailed sections below explain shared algorithms, exactness guarantees, and limitations. It is maintained together with `docs/reference/primary_api_manifest.toml`, and documentation tests require every root-level primary API to map here rather than merely appearing in the generated public index.

| API | Kind | Role / return |
|---|---|---|
| `classify_real_roots` | function | Classify real roots of a univariate polynomial or polynomial family. |
| `associated_primes` | function | Return certified associated primes from a primary decomposition. |
| `certified_independent_localization` | function | Choose a maximal independent set and certify zero-dimensional localization over QQ(U). |
| `certified_radical_minimal_prime_decomposition` | function | Compute a certified radical/minimal-prime decomposition when the exact engine completes. |
| `compress_primitive_element` | function | Compress an algebraic function-field tower to one certified primitive extension. |
| `contract_localized_ideal` | function | Contract an ideal from QQ(U)[X] by canonical denominator clearing and saturation. |
| `maybe_compress_primitive_element` | function | Apply primitive-element compression only when the exact cost heuristic predicts a benefit. |
| `primary_decomposition` | function | Compute certified primary components in the supported exact regimes. |
| `radical_ideal` | function | Compute and certify generators of the radical ideal when decomposition completes. |
| `recursive_regular_chain_decomposition` | function | Compute a proof-producing recursive squarefree regular-chain decomposition. |
| `saturation_stabilization` | function | Compute the exact exponent where I:h^m stabilizes and verify the GTZ split identity. |
| `verify_independent_localization_certificate` | function | Replay a maximal-independent-set localization certificate exactly. |
| `verify_localization_contraction_certificate` | function | Replay localized ideal contraction and re-extension exactly. |
| `verify_minimal_prime_decomposition_certificate` | function | Replay a radical/minimal-prime decomposition certificate. |
| `verify_primary_decomposition_certificate` | function | Replay primaryness, reconstruction, and associated-prime claims exactly. |
| `verify_primitive_element_compression` | function | Replay primitive-element guard and field-isomorphism data exactly. |
| `verify_radical_ideal_certificate` | function | Replay radical reconstruction and containment claims exactly. |
| `verify_regular_chain_decomposition_certificate` | function | Replay the recursive regular-chain branch proof without rerunning search. |
| `verify_saturation_stabilization_certificate` | function | Replay the colon stabilization chain and GTZ reconstruction identity exactly. |
| `verify_triangular_primality_certificate` | function | Replay successive fraction-field irreducibility/primality stages exactly. |
| `zero_dimensional_primary_decomposition` | function | Compute GTZ3 primary components over QQ or a certified finite algebraic extension. |
| `verify_zero_dimensional_primary_certificate` | function | Replay coefficient-field validity, maximal-prime separation, saturation, and transport exactly. |
| `ThomEncoding` | class | Exact derivative-sign certificate identifying one real root of a univariate polynomial. |
| `thom_encoding` | function | Construct the Thom encoding of a specified exact real root. |
| `thom_encodings` | function | Enumerate all distinct real roots by exact Thom encodings in increasing order. |

## `root_of` and `AlgebraicRootFunction`

`root_of(polynomial, variable, index)` denotes an ordered real algebraic root with binder-aware substitution semantics. `AlgebraicRootFunction` carries the root identity as base parameters vary and supports certified specialization, comparison, derivatives, and regularity checks in supported cases.

When a fully specialized low-degree root has a uniquely certified exact radical presentation, semialg may display that radical while retaining the ordered-root identity as the certification basis.

## `classify_real_roots`

```text
classify_real_roots(polynomial, variable, *, parameters=None)
```

Classifies real-root behavior, including parameter-dependent cases. String variables/parameters resolve to the actual expression symbols.

## `solve_zero_dimensional_system`

```text
solve_zero_dimensional_system(
    equations, inequalities=None, variables=None, *,
    backend="rur", real=True, parameter=None,
    max_separating_attempts=64, return_result=False
)
```

Solves supported finite polynomial systems exactly using rational-univariate representation. The default return is the tuple of exact solution-coordinate tuples in the requested variable order. Optional inequalities are exact filters on algebraic candidate points. Set `return_result=True` to retain backend/RUR metadata in `ZeroDimensionalSolveResult`.

## Thom encodings and RUR sign determination

A `ThomEncoding` records the signs of the successive derivatives of a univariate defining polynomial at one exact real root. `thom_encoding` constructs one certificate and `thom_encodings` enumerates the distinct real roots in exact order. This gives a root identity that does not depend on a decimal approximation.

`RationalUnivariatePoint.thom_encoding` applies the same representation to the parameter root of an existing rational-univariate representation. `RationalUnivariatePoint.sign_of()` reduces a multivariate polynomial modulo the RUR to a univariate polynomial and determines its sign at that Thom-certified root. This is the preferred exact path for repeated sign queries on finite polynomial systems.

## Certified polynomial root intervals and separators

`certify_polynomial_root_interval`,
`verify_polynomial_root_interval_certificate`,
`rational_between_algebraic_reals`, and
`certified_sign_stable_root_neighborhood` provide the advanced exact root-interval
contracts used by CAD lifting and direct algebraic callers. The interval verifier
can optionally bind a certificate to an expected polynomial, variable, interval,
and endpoint convention in addition to replaying its stored proof data.

These APIs accept exact univariate polynomials over `QQ`. For algorithms, endpoint
semantics, examples, and the distinction between distinct roots and multiplicity,
see [Certified algebraic roots](../guides/certified_algebraic_roots.md). For the
separate transcendental solver boundary, see
[Transcendental and algebraic solving scope](../guides/transcendental_scope.md).

## Exact comparisons

Certified algebraic ordering uses isolating intervals, minimal/defining polynomials, exact sign determination, and appropriate algebraic representations. Fixed-precision sorting is not a proof mechanism.

## Advanced algebraic APIs

The package also exposes RUR, border-basis, subresultant, and related algebraic utilities for advanced users. These are lower-level than the primary decision/solve interfaces and may require stronger preconditions.

See [Root classification](../root_classification.md) and [Exactness and certification](../concepts/exactness_and_certification.md).

## RUR coefficient domains

Rational univariate representation supports zero-dimensional systems over
`QQ` and exact algebraic-number extensions such as `QQ<sqrt(2)>` or a common
field generated by several algebraic coefficients. For example,

```python
x, y = sp.symbols("x y", real=True)
solve_zero_dimensional_system_with_rur(
    [x + sp.sqrt(2) * y, y - sp.sqrt(2)],
    [x, y],
)
# ((-2, sqrt(2)),)
```

A coefficient domain containing symbolic parameters or transcendental constants
is not treated as an algebraic number field. The RUR planner declines that
branch cleanly so a higher-level exact planner can choose another backend.
Low-level SymPy domain-conversion exceptions such as `CoercionFailed` are not
part of the public RUR applicability contract.

## Reduced polynomial-locus decomposition

`equidimensional_decomposition(equations, variables)` computes an exact finite
cover of a polynomial zero set using reduced real-set semantics. Multiplicity in
the input ideal does not create additional geometric components, and every
returned `AlgebraicLocusPiece` records both its exact real semialgebraic
dimension and the Krull dimension of the corresponding algebraic branch.

The decomposition engine combines several exact reductions according to the
structure exposed by the normalized ideal:

- monomial Groebner bases are reduced to minimal coordinate primes using minimal
  vertex covers of the monomial supports;
- factorizable Groebner generators are represented by a factor-support
  hypergraph whose inclusion-minimal hitting sets give an exact reduced-set
  cover;
- lexicographic triangular slices are refined through exact initial and separant
  saturations; terminal branches can then be certified by a squarefree regular
  chain.  The chain certificate requires distinct leaders, regular initials and
  separants under exact saturation, and exact equality between the saturated
  chain ideal and the branch ideal.  In characteristic zero this proves the
  saturated ideal radical and unmixed, hence equidimensional;
- Rabinowitsch saturation may split an ideal through
  `V(I) = V(I + <h>) union V(I : h^infinity)` when radical-containment checks
  prove that both branches are strict reduced subloci;
- definite sums of even monomials and positive-semidefinite quadratic
  polynomials with exact minimum zero admit certified real-radical reductions.

These reductions cooperate rather than define separate public modes. For
example, `(f*x, f*z)` is represented as `V(f) union V(x, z)`, while `(f**2,
f*z)` reduces to `V(f)` after contained branches are removed. Over the reals,
`x**2 + y**2 = 0` in three variables is represented by the line `V(x, y)`, and
`x**2*y**2 + z**2 = 0` is represented by `V(x, z) union V(y, z)`. A definite
equation such as `x**2 + y**2 + 1 = 0` is certified empty. When real reduction
lowers dimension, `algebraic_dimension` preserves the Krull dimension of the
pre-reduction branch so the distinction is explicit.

The regular-chain layer also supplies certified recursive minimal-prime
refinement.  Primality is tested in the successive quotient/fraction fields of
the triangular prefix rather than by factoring each polynomial in the ambient
ring.  Reducible leader polynomials produce exact child branches, which are
re-normalized and recursively analyzed before the final radical-intersection
reconstruction check.  The exact coefficient layer supports recursive monogenic extensions over
rational-function fields.  Factorization descends by norms/resultants and
recovers factors by gcd in the extension; quotient-ring/saturation splitting is
retained as a conservative fallback when the norm route cannot decide.

The regular-chain layer substantially broadens general equidimensional
certification beyond hypersurfaces and generator-count complete intersections;
for example, the monomial curve `<x**2-y, x*y-z>` is certified through the
squarefree chain `(y-x**2, z-x**3)` even though its reduced Groebner basis has
three generators at codimension two.  The implementation now also provides a certified primary-decomposition layer.
Prime, principal, monomial, and zero-dimensional regimes use specialized fast
paths; arbitrary positive-dimensional nonmonomial ideals fall through to the
recursive GTZ localization/contraction engine.  A general real-Nullstellensatz
algorithm remains outside this layer.  A positive-dimensional
terminal ideal that cannot be proved equidimensional remains represented exactly with
`equidimensional=False`, and the enclosing decomposition has `complete=False`.
No unsupported case is promoted to a componentwise geometric conclusion.

Every decomposition result carries a `DecompositionCertificate`. The certificate
records the normalized input, canonical reduced-real pieces, dimensions,
equidimensionality claims, the exact reduction families exercised, and the
piece budget. `verify_decomposition_certificate(certificate)` deterministically
replays the decomposition with the same budget and rejects any mismatch or
mutation. Geometry consumers require both `complete=True` and successful
certificate replay before using component dimensions.

This certificate boundary is the contract used by component-relative regularity
analysis. `region_singular_locus` tests each certified piece with its own real
codimension, keeps intrinsic component singularities distinct from
singularities caused by intersections of distinct pieces, preserves reduced-set
semantics for repeated generators, and restricts the result to the realized
semialgebraic boundary where appropriate.

## Associated primes and primary decomposition

`primary_decomposition(equations, variables)` returns exact primary components
only after their intersection is verified to equal the original ideal.  The
current proof-producing algorithms are complete in four important regimes:

- a certified prime ideal is itself a primary component;
- a principal ideal is factored over `QQ`, and each irreducible prime power
  `<p**e>` is certified `<p>`-primary;
- an arbitrary monomial ideal is Artinianized, its maximal standard monomials
  are enumerated exactly, and the corresponding irreducible monomial ideals are
  de-Artinianized and intersected back to the source ideal; this includes
  embedded primes;
- a zero-dimensional ideal with certified maximal minimal primes is decomposed
  by exact localization/saturation separators.  Since an ideal with maximal
  radical is primary, each isolated component receives a theorem-level primary
  certificate.

`associated_primes(...)` returns the distinct radicals of an irredundant
certified primary decomposition.  For example,
`<x**2, x*y> = <x> intersect <x**2, y>`, so both `<x>` and the embedded prime
`<x, y>` are returned.  `verify_primary_decomposition_certificate` rechecks the
primary criterion for every component, exact ideal reconstruction, dimensions,
degrees, distinct radicals, and irredundancy; `replay_certificate` dispatches to
that verifier without trusting the search result.

For GTZ localization, `zero_dimensional_primary_decomposition(generators, variables,
coefficient_field=...)` supplies the coefficient-field-generic zero-dimensional
stage.  Over `QQ` it works directly.  For a finite monogenic tower over `QQ`,
semialg first certifies every defining polynomial irreducible over the preceding
field, lifts the ideal to `QQ[alpha_1, ..., alpha_r, X]`, certifies the maximal
minimal primes there, isolates each local component by an exact separator and
saturation, and transports the component back to the requested coefficient
field.  Because a zero-dimensional quotient is Artinian, all of its prime ideals
are maximal; the isolated local component is therefore primary.
`verify_zero_dimensional_primary_certificate` independently rechecks the field
tower, minimal-prime certificate, separators, saturation certificates, exact
intersection reconstruction, and lift/transport identities.

For arbitrary positive-dimensional nonmonomial ideals, the recursive GTZ driver
now localizes at a certified maximal independent set `U`, decomposes the finite
`QQ(U)`-algebra intrinsically, contracts its local primary factors, and recurses
on the exact companion branch from the saturation identity.  The localized
Artinian decomposition is certified from multiplication matrices: the trace-form
rank equals the dimension of the reduced algebra in characteristic zero, and a
deterministically searched linear form is accepted only when the squarefree
degree of its characteristic polynomial reaches exactly that rank.  Its
irreducible factors then separate the maximal factors of the reduced quotient.
The trace-form kernel supplies the nilradical, so the corresponding localized
maximal radicals and primary factors are reconstructed exactly.

Contraction back to `QQ[U,X]` uses the existing monic-local-basis denominator
certificate.  To recover the residual branch, GTZ computes for each contracted
hull generator `g` a nonzero parameter annihilator in `(I:g) intersect QQ[U]`;
their product `h` is then certified by exact colon stabilization.  Replay checks
`I:h^infinity` equals the localized hull and verifies
`I = (I:h^infinity) intersect (I + <h^m>)`.  The residual branch is required to
make strict ideal progress and never increase dimension; it need not lower
dimension at every step, which is essential for examples such as `<x*y>`.

`verify_gtz_primary_decomposition_certificate` replays the recorded tree rather
than repeating search.  It rechecks every localization certificate, trace-form
Artinian split, primary/radical contraction, saturation exponent, recursive
companion branch, final same-prime merges, irredundancy, dimensions/degrees, and
exact intersection reconstruction.

## Recursive regular chains and Hilbert invariants

`recursive_regular_chain_decomposition(equations, variables)` is the exact
triangular decomposition layer used to strengthen positive-dimensional
algebraic reasoning.  It is obstruction-driven: when squarefree-regular-chain
certification fails, the recursion identifies the first prefix where an initial
or separant is a zero divisor, and uses that exact condition for a certified
closed/saturation split.  Competing equations with the same leader are analyzed
by subresultant PRS data, principal coefficients, defective-entry coefficients,
and bidirectional pseudo-remainders; this regular-gcd information is collected
from both the lex basis and the current branch generators so Gröbner reduction
does not hide a useful split.  A systematic saturation fallback also considers
triangular initials, separants, lower-variable coefficients, and their exact
squarefree factors.  Quotient-field factorization remains the final algebraic
factor splitter.  Every proposed condition is accepted only after proving
``V(I) = V(I + <h>) union V(I : h**inf)`` with two strict children.  A terminal
component is marked complete only when exact saturation proves that its
triangular equations form a squarefree regular chain.  Prefix regularity is
checked relative to the prefix already saturated by every previous initial.

The result also carries independent degree accounting.  `source_degree` is the
scheme-theoretic degree of the input ideal; `reduced_degree` is the degree of
the exactly reconstructed reduced union.  These can differ for inputs with
multiplicity, e.g. `x**2*y` has source degree 3 and reduced degree 2.
`degree_complete=True` means the sum of degrees of the top-dimensional regular
chain components equals the Hilbert degree of that reconstructed reduced
union.

The Hilbert infrastructure is available from `semialg.algebraic` (and the
`semialg.algebraic_geometry` module):

```python
ideal_hilbert_data(equations, variables)
hilbert_function(equations, variables, degree)
hilbert_polynomial(equations, variables)
ideal_degree(equations, variables)
```

These invariants are computed independently from a degree-compatible exact
Groebner basis.  The leading monomial ideal has the same Hilbert function as
the source ideal; its Hilbert-series numerator is obtained by exact
inclusion--exclusion on minimal monomial generators.  Cancellation at `t=1`
then gives Krull dimension and affine multiplicity/degree.  No regular-chain
assumption is used in this calculation.

## Certified modular Groebner acceleration

For polynomial ideals over `QQ`, the shared exact Groebner constructor can use
multimodular reconstruction before falling back to SymPy's direct rational
algorithm.  `modular_groebner_basis_qq(generators, variables, order=...)`
computes reduced bases over a deterministic sequence of finite fields, rejects
unlucky images whose leading-monomial/support signature disagrees, combines
matching coefficients with CRT, and applies rational reconstruction as the
modulus grows.

For localized GTZ arithmetic the same idea now extends to rational-function
coefficient fields.  `modular_groebner_basis_fraction_field(generators,
variables, parameters, order=...)` computes images over
`GF(p)(parameters)`, normalizes each rational-function coefficient by a monic
denominator, rejects changing numerator/denominator supports, reconstructs both
numerator and denominator coefficients with CRT plus rational reconstruction,
and then verifies the reconstructed basis exactly over `QQ(parameters)`.
The shared Groebner dispatcher attempts this path automatically for eligible
fraction-field domains, so GTZ localization and trace-algebra work benefit
without changing their APIs.

The finite-field result is never a proof.  A reconstructed candidate is
accepted only after all of the following exact checks over `QQ` succeed:

1. recomputing the Groebner basis of the candidate leaves the candidate
   unchanged (Buchberger/reduced-basis verification);
2. every original generator reduces to zero by the candidate, proving
   `I <= <G>`; and
3. exact Macaulay linear algebra finds polynomial multipliers expressing every
   element of `G` in the original generators, proving `<G> <= I`.

The resulting `ModularGroebnerCertificate` stores those exact membership
representations and can be replayed with
`verify_modular_groebner_certificate`.  If rational reconstruction, support
stabilization, or exact membership certification is inconclusive, the shared
Groebner helper transparently uses the ordinary exact path instead.  Modular
computation therefore changes performance only, never the mathematical
contract.


## Certified modular resultants and subresultants

The same proposal/reconstruction/proof boundary is used for elimination
polynomials.  `modular_resultant_qq(f, g, x)` supports coefficients in
`QQ[parameters]`: finite-field images are computed over
`GF(p)[parameters][x]`, images whose eliminated-variable degrees or polynomial
supports change are treated as unlucky, and matching coefficients are merged
by CRT before rational reconstruction.

A reconstructed resultant is *not* verified by simply recomputing the symbolic
resultant.  If `m = deg_x(f)` and `n = deg_x(g)`, each parameter `a` has the
conservative exact degree bound

```text
deg_a Res_x(f,g) <= n*deg_a(f) + m*deg_a(g).
```

`ModularResultantCertificate` records these bounds and a deterministic tensor
grid of `bound + 1` rational values per parameter.  At every grid point the
verifier computes the fixed-`m,n` Sylvester determinant exactly over `QQ`.
Equality on the full grid proves polynomial identity by interpolation
uniqueness, even at points where specialization lowers the apparent degree in
`x`.  Large grids are rejected and transparently fall back to direct exact
resultant computation.

`modular_subresultants_qq(f, g, x)` reconstructs the entire native
subresultant sequence coefficientwise with the same lucky-prime/CRT/rational
reconstruction discipline.  Its certificate additionally verifies the
resultant certificate and replays the exact characteristic-zero subresultant
normalization.  This last check is stronger than verifying a
pseudo-remainder recurrence: a recurrence identifies a PRS only up to scale,
whereas callers use the principal subresultant coefficients and therefore
need the exact normalization.

`subresultant_prs` and the CAD resultant helpers attempt these certified
modular paths automatically for eligible rational polynomial input and use the
existing exact implementation whenever modular reconstruction or certification
is inconclusive.

### Exact algebra performance

Primitive-element compression now avoids symbolic determinant expansion in the
primitive parameter once the compositum has degree greater than three.  If the
power-basis matrix has size `d`, its determinant has primitive-parameter degree
at most `d(d-1)/2`; semialg evaluates the exact determinant at `bound + 1`
integer parameters, reconstructs the whole guard by interpolation over the
rational-function base, and verifies the reconstruction at one additional
point.  This is an exact reconstruction algorithm, not a numerical heuristic.

GTZ independent-set selection now scores every certified maximal independent
set by estimated parameter-degree and coefficient occurrence, preferring the
localization expected to cause less `QQ(U)` coefficient growth.  Regular-chain
PSC/saturation splitter proposals are likewise ordered by a structural cost
score favoring low-degree, sparse conditions involving fewer variables.  In
both cases the planner only changes proposal order: the existing exact
independence, zero-dimensionality, saturation, and reconstruction checks remain
mandatory.

`scripts/benchmark_exact_algebra_p2.py` times direct versus certified modular
`QQ(U)` Groebner computation and primitive-element guard reconstruction on a
deterministic corpus.  Benchmark cases are accepted only after the resulting
certificate/compression replay succeeds.

### GTZ execution planning and caching

The recursive GTZ driver uses the same exact certificate contract throughout,
but its expensive canonical `QQ` Groebner bases now pass through semialg's
certified modular dispatcher.  Finite-field work is therefore only a proposal:
the reconstructed basis is accepted only after the existing exact Buchberger
and two-way ideal-membership certificate succeeds.  Inconclusive modular work
falls back to direct exact arithmetic.

Repeated canonical bases over both `QQ` and localized fraction fields `QQ(U)`
are held in bounded process-local LRU caches.  Recursive GTZ nodes additionally
memoize canonical source ideals during a single decomposition.  These caches
store exact canonical algebraic objects, not truth values or unverified search
results, so certificate replay and final reconstruction remain exact.

`plan_gtz_primary_decomposition(...)` exposes a deterministic execution plan
with generator count, Krull dimension, degree when available, maximum input
degree, and acceleration policy.  The planner changes cost policy only; it
cannot bypass primary/reconstruction verification.  `gtz_cache_info()` and
`clear_gtz_caches()` support profiling and long-running-process hygiene.

`scripts/benchmark_gtz.py` provides a deterministic corpus covering isolated,
embedded, localized, and three-variable cases.  It times the computation but
accepts a benchmark result only when the resulting GTZ certificate replays.
The same printed examples can be used for external comparison with Singular's
`primdecGTZ`; Singular is not required at runtime.


## Polynomial-map implicitization and Zariski closures

`implicitize_polynomial_map(mapping, parameters, ...)` computes the exact
elimination ideal of a polynomial map graph.  The optional `domain_equations`
restrict the parameter space algebraically; semialgebraic inequalities are not
accepted because their effect on the Zariski closure requires a separate
density argument.  `return_result=True` exposes both the graph Groebner basis
and the reduced image equations.

`zariski_closure(mapping, parameters, ...)` presents the same exact elimination
as a conjunction of polynomial equalities in the image coordinates.  This is
the algebraic closure of the polynomial image, not the exact real
semialgebraic image; use `region_image` when inequalities and real image
semantics must be retained.

The result of `implicitize_polynomial_map(..., return_result=True)` can be
passed directly to `singular_locus`, which applies the Jacobian criterion in
the implicit image coordinates.  For example, `(t**2, t**3)` implicitizes to
the cusp and its singular locus contains the origin.


## Certified reduced varieties and component singularities

`certified_radicalization(equations, variables)` computes `sqrt(I)` through the certified regular-chain radical engine and replays its certificate before exposing the reduced generators. `irreducible_components(...)` is stricter: it returns components only when the radical cover is certified complete and every irredundant component has a primality certificate, so the result is a complete minimal-prime decomposition.

`singular_locus(...)` now radicalizes by default before applying the Jacobian criterion. Thus `x**2 = 0` and `x = 0` define the same smooth reduced line. `reduced_component_singular_loci(...)` computes the Jacobian singular locus on each certified irreducible component independently, distinguishing singularities intrinsic to a component from intersections between otherwise smooth components.

## Stratified singular geometry

`MinimalPrimeIntersection`, `LocalDimensionStratum`, `SingularGeometryStratum`, and `StratifiedSingularGeometry` are the structured results for this layer. `minimal_prime_intersections()` computes exact intersections of certified minimal-prime components by ideal sums. `local_dimension_strata()` returns constructible loci on which local algebraic dimension is constant. `stratified_singular_geometry()` combines those results with each reduced component's intrinsic Jacobian singular locus and refines component membership into disjoint constructible singular strata. Component crossings are represented separately from intrinsic singularities.

## Local branch and reduced variety geometry

`ReducedAlgebraicVariety`, `IrreducibleAlgebraicComponent`, and
`ReducedComponentSingularLocus` make the reduced algebraic variety explicit before
singularity analysis. `PolynomialMapImplicitizationResult` and `ZariskiClosureResult`
record exact elimination results without identifying a Zariski closure with a real
semialgebraic image.

`LocalBranchGeometry`, `BranchTangentGeometry`, and `ComponentIntersectionGeometry`
separate intrinsic branch singularities from intersections of distinct minimal-prime
components. `local_branch_geometry()` computes tangent-cone and tangent-space data,
tangent-dimension excess, local multiplicity/degree information, and certified
transversality information at the supplied exact point.
