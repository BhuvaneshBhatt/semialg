# Limitations

`semialg` is conservative by design. Exact APIs should return certified results for supported cases and decline unsupported or uncertifiable cases rather than silently promote a fixed-precision approximation into an exact symbolic claim.

## Fundamental computational limits

- CAD/QE has doubly exponential worst-case complexity and can become expensive after only a few additional variables.
- Variable ordering can change CAD cost dramatically.
- Gröbner, resultant, RUR, tangent-cone, and active-set computations can also grow rapidly.
- General transcendental quantifier elimination lies outside real-closed-field CAD.
- Exact geometric decomposition can succeed even when SymPy cannot express the resulting integral in closed form.

## CAD and quantifier elimination

The complete path uses structural presolve, safe affine substitution, exact Fourier–Motzkin elimination, quantifier-aware variable ordering, and specialized exact backends before general CAD. Reduced/equational-constraint CAD is used only when its side conditions are certified; otherwise the complete path remains authoritative.

Pilot lifting is a bounded cost probe, not a proof procedure. It deliberately avoids full CAD sign/provenance work. A cheap heuristic score can therefore mispredict runtime without affecting correctness.

## Parameter-dependent results

Parameterized solvability, root counts, optimization, ranges, and selected integration problems can return exact guarded strata. Arbitrary-degree univariate root-count families use subresultant/Sturm sign stratification when an exact parameter CAD can be constructed; low-degree closed forms remain preferred for readability and speed. General multidimensional nonlinear parameter-dependent CAD bounds, intrinsic parameterized geometry, and compact closed-form reconstruction are not complete.

An exact result may remain a quantified first-order relation. Quantifier-free presentation is not itself part of the exactness guarantee.

## Algebraic root functions

`AlgebraicRootFunction` supports exact specialization, certified same-stack comparison, implicit differentiation, and regularity checks. Comparisons without sufficient shared certification may remain unknown. Radical forms are presentation conveniences unless branch identity is certified.

## Intrinsic integration

Regular triangular CAD graph strata can be integrated with the induced Hausdorff metric. General singular/non-graph stratification, multiple-chart geometry, and arbitrary singular manifolds remain incomplete. Uncertified algebraic sections are not silently treated as regular.

## Formula simplification

The supported polynomial fragment has a stable deterministic normal form: polynomial atoms are primitive and consistently oriented, repeated zero-set multiplicities are removed, scalar bounds and Boolean structure are normalized to a fixed point, and guarded shared-CAD implication checks remove provable redundant literals/branches. This makes repeated simplification idempotent for the supported fragment and greatly reduces representation-dependent output. The package does not claim a globally smallest formula under every possible Boolean/algebraic cost model; semantic equality should still be checked with `equivalent`/region equality predicates when required.

## Polynomial optimization

The optimizer handles many exact polynomial problems through structural reduction, KKT/active-set systems, singular and positive-dimensional loci, RUR solving, exact comparison, and CAD certification. Hard cases include high-dimensional open or unbounded regions, large Boolean expansions, and problems with no useful finite candidate structure.

## Region integration

Ambient-measure integration supports standard regions and typed CAD-cell bounds, including disjoint decomposition of composite Boolean regions. It no longer assumes a fixed vertical orientation: a bounded coordinate-order search can select simpler certified CAD bounds, including arbitrary-dimensional full-dimensional CAD cells. Symbolic antiderivative evaluation may still remain unevaluated even when the geometric decomposition is exact. Singular intrinsic strata and general nonlinear multidimensional parametric bounds may require unsupported decomposition.

## Topology and geometry

- Connected components, connectivity/path-connectivity decisions, Euler characteristic, and certified CAD connector chains are available for supported semialgebraic sets.
- `path_between` does not construct a general Canny/Basu–Pollack–Roy algebraic roadmap parameterization in dimensions above one.
- Homology and Betti numbers are not implemented. `CADCellComplex` provides unsigned codimension-one incidence, not oriented cellular boundary operators.
- The optimization backend can analyze a positive-dimensional projected KKT or singular active-boundary locus through exact function-range/image-CAD computation when recursive coordinate elimination stops, and can recover an exact witness for an attained endpoint. The higher-level `critical_values` API still does not return a complete symbolic decomposition of every varying critical-value image.
- `SemialgebraicRegion.singular_locus()` uses reduced-real component-relative Jacobian semantics and restricts inequality contributions to residuals realized on the actual boundary. Use `nonsmooth_locus()` for transverse corners/ridges on multi-active boundary strata and `active_boundary_strata()` for exact realized defining-inequality activity strata; these are not full Whitney stratifications.
- CAD-cell triangulation is numerical. Conforming assembly checks sampled common-face triangulations, but ambient isotopy to the exact semialgebraic set is not certified.
- Generic convexity can still be expensive because the complete fallback is a quantified semialgebraic problem with duplicated ambient variables and a segment parameter. The staged convexity engine handles complete one-dimensional connectivity, affine/polyhedral sets, basic quadratics, global and domain-relative polynomial Hessian certificates, selected topology-based rejections, and exact segment witnesses before that fallback. Representation-dependent cases that evade these certificates may still require full QE.
- Derived operations such as singular affine images, Minkowski sums, distance sets, diameter, and support functions may reduce to full QE or exact global optimization.

## Exact sampling and numerical presentation

Exact witness generation may use numerical evaluation to **propose** rational candidates, but acceptance is based on exact comparison or exact formula validation. In exact mode, interval membership and bound certification do not depend on binary-float equality or a fixed decimal epsilon. Explicit plotting, discretization, and `exact=False` sampling are numerical by contract.

## Failure and unsupported cases

Expected strategy failure, malformed input, certification failure, and resource limits have distinct exception/result paths. Non-core helpers may use conservative fallback behavior, but unexpected implementation failures should not be converted into mathematical claims.

See [Errors and failure modes](guides/errors_and_failure_modes.md) for user-facing behavior and [Exactness and certification](concepts/exactness_and_certification.md) for the proof/representation distinction.


## Limitation types and extensibility

The remaining limitations fall into different categories and should not be interpreted as equally tractable:

| Limitation | Category | Near-term status |
|---|---|---|
| General CAD/QE scaling | Fundamental/algorithmic | Incremental planning, equational-constraint, decomposition, and cache improvements only; no general polynomial-time fix is expected. |
| Semialgebraic function graph coverage | Missing implementation | Actively extensible. Current exact coverage includes nested rational/algebraic graphs, `Abs`, `sign`, `Heaviside`, `Min`, `Max`, finite `Piecewise`, and broader Boolean formula composition. |
| Parameterized real-root counts above quartic | Engineering/algorithmic | Improved: general subresultant/Sturm sign stratification is available, with conservative fallback when exact CAD construction fails. |
| Closed-form antiderivatives after exact geometric decomposition | External symbolic-integration limitation | Exact geometry can be retained even when presentation as an elementary closed form fails; broader exact-unevaluated and certified-numerical result modes are feasible extensions. |
| Positive-dimensional critical-value images | Partial implementation | Optimization uses exact range/image-CAD analysis as a terminal path on projected KKT and singular active-boundary loci and recovers attained endpoint witnesses; a complete symbolic image decomposition in `critical_values` is not implemented. |
| Corners/nonsmooth boundary strata | Partial implementation | `region_nonsmooth_locus` detects algebraic singularities plus transverse corners/ridges, and `region_active_boundary_strata` classifies active inequality sets. Full Whitney/regular stratification remains open. |
| Betti numbers/homology | Missing subsystem | Oriented CAD incidence is a medium-sized extension; full roadmaps are substantially harder. |
| Globally minimal formula representation | Fundamental/representation-dependent | Local deterministic simplification can improve, but no globally minimal representation is promised. |

## Parametric geometry and affine specializations

Bounded structural charts are currently provided for recognized points,
intervals, boxes, simplices, polygons, polyhedra, parallelograms,
parallelepipeds, and bounded `ParametricRegion` objects. A formula region needs
explicit finite clipping bounds to obtain an identity-chart cover without CAD.
If parameter-domain dimension or map rank is not certified uniformly, dimension
queries fall back to the complete CAD implementation rather than using a
generic rank as a universal statement.

`parametric_map_degree` is the generic complex algebraic degree of a rational
map. It does not count only real fibers or impose semialgebraic parameter-domain
restrictions; those are separate data needed by change-of-variables and
intrinsic-integration code.

Affine box clipping is intentionally restricted to one- and two-dimensional
affine subspaces. General polyhedral intersections continue to use the ordinary
region/CAD machinery.
