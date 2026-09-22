# Optimization and range reference
## Family contract

**Mathematical return.** Optimization returns exact infimum/supremum information, attainment and witnesses where available; range APIs characterize all attainable scalar values.

**Exactness and certification.** Exact candidate generation is distinct from global certification. Structured results expose certification/attainment information where relevant.

**Algorithm.** The implementation combines structural preprocessing, exact stationary/KKT and boundary candidates, algebraic solving, exact global comparison, and parameter stratification. Supported affine parametric relations may be reconstructed directly before generic QE.

**Complexity and limitations.** Nonconvex and parameterized problems may require CAD/QE and can be expensive. Use `argmin_set`/`argmax_set` when the entire optimizer locus matters.




## Primary API overview

This table is the substantive coverage target for the primary APIs assigned to this reference page. Each entry states the API's primary role; the family contract and detailed sections below explain shared algorithms, exactness guarantees, and limitations. It is maintained together with `docs/reference/primary_api_manifest.toml`, and documentation tests require every root-level primary API to map here rather than merely appearing in the generated public index.

| API | Kind | Role / return |
|---|---|---|
| `function_range` | function | Return a quantifier-free formula describing a real function range. |
| `FunctionRangeResult` | class | Exact range summary for a supported semialgebraic image problem. |
| `OptimizationResult` | class | Exact optimum summary for supported semialgebraic problems. |
| `CertifiedDecisionAttempt` | class | One backend attempt in a certified polynomial-decision portfolio trace. |
| `CertifiedDecisionResult` | class | Final certified decision plus ordered backend attempts, witness, and certificate. |
| `SOSCertificate` | class | Exact Gram-matrix sum-of-squares certificate. |
| `SOSSearchPlan` | class | Sparse Newton-basis/Gram-size estimate controlling automatic SOS search. |
| `PSDVerification` | class | Exact LDL/congruence PSD verification result. |
| `SOSSearchResult` | class | Optional SOS-search result whose certificate is accepted only after exact verification. |
| `plan_sos_search` | function | Estimate sparse Newton Gram basis and decide whether auto mode should launch SOS search. |
| `sparse_sos_monomial_basis` | function | Compute the exact Newton-polytope-filtered SOS Gram basis. |
| `search_sos_certificate` | function | Ask optional `symbopt` for an SOS candidate and verify it exactly. |
| `verify_psd_exact` | function | Verify exact PSD using symmetric LDL/congruence elimination. |
| `verify_sos_certificate` | function | Verify an exact Gram identity and positive-semidefinite Gram matrix. |
| `PolynomialNegativityResult` | class | Certified negative-point/nonnegativity/incomplete result for a polynomial. |
| `find_negative_point` | function | Return an exact point where a polynomial is negative, or certified nonnegativity. |
| `polynomial_nonnegative` | function | Return a certified global polynomial nonnegativity Boolean, or the SOS → Zeng/ARS → CAD portfolio trace with `return_result=True`. |
| `zeng_negative_point` | function | Specialized exact critical-value decision backend for polynomial negativity. |
| `semialg.parameters.root_count_conditions` | function | Return parameter conditions grouped by distinct real-root count. |
| `semialgebraic_maximize` | function | Return an exact maximum/supremum for a polynomial semialgebraic problem. |
| `semialgebraic_minimize` | function | Return an exact minimum/infimum for a polynomial semialgebraic problem. |
| `solvability_conditions` | function | Return parameter conditions for real solvability of a constraint system. |


## `polynomial_nonnegative`

```text
polynomial_nonnegative(
    polynomial, variables, *, strategy="auto", sos_backend="auto",
    return_result=False, random_lines=8, seed=1234
)
```

This is the single public entry point for certified global polynomial
nonnegativity. By default it returns the mathematical Boolean. Set
`return_result=True` to receive `CertifiedDecisionResult`, which records the
selected backend, ordered attempts, witness, and certificate. `strategy` may be
`"auto"`, `"sos"`, `"zeng"`, `"ars"`, or `"cad"`. Automatic mode tries an
exactly verified SOS certificate, then the Zeng/ARS specialized route, then
complete CAD. `random_lines` and `seed` tune the Zeng witness-search stage,
including when that stage is reached through `strategy="auto"`.

## `semialgebraic_minimize`

```text
semialgebraic_minimize(
    objective, constraints=None, variables=None, *, domain="reals",
    return_result=False, certification="auto", range_cost_limit=2500,
    recursion_limit=4, parameters=None, return_stratified=False, eliminate_quantifiers=False
)
```

Computes an exact minimum/infimum for supported semialgebraic problems. The default return is ``[extremum, optimizer_points]``, where ``optimizer_points`` is a list of exact variable-to-value mappings. Set ``return_result=True`` for certification, attainment, method, and diagnostic metadata. The pipeline can use equality reduction, stationary/KKT systems, active-set pruning, finite RUR solving, positive-dimensional critical-locus recursion, exact candidate comparison, and CAD certification.

`OptimizationResult` distinguishes:

- `value`: exact optimum/infimum value when obtained;
- `points`: known attaining optimizer points;
- `attained`: whether the value is attained;
- `certified`: whether global optimality was established exactly;
- `method`: diagnostic method information.

## `semialgebraic_maximize`

The corresponding maximum/supremum API with the same direct ``[extremum, optimizer_points]`` default and certification model.

## Certification policy

- `"candidate"`: avoid the expensive full range fallback;
- `"auto"`: use a symbolic cost model to decide whether that fallback is reasonable;
- `"complete"`: permit complete exact range certification regardless of the automatic cost estimate.

## `polynomial_locus_dimension(equations, variables)`

Returns exact dimension information used to distinguish empty, zero-dimensional, and positive-dimensional polynomial critical loci.

## `function_range`

```text
function_range(
    expression, constraints=None, variables=None, *, value_symbol=None,
    domain="reals", method="qe", return_result=False,
    parameters=None, return_stratified=False, eliminate_quantifiers=False
)
```

Computes an exact real image/range condition by introducing a graph relation and eliminating source variables.  Polynomial and rational expressions use direct graph equations; `Abs`, `sign`, rational powers, `Min`, `Max`, and finite `Piecewise` expressions use the shared semialgebraic function-graph engine.

Before range projection, uniquely solved algebraic equalities may be back-substituted to reduce dependent variables.  The range engine also recognizes two exact finite algebraizations of transcendental-looking expressions:

- **commensurate trigonometric polynomials/rational expressions** in `sin(a*x)` and `cos(a*x)`, when all frequencies are rational multiples of one real fundamental frequency and no bare source-variable dependence remains.  They are mapped to polynomial/rational expressions on `s**2 + c**2 = 1`;
- **commensurate exponential/hyperbolic polynomials/rational expressions** in `exp(a*x)`, `sinh(a*x)`, and `cosh(a*x)`, when all rates are rational multiples of one real fundamental rate.  They are mapped through `t = exp(g*x)` with the exact side condition `t > 0`.

These transformations are accepted only when the transformed expression and constraints no longer depend on the original source variable.  Noncommensurate frequencies, mixed transcendental families, phase shifts not covered by the current rules, and expressions such as `sin(x) + x` are declined rather than relaxed to independent auxiliary variables.  This preserves exactness: CAD/QE always receives a genuinely equivalent semialgebraic problem.

Ordinary rational `Pow` retains SymPy principal-branch semantics. Use `sympy.real_root` when a real odd root is intended; integer powers of canonical odd `real_root` expressions are recognized directly and encoded with a single real-root graph relation.

## Parametric results

With `parameters=[...]` and `return_stratified=True`, optimization/range APIs can return guarded exact parameter-dependent relations rather than forcing a second QE solely to create a compact presentation.

See [Optimization](../optimization.md), [Function range](../function_range.md), and [Performance](../guides/performance.md).


When `return_stratified=True`, `eliminate_quantifiers=True` requests a second complete-CAD QE pass on each exact parametric optimum/range relation. The default is `False` to preserve the cheaper exact first-order relation.

## Critical values and metric optimization

`critical_values(expression, region, variables)` exposes exact values at isolated KKT/singular candidates, certifies constant objective values on positive-dimensional connected KKT components when possible, and includes attained global extrema. For metric geometry, `distance_to_region` and `distance_between_regions` minimize squared Euclidean distance through the same exact semialgebraic optimization pipeline; use `return_result=True` to retain closest-point assignments and the underlying `OptimizationResult`.

`bounding_box(region, variables)` is coordinate-wise exact optimization and returns `{variable: (infimum, supremum)}`; `return_result=True` additionally records endpoint attainment.
