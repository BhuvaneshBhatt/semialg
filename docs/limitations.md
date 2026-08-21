# Limitations

`semialg` is conservative by design. It returns exact results for supported cases and should decline rather than turn an uncertified numerical approximation into an exact symbolic claim.
 Core CAD root isolation, root ordering, isolating intervals, sector sampling, topology incidence, and exact interval simplification follow this rule as well: fixed-precision `nroots`/`evalf` values are not accepted as exact certificates.

## Fundamental limits

- CAD/QE has doubly exponential worst-case behavior and can become expensive with only a few additional variables.
- Variable ordering strongly affects projection and lifting cost.
- General transcendental quantifier elimination lies outside real-closed-field CAD.
- Symbolic evaluation of an exact algebraic integral may fail even after the geometric decomposition is complete.

## CAD and quantifier elimination

The planner uses Brown/SOTD/NDRR-style heuristics, actual projection complexity on small shortlists, estimated lifting/root counts, coefficient-height and algebraic-degree estimates, and bounded pilot lifting for the leading candidate orders. Projection towers and exact algebraic subcomputations are cached through a solve-scoped `ExactComputationContext` backed by bounded process-wide LRUs. Lazy CAD exploits prefix truth, necessary equational constraints, derived resultant constraints, and existential section lifting where logically safe. These reduce practical cost but do not remove CAD's fundamental complexity.

Multi-level equational-constraint propagation is deliberately conservative: only logically necessary derived constraints are used for pruning, and universal variables are not restricted to EC sections.
 Pilot lifting deliberately prefers rational sample paths when available; exact algebraic-coefficient fibers are still supported in the main CAD path for low degrees, but the pilot remains a bounded cost probe rather than a complete test of every algebraic stack.

## Parameterized results

First-class conditional/stratified results are available for parameter-dependent answers, solvability conditions, root counts, parameter CAD strata, exact polynomial optimization, and function ranges. Parametric optimization and range branches carry exact first-order relations with explicit quantifier prefixes rather than automatically forcing an expensive second QE pass. Automatic quantifier-free/Piecewise elimination of those relations and fully parametric region integration remain incomplete.

## Algebraic root functions

`AlgebraicRootFunction` supports exact evaluation, specialization, certified same-stack comparisons, implicit derivatives, and cell-wide regularity checks. Comparisons outside a shared certified CAD stack may deliberately remain unknown. Radical expressions are presentation-only unless their branch identity is certified.

## Intrinsic integration

Regular triangular CAD graph strata can be integrated with the induced Hausdorff metric `sqrt(det(J.T*J))`. Singular strata are exposed explicitly by `stratify_intrinsic_solution`; uncertified algebraic sections are never treated as regular. General singular/non-graph semialgebraic stratification, multiple-chart geometry, and arbitrary singular manifolds remain incomplete.

## Polynomial optimization

The exact optimization pipeline includes equality-based dimension reduction, pruned active-set/KKT systems, singular and positive-dimensional active loci, RUR-backed zero-dimensional solving, recursive locus optimization, exact candidate comparison, and CAD global certificates. The hardest remaining cases are high-dimensional open/unbounded problems with no useful finite candidate, large Boolean expansions, and problems outside the polynomial/semialgebraic setting.

## Region integration

Arbitrary-dimensional typed CAD bounds and full-dimensional iterated-integral adapters are available. Exact symbolic integration can nevertheless remain unevaluated, and general non-disjoint unions or singular intrinsic strata may require decomposition that is not yet implemented.

## Engineering limits

Broad fallback handlers remain in some older subsystems and are being narrowed incrementally. A shared symbol resolver now covers the public subsystems found to be vulnerable to assumption-distinct duplicate symbols; lower-level helpers should continue migrating to it when they begin accepting user-facing string names.
