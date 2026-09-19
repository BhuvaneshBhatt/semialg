# Exactness and certification

semialg follows an **exact-first** policy. This does not mean every symbolic-looking object is automatically a proof. The package distinguishes four ideas that users should keep separate.

## Public exactness contract

For exact/certified APIs, `semialg` does **not** silently replace an unsupported exact computation with a fixed-precision numerical approximation. A specialized exact backend may decline; the planner may then try another exact backend or a complete CAD path. If the requested exact operation is unsupported or cannot be certified by the implemented method, the operation should fail explicitly or return a result whose incomplete/candidate status is explicit in that API's contract.

**Exact does not mean inexpensive.** Complete real quantifier elimination can be extremely costly even when the final answer is simple.

The package uses several exact computational representations: ordinary SymPy rational/algebraic expressions, exact algebraic-root objects and isolating data, rational-univariate representations for finite systems, and semialgebraic Boolean formulas for sets and relations. Quadratic virtual substitution is also an exact low-degree QE/witness backend; it is not a numerical approximation.


## Exactness contract by operation

Numerical computation is permitted in different roles depending on the API. The decisive rule is whether a numerical value is allowed to establish the mathematical claim.

| Operation family | Numerical candidates allowed? | Exact validation required for an exact result? |
|---|---:|---:|
| CAD projection/lifting sign decisions | No fixed-precision proof substitute | Yes |
| Exact `find_instance` / exact sampling | Yes, as candidate generation only | Yes |
| Algebraic bound ordering | A numerical proposal may seed a rational candidate | Yes, by exact comparison |
| RUR / zero-dimensional solving | Heuristics may choose representation/order | Yes |
| Exact optimization | Candidate generation may use structural heuristics | Yes, for feasibility/value/global certification |
| Exact range / QE reconstruction | Heuristics may choose a backend or presentation | Yes |
| Plotting / discretization | Yes | No; these are presentation APIs |
| `exact=False` sampling/integration modes | Yes | No; inexactness is explicit |

A useful implementation principle follows: **numerics may suggest; they may not silently certify an exact API.**

## Exact representation

An object is exact when it denotes a mathematical value or set without finite-precision approximation. Examples include rational numbers, `sqrt(2)`, exact `RootOf`/algebraic-root objects, polynomial relations, and certified isolating intervals.

Exact representation answers **what object is being represented**. It does not by itself establish a global claim such as optimality.

## Certified conclusion

A conclusion is certified when semialg has established the relevant mathematical property by an exact procedure. Examples include:

- a CAD/QE proof that a formula is true or false;
- a validated exact witness for satisfiability;
- a proof that no feasible point has a better objective value;
- a verified parameter-stratification partition;
- a certified ordered algebraic root.

Structured result objects expose certification information where the distinction matters. For example, `OptimizationResult.certified` records whether global optimality was established exactly.

## Candidate or heuristic information

Some algorithms first generate candidates using KKT systems, active sets, structural heuristics, variable-order estimates, or pilot lifting. A candidate may be exact as a point while the statement “this is the global optimum” is not yet certified.

semialg should not silently upgrade candidate generation into a proof. Certification is a separate operation.

## Numerical approximation

Numerical methods are useful for plotting, exploratory sampling, diagnostics, and explicitly inexact workflows. They are not used as hidden proof substitutes in certified CAD/root-ordering paths.

For example, public sampling can opt into numerical random sampling with `exact=False`. That is different from an exact representative sample.

## Conservative failure

If semialg cannot establish an exact comparison, root order, CAD invariant, or global certificate within the supported method, the preferred behavior is to decline that exact step rather than make a fixed-precision guess.

This matters for very close algebraic values. A comparison based on 50 or 100 decimal digits can still be wrong; certified root ordering therefore uses exact algebraic comparison and isolating information.

## Exact optimum versus attained optimum

Consider

$$
\inf_{-1 < x < 1} x = -1.
$$

The value $-1$ is exact, but it is not attained. `OptimizationResult` keeps these facts separate through its `value` and `attained` fields.

Likewise, an exact candidate value can exist without a completed global certificate. Check `certified` when the distinction matters.

## Parameter-dependent exactness

For parametric optimization or range problems, a mathematically exact result may naturally be a first-order relation with explicit quantifiers rather than a compact quantifier-free `Piecewise` expression. semialg can retain that exact relation instead of automatically triggering a second expensive QE merely for presentation.

## Practical rule

When consuming a structured result, treat these questions independently:

1. **Representation:** Is the value/formula exact?
2. **Validation:** Has the returned witness/cell/branch been checked?
3. **Certification:** Has the requested global statement been proved?
4. **Approximation:** Was an inexact mode explicitly requested?

For what happens when one of these stages cannot be completed, see [Errors and failure modes](../guides/errors_and_failure_modes.md).
## Solver scope and exact transformations

Exactness does not imply that every symbolic function belongs to one solver. Polynomial
root certification, equivalence-preserving algebraization, and the existing transcendental
solver have distinct contracts. See
[Transcendental and algebraic solving scope](../guides/transcendental_scope.md) for the
canonical boundary.

