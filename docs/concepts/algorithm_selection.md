# Algorithm-selection map

`semialg` is exact-first: specialized algorithms are used when their hypotheses can be established, while complete exact fallbacks remain available where the package supports them. A fast path may decline; it may not silently weaken the mathematical claim.

## One map of the major operations

| Goal | Preferred specialized backends | Exactness / completeness boundary | Conservative fallback | Optional dependencies |
|---|---|---|---|---|
| Feasibility / witnesses | affine presolve, exact Fourier–Motzkin, incidence decomposition, RUR for finite equality systems, ARS polar reduction, validated exact samples | A witness proves SAT; specialized equality solvers report completeness explicitly. Failure to find a witness never proves UNSAT. | lazy/partial CAD with exhaustive traversal for UNSAT; complete CAD where required | `python-flint` can accelerate exact algebra but is not required for correctness |
| Quantifier elimination | safe substitution, quadratic virtual substitution, zero-dimensional algebraic elimination, quantifier-aware variable ordering | Virtual substitution applies only to its supported low-degree fragment. Reduced CAD is used only when its side conditions are certified. | complete Collins-style CAD | none required |
| Optimization / ranges | affine reductions, stationary/KKT systems, active-boundary recursion, exact algebraic candidate solving, parameter stratification | Candidate generation alone is not a global proof; certification compares all required exact candidates/branches. | complete semialgebraic decision/QE machinery for unresolved global comparisons | optional external numerical tools may propose candidates, never certify them |
| Algebraic decomposition | regular chains, modular Gröbner reconstruction, independent localization, zero-dimensional local decomposition, recursive GTZ contraction/saturation | A decomposition is complete only after exact ideal reconstruction, radical/primary checks, and replay data succeed. | broader recursive GTZ path; unsupported coefficient domains remain explicit | `python-flint` is an optional exact accelerator; Singular is useful only as an external differential oracle |
| Polynomial positivity | structural tests, coercive critical-value reasoning, sparse Newton-polytope SOS search | Numerical SDP output is never proof. SOS is accepted only after exact Gram identity reconstruction and exact PSD verification. SOS incompleteness is not negativity. | Zeng/ARS negativity search followed by complete CAD in the certified decision portfolio | optional `symbopt`; optional CVXPY-backed Clarabel/SCS/MOSEK proposal backends |
| Integration / measure | standard-region conversion, exact univariate bounds, CAD cells, algebraic root-function endpoints for supported parameter fibers | Integration is exact only when bounds/cells and antiderivative/evaluation steps are supported exactly. General singular or multidimensional parametric families may remain unsupported. | more general CAD-based region decomposition when available | none required |

## Shared planning rule

The planner first normalizes the problem and extracts structural information: variables, polynomial degrees, equality structure, sparsity, quantifier blocks, finite-dimensionality, and estimated algebraic/CAD/SOS cost. It then tries the cheapest backend whose hypotheses imply the requested contract.

```text
input
  -> exact normalization / structural analysis
  -> specialized exact backend if certified applicable
  -> exact reconstruction and verification
  -> conservative broader backend when the specialized method declines
  -> explicit unsupported/unknown when no certified route applies
```

This distinction matters most for incomplete search procedures. An SOS search that finds no certificate, a witness search that finds no point, or a specialized decomposition routine that declines is **not** allowed to manufacture a mathematical negative answer.

## CAD-specific planning

CAD cost is highly order-dependent. Automatic ordering respects logical structure: free variables may be reordered among themselves and variables inside one homogeneous quantifier block may be reordered, but variables do not cross an `exists`/`forall` boundary. Reduced/equational-constraint projection is used only when its correctness conditions are certified; otherwise the implementation falls back to complete projection.

The projection representation is factorized. A reducible source polynomial may be replaced internally by its exact squarefree factors, because sign-invariance of those factors plus the required resultants determines the source polynomial's sign. Internal projection storage is therefore not a promise to retain every input polynomial verbatim.

## Algebraic and SOS cost controls

Exact algebra uses modular reconstruction when coefficient growth is likely to justify the proposal/certification overhead; small systems stay on direct exact arithmetic. Function-field towers can be compressed through certified primitive elements when depth/degree predicts expression growth.

SOS planning uses the Newton polytope of the polynomial to remove Gram-basis monomials that cannot occur in any SOS representation. External SDP solvers, when installed, are candidate generators only. Their floating Gram matrices are projected onto the exact coefficient-affine space and then checked by exact LDL/congruence PSD verification.

See [Exactness and certification](exactness_and_certification.md), [Certificate model](certificate_model.md), and the [Performance guide](../guides/performance.md).
