# Native CAD proofs

Use CAD when a polynomial question must hold on an entire real region, including exceptional parameter values. A root at a sample point is only part of the argument: projection must justify moving that root throughout its base cell, and lifting must retain every relevant section and sector. This guide explains how to construct and replay those obligations.

| Task | Interface | Guarantee and scope |
| --- | --- | --- |
| Ambient sign-invariant decomposition | `decomp_lazard_complete` or `cad(..., strategy="lazard-complete")` | Full Lazard projection and native valuation lifting over supported exact coefficient fields |
| Automatic nullification fallback | `decomp_collins_complete` | Starts with coefficient projection and isolated-point delineation; restarts with full Lazard when needed |
| Formula-oriented EC reduction | `decomp_form_reduced_safe` | Recursive degree/order proof on admissible EC sections, with ambient fallback |
| Independent native replay | `verify_cad` | Reconstructs projection and stack geometry; does not accept completion flags as evidence |
| Collision-fiber branch limits | `certify_branch_correspondence` | Degree-preserving limits from delineable source cells; ambiguous one-dimensional bases use exact Sturm germs |

[Executable application](https://github.com/BhuvaneshBhatt/semialg/blob/main/examples/native_cad_proofs.py) and [regression tests](https://github.com/BhuvaneshBhatt/semialg/blob/main/tests/test_cad_geometric_replay.py) exercise these interfaces against exact mathematical answers.

## Full projection and moving Taylor selectors

For an irreducible factor basis in the highest variable, the full Lazard operator projects contents, leading coefficients, trailing coefficients, discriminants and all pairwise resultants. It applies the same operator recursively. The algorithm does not discard a positive-dimensional stratum where a projection polynomial becomes identically zero.

At a parent sample, native lifting finds the lexicographically first nonzero Taylor coefficient by iterated Hasse derivatives. A positive valuation records both the orders and the mixed-derivative residue in the original coordinates. That polynomial remains a moving selector on a base cell with constant valuation. Substituting an algebraic sample into an unrelated global number field is unnecessary. Exact real algebraic ground coefficients are retained in their coefficient field separately from the selected parent tower.

<!-- examples: executable -->

```python
import sympy as sp
from semialg.cad_algorithms import decomp_lazard_complete, verify_cad

x, y, t, z, w = sp.symbols("x y t z w", real=True)
variables = (x, y, t, z, w)
polynomial = w**2 + y * z - x**2
source = decomp_lazard_complete((polynomial,), variables)
proof = verify_cad(source, polynomials=(polynomial,), return_result=True)
assert proof.valid, proof.failures
assert proof.checked_cells == sum(map(len, source.cells_by_level.values()))
assert proof.checked_valuation_claims > 0
```

The free `t` coordinate makes the nullifying base stratum positive-dimensional. The relevant intermediate Taylor selector still separates `z=0`. `source.backend` names the effective backend, including an automatic restart from ordinary coefficient projection. Cooperative limits on the direct decomposition functions can stop a computation; an exhausted budget does not establish infeasibility.

## Recursive equational constraints

A formula-wide equation can restrict lifting to its sections. Derived ECs must be logical consequences of the formula: the resultant of two implied equations is necessary, but an individual factor of an implied product need not vanish. The implication records retain parent equations and the eliminated variable.

The highest EC step needs sign invariance on its EC sections. Intermediate steps need order invariance there. The stronger intermediate operator retains every active leading coefficient and discriminant, every pairwise resultant, and the designated EC's coefficients. `ECLevelHypotheses` records these families and the degree guards. The driver checks guards only where the lower designated ECs hold; false lower EC cylinders need no higher cuts for this formula.

```python
from semialg.cad_algorithms.reduced import decomp_form_reduced_safe
from semialg.formula import parse_formula

x, y, z = sp.symbols("x y z", real=True)
matrix = parse_formula(sp.And(sp.Eq(y - x, 0), sp.Eq(z**2 - y, 0)))
reduced = decomp_form_reduced_safe(matrix, (x, y, z), backend="mccallum")
assert not reduced.used_fallback
assert reduced.certificate.checked_refinement_cells == 0
assert verify_cad(reduced.cad, formula=matrix)
```

The projection admits a real solution exactly when `x>=0`. A positive Taylor valuation does not prove that an EC-only lift preserves the original equation, and a reference CAD's sample agreement cannot establish missing delineability hypotheses. Failed degree/order induction therefore selects the ambient fallback. Sign-reduced and TTICAD drivers use separate reference-certification formats; those formats are outside `verify_cad`.

## Collisions and complex-to-real births

An endpoint's multiplicities do not identify which roots existed on the source cell. For `((y+1)**2+x)*(y-2)` with `x>0`, the sole real branch is `y=2`. The double root born at `x=0, y=-1` must not capture that branch's rank.

```python
from semialg.algebraic.samples import RationalSample
from semialg.topology.branch_correspondence import (
    certify_branch_correspondence,
    replay_branch_correspondence,
)

x, y = sp.symbols("x y", real=True)
limit = certify_branch_correspondence(
    ((y + 1) ** 2 + x) * (y - 2),
    y,
    0,
    (x,),
    (RationalSample(sp.S.One),),
    (sp.S.Zero,),
)
assert limit.value == 2
assert replay_branch_correspondence(limit)
```

The proof enumerates ordered collision allocations that conserve source multiplicities. Missing endpoint multiplicity must occur in conjugate pairs. A limit is certified when every allocation assigns the selected branch to the same endpoint. If allocations disagree and the base is one-dimensional, a segment from the endpoint to the source sample stays inside its base interval. Exact Sturm signs over the resulting algebraic rational-function field determine the eventual root ranks at rational barriers as the path parameter tends to zero from above. There is no sampled epsilon or numerical branch matching.

The source-cell assumptions matter: the selector must be delineable with constant multiplicities on that cell, and the target must lie in its closure. The incidence engine supplies this context. Standalone certificate replay recomputes the algebraic claims for the recorded statement; it does not independently prove arbitrary source-cell hypotheses. One-parameter degree loss uses reciprocal charts; exterior Sturm ranks select negative or positive infinity without attaching an escaping section to a finite endpoint. Multibase degree loss, ambiguous higher-dimensional base limits and path-dependent collisions still require further refinement. The solid-torus component example is covered by degree-preserving collision correspondence.

## What replay proves

`verify_cad` bypasses native memoization, reconstructs the required projection and rebuilds each stack over verified parent samples. It checks projection provenance, local delineating families, original and Taylor-residue selectors, absolute root ranks, isolating intervals, cached root expressions, neighboring bounds, fiber coverage and exact recorded signs. Recursive truth proofs additionally require the actual formula and replay their EC implications and degree/order hypotheses. `polynomials=` binds an ambient certificate to an external polynomial family.

The default result is a Boolean; `return_result=True` returns the failure reasons and counts of checked cells, selectors and valuation claims. These counters are derived outcomes. The verifier shares the exact arithmetic kernel with construction, and its supported formats exclude finite-variety and reduced/TTICAD certificate formats. A successful sample-sign check alone has a narrower guarantee.

Exact image certification and implication minimization also retain their proofs when reusing work. Diagonal-quadratic translations transform the objective and domain by the same real bijection before image caching. Related Boolean queries reuse truth bitsets only on a CAD containing the same polynomial boundaries within one computation context. A new boundary requires a new context.

See [algorithm references](references.md#cylindrical-algebraic-decomposition), [limitations](limitations.md) and [certificate-field audits](certificate_audits.md).
