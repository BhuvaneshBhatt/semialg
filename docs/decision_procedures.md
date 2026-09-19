# Decision procedures and inequality proving

The decision layer answers exact real-domain questions about formulas.

## Satisfiability and equivalence

```python
import sympy as sp
from semialg import is_satisfiable, is_tautology, implies, equivalent

x, y = sp.symbols("x y", real=True)

is_satisfiable(sp.And(x**2 + y**2 <= 1, x > 0, y > 0), [x, y])
# True

is_tautology(sp.Or(x < 0, x >= 0), [x])
# True

implies(x > 1, x**2 > 1, [x])
# True

equivalent(x**2 <= 1, sp.And(x >= -1, x <= 1), [x])
# True
```

## Inequality proving

The inequality provers reduce sign claims to satisfiability checks:

- `prove_nonnegative(f)` checks unsatisfiability of `f < 0`.
- `prove_positive(f)` checks unsatisfiability of `f <= 0`.
- `prove_nonpositive(f)` checks unsatisfiability of `f > 0`.
- `prove_negative(f)` checks unsatisfiability of `f >= 0`.

```python
from semialg import prove_positive, prove_nonnegative

prove_nonnegative((x - 1) ** 2, [x])
# True

prove_positive(x**2 + 1, [x])
# True

prove_nonnegative(x * y, [x, y], assumptions=sp.And(x >= 0, y >= 0))
# True
```

## Notes

These functions work over the real domain and are intended for polynomial and semialgebraic formulas. They are often used internally by simplification, region predicates, and optimization routines.

## Specialized algebraic feasibility and polynomial sign backends

For rational polynomial equality systems, `real_algebraic_feasibility` provides a
certified positive-dimensional alternative to CAD. It first obtains a replayable
certified equidimensional decomposition, constructs Aubry--Rouillier--Safey El Din
([ARS2002](references.md))
polar/critical systems over a rational function field in symbolic generic-point
parameters, proves the generic dimension drop, extracts an exact exceptional
parameter hypersurface from the generic Groebner basis, and constructs a rational
specialization outside that bad locus. The specialized system is independently
checked for the predicted dimension before recursion to the exact rational-univariate
solver. Elimination or decomposition failure is reported as incomplete; it is never
converted into an emptiness claim.
`solve_real_algebraic_set` is the convenience form that returns one exact witness,
`None` for certified emptiness, and raises when certification is incomplete.

The decomposition engine also uses triangular initial and separant splitters and
a squarefree regular-chain certification layer ([Kalkbrener1993](references.md), [ChenEtAl2013](references.md)).  A terminal triangular branch is
accepted as equidimensional only after exact saturation proves regular initials,
regular separants, and equality of the saturated chain ideal with the branch.
This allows ARS to consume prime/radical curve and higher-codimension branches
whose Groebner basis has more generators than the height. Each split is implemented by the exact identity
`V(I) = V(I + <h>) union V(I : h^infinity)` and is accepted only after radical
containment proves both branches strict. `initial_split` and `separant_split` are
recorded in the same replayable decomposition certificate as factor-incidence,
monomial, real-radical, and saturation strategies.

`zeng_negative_point` is the specialized exact polynomial-negativity backend. Fast
witness searches may prove that `f < 0` at an exact point. For the certified
coercive critical-value fragment of the Zeng-family strategy ([ZZ2004](references.md), [ZX2012](references.md)), the
backend proves attainment of the global minimum and solves zero-dimensional
gradient loci exactly. Positive-dimensional gradient loci are no longer an automatic
incomplete case: negativity on the critical locus is encoded by introducing a real
slack variable `u` and the equality `f*u**2 + 1 = 0`, then certified with the ARS
positive-dimensional equality-feasibility backend. Coercivity is certified separately from the literature-specific reduction; see the coercivity section below. Unsupported cases remain incomplete. The current backend is Zeng-family inspired rather than a verbatim implementation of [ZZ2004]. `polynomial_nonnegative` and `find_negative_point` provide
convenience APIs. `polynomial_nonnegative_decision` exposes the certified portfolio
trace. Automatic global nonnegativity dispatch tries an optional SOS search whose
Gram certificate must verify exactly, then Zeng (including its ARS handoff), and
finally complete semialgebraic decision/CAD. `prove_nonnegative` and
`prove_nonpositive` consume this same portfolio before general implication/QE.

For pure polynomial equality feasibility, `is_satisfiable`, `find_instance`, and
existential sentence reduction now try ARS after finite RUR solving and before
general CAD. An incomplete ARS result has no truth value and simply falls through.

`semialg` owns the exact SOS proof boundary. `SOSCertificate` stores an exact Gram
representation following [Parrilo2000](references.md) and [Parrilo2003](references.md),
and `verify_sos_certificate` checks both the polynomial identity and positive
semidefiniteness exactly. PSD verification uses exact symmetric LDL/congruence
elimination rather than enumerating all principal minors. Singular PSD matrices are
handled by the exact zero-pivot rule: a zero diagonal in a PSD matrix must have a
zero row and column.

SOS planning uses the Newton polytope. A monomial `x**alpha` is retained in the
Gram vector only when `2*alpha` lies exactly in `Newton(p)`. This is a necessary
condition for every SOS Gram representation, so the reduction is lossless and can
be much smaller than the dense `binomial(n + d/2, d/2)` basis. The planner records
both sparse and dense dimensions and applies its automatic budget to the sparse
problem.

Search can still be delegated to `symbopt`; the SOS layer also provides optional
CVXPY-compatible external SDP proposal backends (`cvxpy`, `clarabel`, `scs`, and
`mosek`). Numerical Gram matrices are not trusted: their free affine coordinates are
rationally reconstructed, the exact coefficient equations are solved again, and the
result is accepted only after exact Gram-identity and LDL PSD verification. Missing
solvers or unrecoverable numerical candidates remain incomplete results. Explicit
SOS requests bypass the automatic performance planner.
When current `symbopt` returns its richer exact `SOSCertificate`, `semialg` prefers
that recovered certificate over the flattened compatibility `gram_matrix` field. Rich-certificate
metadata is treated only as a conservative gate: `semialg` reconstructs the Gram proof
and independently rechecks the exact polynomial identity and PSD condition. Older
`symbopt` integrations using the flattened `monomial_basis`/`gram_matrix`
contract.

This keeps the proof boundary inside `semialg` while optimization-oriented search
remains in `symbopt`.

`find_negative_witness_fast` contains the one-sided accelerators. Odd-degree/ray
searches and deterministic pseudo-random rational lines are allowed to establish
satisfiability only after exact substitution verifies the returned witness. Failure
to find a witness has no logical meaning and never proves infeasibility. The
univariate step chooses rational points from certified isolating intervals; it
does not manufacture witnesses from floating-point root approximations. Search
counts and coefficient bounds are validated at the public boundary.

### Coercivity certification

Zeng coercivity is certified from exact positivity of the leading homogeneous form. The backend first uses the cheap even-monomial criterion, then optionally seeks an exact SOS radial margin `H - eps*(sum(x_i**2))**k`, and finally decides positivity on the real unit sphere by checking exact unsatisfiability of `sum(x_i**2) = 1` together with `H <= 0`. A cheap exact nonpositive-direction probe avoids unnecessary CAD work for obviously noncoercive or degenerate forms.

### Certified radical and minimal-prime decomposition

`semialg.algebraic_geometry.certified_radical_minimal_prime_decomposition`
constructs an exact reduced decomposition and separately certifies two claims.
`radical_complete=True` means the reconstructed intersection `J` of the returned
radical component ideals satisfies `I <= J <= sqrt(I)`, hence `J = sqrt(I)`.
`minimal_primes_complete=True` additionally means every component has an exact
primality certificate and the prime cover is irredundant, so those components
are precisely the minimal primes of `I`.  Squarefree regular chains are used as
radical/unmixed certificates but are never treated as prime merely from their
triangular shape.  Terminal regular chains are tested successively in the
quotient/fraction field of each preceding prefix.  `semialg` now represents
finite monogenic algebraic function-field towers exactly and factors leader
polynomials by recursive Trager-style norm/resultant descent.  Arbitrary nonzero
inputs are first decomposed exactly into monic squarefree factors with
multiplicities over the function-field tower; the Trager core then factors each
squarefree part. Shift selection is exact rather than budgeted: a shift is
exceptional precisely when the ``z``-discriminant of
``Norm(f(z - s*alpha))`` vanishes. For each squarefree part in characteristic
zero this is a nonzero polynomial in ``s``, so semialg tests
integer specializations exactly (via ``gcd(N_s, dN_s/dz)``) until it finds one
outside the finite bad set.  A certified
factorization recursively splits the branch; an irreducible tower proves
primality.  This is essential because ambient irreducibility can be
misleading: modulo `y - x**2`, the ambient-irreducible polynomial `z**2 - y`
becomes `(z-x)*(z+x)`.  Nonlinear algebraic towers are handled internally rather than delegated to
SymPy's coefficient-domain layer.  When norm factorization is inconclusive,
quotient-ring/saturation splitting remains a one-sided exact fallback and an
unresolved case stays explicitly incomplete.  The recursive regular-chain driver
also turns failed prefix regularity into exact zero-divisor splits and uses
regular-gcd/subresultant strata for competing leaders; every such branch is
validated by the closed/saturation identity before recursion continues.  This
follows the regular-chain decomposition tradition cited in `docs/references.md`.

### Radical and triangular-primality certificates

`radical_ideal` computes a radical from the unbounded-by-default recursive
squarefree regular-chain decomposition.  It returns generators only after the
component intersection has been checked both as containing the input ideal and
as lying in the input radical.  Unsupported branches remain incomplete rather
than being promoted to a radical result.

`certify_triangular_primality` records every successive leader extension,
including exact factor degrees and independent factorization reconstruction.
A regular chain is certified prime only when each leader polynomial is
irreducible over the fraction field of the certified quotient at the preceding leader;
a reducible stage instead returns exact ambient polynomial splitters.

### Replayable algebraic-decomposition certificates

Recursive regular-chain, radical, and minimal-prime results carry proof payloads
that are verified independently of the decomposition search.  The verifier
replays each recorded parent/child variety identity, recomputes closed and
saturated branches when a splitter is available, verifies squarefree regular
terminal chains, reconstructs the radical union by exact ideal intersection,
and independently checks Hilbert degree accounting.  Minimal-prime replay also
rechecks primality, irredundancy, and successive algebraic-function-field
factorization certificates.  Tampering with branch equations, splitters,
degrees, radical generators, or prime claims therefore invalidates replay.
