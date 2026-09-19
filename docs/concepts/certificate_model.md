# The certificate model

Many `semialg` algorithms look different computationally but share one architectural rule: **search is separated from proof**.

```text
candidate generation
    -> exact reconstruction
    -> exact mathematical verification
    -> replayable proof payload
```

A candidate may come from modular arithmetic, a numerical solver, a structural search, or an algebraic decomposition heuristic. None of those origins is trusted merely because the search succeeded.

## 1. Candidate generation

The first stage is allowed to be opportunistic. Examples include finite-field Gröbner images, a numerical SDP Gram matrix, a chosen primitive element, a regular-chain split, an equational constraint for reduced CAD, or stationary/KKT candidates for optimization.

This stage may be incomplete. Its failure normally means only "try another route."

## 2. Exact reconstruction

Candidate data is converted back into exact mathematical objects. Modular Gröbner computations use CRT and rational reconstruction. Numerical SOS candidates are projected onto the exact Gram coefficient equations. Algebraic candidates are represented by exact polynomials, ideals, algebraic numbers, or rational-function field elements.

Reconstruction is distinct from verification: producing an exact-looking object does not prove that it is correct.

## 3. Exact verification

The reconstructed object is checked against the defining mathematical identities. Typical checks include ideal containment/equality, exact polynomial reduction, exact reconstruction of an intersection, exact Gram identity `p = z.T*Q*z`, exact LDL/congruence PSD verification, sign-invariance/CAD invariants, and exact objective comparison.

Floating tolerances do not cross this boundary.

## 4. Replay

Where a result carries enough proof payload, `replay_certificate(result)` rechecks it through a public independent path. Replay does not trust stored `complete`, `verified`, or `certified` flags. A replay result with `verified=False` means the payload failed rechecking; `verified=None` means that result family does not expose enough replay information.

## How this appears in each subsystem

| Subsystem | Candidate/search | Exact reconstruction | Verification / replay |
|---|---|---|---|
| Modular Gröbner | finite-field bases and lucky-prime/signature filtering | CRT + rational reconstruction over `QQ` or `QQ(U)` | exact reduced Gröbner checks, source reduction, exact membership representations |
| Decomposition | regular-chain splits, localization, primitive separators, recursive GTZ branches | exact ideals over rational/function fields and contractions | exact branch identities, saturation, primary/radical checks, final intersection reconstruction |
| SOS | internal or external SDP Gram proposal | rational reconstruction plus exact solution of Gram coefficient equations | exact polynomial identity and exact LDL/congruence PSD verification |
| CAD / QE | projection choice, variable order, reduced/partial traversal | exact projection polynomials, algebraic sample points, cells | sign/truth invariance, reduced-CAD side conditions, complete fallback; CAD/QE replay validates recorded invariants |
| Optimization | KKT/stationary/boundary candidates or specialized reductions | exact algebraic candidate points/values and parameter branches | feasibility, objective identity, global exact comparison; replay recomputes the certified optimization problem when payload permits |

## Why external software does not become a proof dependency

Optional tools can make candidate generation dramatically faster. They do not change the proof boundary. For example, Clarabel, SCS, or MOSEK may suggest an SOS Gram matrix, and Singular may be used in development as a differential oracle, but `semialg` accepts its own mathematical result only after its exact verification path succeeds.

This makes optional acceleration compatible with reproducibility: removing an optional proposal backend can change performance or which fast path succeeds first, but it must not turn an uncertified numerical answer into an exact result.
