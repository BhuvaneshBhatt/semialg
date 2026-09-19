# Tutorial: sparse SOS search and exact certification

For a polynomial \(p\), an SOS Gram certificate has the form \(p=z^TQz\) with \(Q\succeq0\). `semialg` separates finding `Q` from proving that it is valid.

## Newton-polytope basis reduction

```python
import sympy as sp
from semialg.sos_certificates import (
    plan_sos_search,
    search_sos_certificate,
    sparse_sos_monomial_basis,
    verify_psd_exact,
)

x, y = sp.symbols("x y")
p = x**8 + y**8

basis = sparse_sos_monomial_basis(p, (x, y))
plan = plan_sos_search(p, (x, y))
print(basis)
print(plan.gram_dimension, plan.dense_gram_dimension)
```

A Gram monomial \(x^\alpha\) can occur only if \(2\alpha\) lies in the Newton polytope of `p`. The sparse basis applies that condition exactly, so discarded monomials cannot participate in any SOS Gram representation.

## External SDP solvers are proposal engines

When an optional CVXPY-backed solver is installed, it can be requested explicitly:

```python
candidate = search_sos_certificate(p, (x, y), backend="clarabel")
```

The external solver's floating matrix is **not** a certificate. `semialg` rationally reconstructs free coordinates, solves the Gram coefficient equations exactly, and accepts a certificate only if the exact polynomial identity holds.

## Exact PSD verification

```python
Q = sp.Matrix([[1, 1, 0], [1, 1, 0], [0, 0, 0]])
check = verify_psd_exact(Q)
assert check.verified
print(check.pivots)
```

PSD is checked by exact LDL/congruence elimination, including singular zero-pivot cases. No eigenvalue tolerance is used.

This yields the common certificate pipeline:

```text
numerical or symbolic Gram proposal
  -> exact affine reconstruction
  -> exact p = z.T*Q*z identity
  -> exact Q >= 0 verification
  -> SOS certificate
```

Failure to find an SOS certificate does not prove that a polynomial is negative—or even that it is not SOS. For a certified nonnegativity decision, use the higher-level positivity portfolio, which can continue to exact critical-value/ARS and CAD backends.
