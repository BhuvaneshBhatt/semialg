# Tutorial: certified primary decomposition

Primary decomposition writes an ideal as an intersection of primary ideals and exposes its associated primes. `semialg` returns exact components together with replayable certification when the supported decomposition completes.

Consider the embedded-prime ideal

\[
I=\langle (x+y)^2,\;y(x+y)\rangle.
\]

```python
import sympy as sp
import semialg
from semialg.algebraic_decomposition import associated_primes, primary_decomposition

x, y = sp.symbols("x y")
I = ((x + y) ** 2, y * (x + y))

result = primary_decomposition(I, (x, y))
assert result.complete
assert result.irredundant

for component in result.components:
    print("Q =", component.equations)
    print("sqrt(Q) =", component.radical)
    print("dimension =", component.dimension)

assert semialg.replay_certificate(result).verified
```

computer algebra systemlly this has the decomposition

\[
I=\langle x+y\rangle\cap\langle x^2,y\rangle,
\]

so the associated primes are the minimal prime \(\langle x+y\rangle\) and the embedded maximal prime \(\langle x,y\rangle\).

```python
primes = associated_primes(I, (x, y))
assert primes.complete
for prime in primes.primes:
    print(prime)
```

The important contract is not the particular internal route. Principal, monomial, zero-dimensional, and recursive localization/GTZ paths may all be selected depending on structure. A result is marked complete only after the exact component intersection reconstructs the source ideal and the required radical/primary facts have been certified.

For development or diagnostics, Singular can be used as an independent differential oracle, but it is not part of the proof: `replay_certificate` checks `semialg`'s own exact payload.
