# Toric and lattice algebra

The `semialg.algebraic` namespace contains exact integer-lattice operations used by algebraic statistics and other monomial models. These functions are generic algebra; they do not impose probability normalization or statistical semantics.

```python
import sympy as sp
from semialg.algebraic import integer_kernel, markov_basis, toric_ideal

A = (
    (1, 1, 1, 1),
    (1, 1, 0, 0),
    (1, 0, 1, 0),
)
p = sp.symbols("p0:4")

integer_kernel(A)
toric_ideal(A, p)
markov_basis(A)
```

`integer_kernel(A)` uses a Smith-normal-form decomposition over `ZZ`, so the returned vectors form a saturated integer lattice basis rather than independently scaled rational nullspace vectors.

`toric_ideal(A)` eliminates a Laurent monomial parameterization. Inverse parameter variables are included explicitly, which makes the algorithm valid for negative as well as nonnegative integer exponents.

`binomial_ideal(moves)` constructs the binomials associated with supplied lattice moves. `lattice_ideal(moves)` saturates those binomials by the product of the coordinate variables, so a lattice basis generates the full lattice ideal rather than only its unsaturated basis ideal.

`markov_basis(A)` extracts exponent-difference moves from a binomial Gröbner generating set of the toric ideal. Such a generating set connects every nonnegative integer fiber of `A`.
