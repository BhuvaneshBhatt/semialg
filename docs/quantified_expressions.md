# Quantified expressions

`semialg` provides small symbolic `Exists` and `ForAll` Boolean nodes because
released SymPy versions do not expose general public quantifier expression
classes.

```python
import sympy as sp
from semialg import Exists, ForAll

x, y = sp.symbols("x y", real=True)
phi = Exists(x, sp.Eq(x**2, y))
psi = ForAll((x, y), x**2 + y**2 >= 0)
```

The nodes are representations, not quantifier-elimination algorithms. They
integrate with SymPy Boolean expressions and implement lexical binding:

- bound variables are removed from `.free_symbols`;
- `.subs(...)` does not replace a variable bound by the quantifier;
- substitutions that would capture a free symbol alpha-rename the conflicting
  binder first;
- vacuous binders and quantified Boolean constants simplify away.

The ambient quantifier domain comes from the surrounding semialg problem.
Narrower domains should be written explicitly in the body. For example, an
integer periodic index is represented as

```python
k = sp.Symbol("k", real=True)
periodic = Exists(
    k,
    sp.And(
        sp.Contains(k, sp.S.Integers),
        x - 2 * sp.pi * k > 0,
        x - 2 * sp.pi * k < sp.pi,
    ),
)
```

The transcendental periodic-reconstruction layer uses this representation for
periodic root and interval families instead of encoding quantification through
`Mod` or `ImageSet`.


## Using quantified expressions with semialg solvers

Expression-facing code should prefer `Exists` and `ForAll` over constructing
raw `("exists", x)` / `("forall", x)` tuples.  The tuple/block form remains an
internal normalized representation used by the QE engines.

```python
import sympy as sp
from semialg import Exists, ForAll, apply_quantifiers, split_quantifiers
from semialg.solve import reduce_complete_expr

x, y = sp.symbols("x y", real=True)

statement = ForAll(x, Exists(y, sp.Eq(x + y, 0)))
reduce_complete_expr(statement)
# True
```

`apply_quantifiers(matrix, prefix)` and `split_quantifiers(formula)` provide an
explicit bridge for code that must interoperate with older internal APIs.  A
`ParsedPrenexFormula` also exposes `.quantified_expr`, so parsing text and
reconstructing a first-class quantified expression round-trip cleanly.

The transcendental state builder accepts a leading semialg quantifier prefix
directly:

```python
from semialg.solve.transcendental import build_trans_state

state = build_trans_state(
    ForAll(x, Exists(y, sp.Eq(sp.sin(x), y))),
    free_variables=(),
)
```

Internally this is lowered to quantifier blocks for dispatch, but callers do not
need to construct those blocks themselves.
