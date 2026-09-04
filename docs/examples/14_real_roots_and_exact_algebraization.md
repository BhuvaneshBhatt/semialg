# Real roots and exact transcendental algebraization

Contrast principal rational powers with explicit real roots, then compute exact ranges of commensurate trigonometric and exponential expressions.

## Executable example

```python
import sympy as sp

from semialg import equivalent, function_domain, function_range

x, t = sp.symbols("x t", real=True)

principal_domain = function_domain(x ** sp.Rational(1, 3), [x])
real_root_domain = function_domain(sp.real_root(x, 3), [x])
trig_range = function_range(sp.sin(x) + sp.cos(2 * x), variables=[x], value_symbol=t)
exp_range = function_range(sp.exp(x) + sp.exp(-x), variables=[x], value_symbol=t)

print("principal cube-root domain:", principal_domain)
print("explicit real cube-root domain:", real_root_domain)
print("trigonometric range:", trig_range)
print("exponential range:", exp_range)

assert principal_domain == (x >= 0)
assert real_root_domain is sp.true
assert equivalent(trig_range, (t >= -2) & (t <= sp.Rational(9, 8)), [t])
assert equivalent(exp_range, t >= 2, [t])
```

The matching executable file is `examples/gallery/14_real_roots_and_exact_algebraization.py`; the documentation test suite runs every gallery script.
