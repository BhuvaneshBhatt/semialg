# Real coverage of curve parametrizations

The expert namespace `semialg.curve_parametrization` checks the complete real image of birational normalization maps constructed by `algroots.curves`. Normalization genus and rationality stay in algroots; real coverage, domain restrictions and exceptional branch checks stay in semialg.

```python
import sympy as sp
from algroots.curves import PlaneCurve
from semialg.curve_parametrization import verify_real_curve_parametrization

x, y, t = sp.symbols("x y t", real=True)
p = PlaneCurve.from_polynomial(x**2 + y**2 - 1, (x, y)).rational_parametrization(t).parametrization
result = verify_real_curve_parametrization(p)
assert result.covers_curve and result.certified
assert result.infinity_point == (0, -1)
```

`verify_real_curve_parametrization(p, coordinates=None, domain=True)` first replays the exact birational map. On the open inverse chart every real curve point has a real parameter. Only the exceptional locus needs existential QE. Its finite parameter image is combined with the affine image of parameter infinity, when present. Exact satisfiability on the exceptional locus determines coverage; the generic birational identities prove the other implication. Algebraic constants require a real embedding.

`RealCurveParametrizationVerification` contains the finite graph, affine point at parameter infinity, projected real image, implicit curve, QE provenance and an `EquivalenceResult`. `covers_curve` is the exact coverage decision; `missing_condition` describes omitted real points, or is false when none exist. `coordinates` renames the two output variables and `domain` restricts the comparison to their declared semialgebraic region. A failed coverage check may supply a counterexample. Unsupported or resource-limited exact decisions raise; there is no sampling fallback.

## Isolated real singularities

```python
import sympy as sp
from algroots.curves import PlaneCurve
from semialg.curve_parametrization import verify_real_curve_parametrization

x, y, t = sp.symbols("x y t", real=True)
p = (
    PlaneCurve.from_polynomial(y**2 - x**2 * (x - 1), (x, y))
    .rational_parametrization(t)
    .parametrization
)
result = verify_real_curve_parametrization(p)
assert result.certified and not result.covers_curve
assert result.missing_condition.subs({x: 0, y: 0}) is sp.true
restricted = verify_real_curve_parametrization(p, domain=sp.Ne(x, 0))
assert restricted.covers_curve
```

The isolated real node `(0,0)` has only nonreal preimages on the normalization. Thus rationality over the reals does not imply coverage of every real point by real parameters. The cusp and a node with two real branches do pass complete coverage. These checks concern the real algebraic set; they do not automatically choose a graph branch or make a multivalued relation into a function. Use [functional-decomposition graph checks](function_verification.md) for those separate totality, uniqueness and composition obligations.
