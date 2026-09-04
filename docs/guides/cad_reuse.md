# Reuse one CAD for many region queries

CAD construction is often the expensive part of exact semialgebraic computation. `SemialgebraicRegion` and `CADRegion` are designed so that a decomposition can be built once and reused.

```python
import sympy as sp
from semialg import SemialgebraicRegion

x, y = sp.symbols("x y", real=True)
region = SemialgebraicRegion(x**2 + y**2 <= 4, (x, y))

result = region.ensure_cad()
cad_region = region.as_cad_region()
assert cad_region.result is result
```

## Point location and sign vectors

```python
loc = cad_region.locate_point((1, 0))

loc.selected
loc.dimension
loc.cell_index
loc.sign_vector
```

The point is located in the existing CAD tree. The sign vector is evaluated exactly for the projection-tower polynomials at that resolved point.

Repeated queries do not rebuild the decomposition:

```python
cad_region.locate_point((0, 0))
cad_region.locate_point((1, 1))
assert region.ensure_cad() is result
```

## Topology from the same CAD

```python
complex_ = cad_region.cell_complex()
complex_.f_vector
complex_.euler_characteristic()
```

`CADCellComplex` groups selected cells by dimension and records exact codimension-one closure incidence. For curved algebraic cells, incidence uses the source CAD's recursive closure logic rather than reconstructing radicals and launching a second QE problem.

## Integration from the same CAD

```python
cad_region.measure()
cad_region.integrate(x + y, evaluate=False)
```

The latter returns exact iterated-integral pieces over certified cylindrical bounds.

## Boolean reuse

When operands share the same CAD, `cad_combine()` combines selected leaf sets directly:

```python
from semialg.cad_region import cad_combine

combined = cad_combine(left, right, op="union")
```

`combined.result.diagnostics["cad_reused"]` records whether that zero-rebuild path was used. If operands do not share a decomposition, a common refinement is constructed instead.

## When not to force CAD

Do not call `ensure_cad()` merely because the object can cache one. Projection, low-degree reasoning, linear optimization, and other specialized exact backends may be cheaper without a complete CAD. Build a `CADRegion` when repeated cell/topology/location/integration operations justify the decomposition cost.
