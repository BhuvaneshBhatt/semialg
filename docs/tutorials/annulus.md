# Tutorial: analyze an annulus end to end

```python
import sympy as sp
from semialg import SemialgebraicRegion

x, y = sp.symbols("x y", real=True)
r2 = x**2 + y**2
annulus = SemialgebraicRegion(sp.And(r2 >= 1, r2 <= 4), (x, y))
```

The annulus is bounded, closed, connected, two-dimensional, and has Euler characteristic zero:

```python
assert annulus.is_bounded()
assert annulus.is_closed()
assert annulus.dimension() == 2
assert len(annulus.components()) == 1
assert annulus.euler_characteristic() == 0
```

Its area is exact:

```python
assert sp.simplify(annulus.measure() - 3 * sp.pi) == 0
```

Build one CAD when you want several decomposition-level queries:

```python
cad_region = annulus.as_cad_region()
loc = cad_region.locate_point((sp.Rational(3, 2), 0))
assert loc.selected

complex_ = cad_region.cell_complex()
assert complex_.euler_characteristic() == 0
```

For visualization, triangulate the bounded full-dimensional CAD cells. The returned mesh is numerical but retains source-cell provenance.
