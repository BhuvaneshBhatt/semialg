# Tutorial: project a 3-D semialgebraic set into 2-D

Consider the unit ball:

```python
import sympy as sp
from semialg import SemialgebraicRegion

x, y, z = sp.symbols("x y z", real=True)
ball = SemialgebraicRegion(x**2 + y**2 + z**2 <= 1, (x, y, z))
```

Project away `z`:

```python
shadow = ball.project((z,))
assert shadow.variables == (x, y)
```

computer algebra systemlly this is existential elimination:

```text
(x, y) is in shadow  iff  exists z: x^2 + y^2 + z^2 <= 1.
```

The result is the exact unit disk, not sampled projection geometry:

```python
disk = SemialgebraicRegion(x**2 + y**2 <= 1, (x, y))
assert shadow.equals_region(disk)
```

Projection can be expensive because it is QE. If several later operations need a CAD of the projected region, build it after projection:

```python
shadow_cad = shadow.as_cad_region()
shadow_cad.locate_point((0, 0))
shadow_cad.euler_characteristic()
```
