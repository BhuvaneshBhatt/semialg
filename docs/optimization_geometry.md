# Geometry for optimization and parameter families

Use `bounding_cover` to obtain a cheap exact outer coordinate box, `optimization_cover` to retain a full structural parameterization, and `region_family` to study exceptional parameter fibers. These interfaces supply geometric facts to an optimizer without introducing an optimization dependency into semialg. Symbopt owns incumbent management, reduced optimization, optimizer lifting and extremum relations.

## Structural bounds and covers

Unlike `bounding_box`, `bounding_cover` does not solve coordinate optimization problems. Its box may be loose. The structural optimization facts and exact family analysis require exact coefficients; floating data are rejected by these interfaces. Original region membership remains necessary. Supported structural cases include boxes, balls, ellipsoids, vertex sets, simplices, unions, transformed regions and supported polynomial maps.

```python
import sympy as sp
from semialg import (
    Ball,
    Box,
    ParametricRegion,
    TransformedRegion,
    bounding_cover,
    optimization_cover,
    region_family,
)

x, y, t, a = sp.symbols("x y t a", real=True)
assert bounding_cover(Ball((1, 2), 3), (x, y)) == {x: (-2, 4), y: (-1, 5)}
image = TransformedRegion(Box(((0, 1),)), (t, t * t), (t,))
assert bounding_cover(image, (x, y)) == {x: (0, 1), y: (0, 1)}
cover = optimization_cover(image, (x, y))
chart = cover.charts[0]
u = chart.parameters[0]
assert tuple(expr.subs(u, sp.Rational(1, 2)) for expr in chart.mapping) == (
    sp.Rational(1, 2),
    sp.Rational(1, 4),
)
```

An optimization cover includes the entire image, including singular points, boundaries and zero-dimensional strata. Union charts may overlap. Transformed covers compose the base maps, preserving parameter conditions. They are not integration charts with zero-measure pieces removed. No Jacobian or multiplicity weight is required when pulling back an objective through a covering map. A consumer must still enforce the real definedness of every coordinate map.

## Changing fibers

`RegionFamily` records the original formula, ordered fiber variables, coefficient parameters, strata and empty-fiber condition. Each `RegionFamilyStratum` carries its parameter guard, an exact sample, the fiber dimension and all selected cylindrical fiber cells. The result uses parameter coordinates first in a complete CAD.

```python
family = region_family(sp.Eq(x * x, a), (x,), (a,))
for parameter, expected in [(-1, (-1, 0)), (0, (0, 1)), (1, (0, 2))]:
    strata = [s for s in family.strata if sp.simplify(s.condition.subs(a, parameter)) is sp.true]
    assert len(strata) == 1
    assert (strata[0].fiber_dimension, strata[0].component_count) == expected
assert sp.simplify(family.empty_condition.subs(a, -1)) is sp.true
assert family.hardt.certified
```

For scalar CAD fibers, the result includes Hardt graph/band strata and exact component counts. Thus two points merging into one and disappearing into an empty fiber are separate cases. With several fiber coordinates, the result supplies cylindrical cells and dimensions; `component_count` is `None`. A CAD cell count is not a connected-component count, a Betti number or a general homeomorphism classification. Polynomial ball and sphere families have a separate radial analysis: negative radii give empty fibers, zero gives a point, and positive radii give a ball or sphere with known component counts and `betti_numbers`. A structural `fiber_cover` is provided in place of cylindrical cells for these families.

```python
rank_family = region_family(sp.And(sp.Eq(a * x, 0), sp.Eq(y, 0)), (x, y), (a,))
for parameter, dimension in [(0, 1), (1, 0)]:
    strata = [
        s for s in rank_family.strata if sp.simplify(s.condition.subs(a, parameter)) is sp.true
    ]
    assert len(strata) == 1 and strata[0].fiber_dimension == dimension
```

## Limits and ownership

All coefficient parameters must be declared and disjoint from fiber coordinates. Sign assumptions are retained; the topology interface requires continuous real coordinates. Unsupported structural bounds raise `NotImplementedError`; they do not trigger a coordinate optimizer. Family CAD and membership QE can be expensive. `region_family` accepts `max_cells` and `timeout` controls for its complete Collins CAD and requires completed decomposition. Cell limits count every lifting level, including cached decompositions. Deadlines are cooperative between exact operations; a single symbolic primitive can overrun the requested time. No unfinished CAD is cached as complete. Unsupported or incomplete computations raise an exception; a partial partition is never reported as complete. General higher-dimensional Hardt trivializations and topological invariants require further algorithms.

For exact optimizer loci and parameter-dependent infimum, attainment and unboundedness relations, use [symbopt's region interface](https://github.com/BhuvaneshBhatt/symbopt/blob/main/docs/guide/regions.md). Semialg continues to own the real geometry, CAD and quantifier elimination.
