# Region representations and when to use them

`semialg` exposes several region representations because they serve different performance and usability roles. They are interoperable, but they are not interchangeable implementation details.

| Representation | Best use | Builds CAD? | Structural identity |
|---|---|---:|---|
| Raw Boolean formula | One-off reasoning and QE | Only when an operation needs it | SymPy formula |
| `StandardRegion` | Named familiar shapes | No | Shape parameters |
| `SemialgebraicRegion` | Canonical symbolic region computation | Lazily | `(formula, variables)` |
| `CADRegion` | Repeated exact geometry/topology on one decomposition | Once, then reused | Wrapped `CADResult` |
| `StructuredCADCell` | Cell-level algebraic bounds and numerical meshing | Reuses source CAD | One cylindrical cell |

The normal flow is:

```text
StandardRegion or formula
          |
          v
SemialgebraicRegion
          |
          | ensure_cad() / as_cad_region()
          v
       CADRegion
        /     \
       v       v
CADCellComplex  StructuredCADCell
                    |
                    v
               numerical mesh
```

## Raw formulas

Use a raw SymPy Boolean formula when the computation is genuinely one-off:

```python
import sympy as sp
from semialg.reasoning import region_subset

x = sp.symbols("x", real=True)
a = sp.And(x >= 0, x <= 1)
b = x**2 <= 1
assert region_subset(a, b, (x,))
```

This keeps input lightweight. There is no persistent context or CAD cache attached to the formula.

## `StandardRegion`

Named regions are convenient for geometry and integration:

```python
from semialg import BoxRegion

box = BoxRegion(((0, 1), (-2, 2)))
```

Convert to the unified symbolic model when you want general region algebra:

```python
import sympy as sp
from semialg import as_semialgebraic_region

x, y = sp.symbols("x y", real=True)
region = as_semialgebraic_region(box, (x, y))
```

## `SemialgebraicRegion`

This is the canonical symbolic region object:

```python
from semialg import SemialgebraicRegion

region = SemialgebraicRegion(x**2 + y**2 <= 1, (x, y))
```

It carries the exact formula and ambient coordinate variables structurally, with lazy reusable caches for a `SemialgebraicContext`, quantifier-free lowering, and an optional CAD. Use it for Boolean operations, membership, projection/image/preimage, topology, measure, integration, and exact region relations.

## `CADRegion`

Build this when several queries will use the same decomposition:

```python
cad_region = region.as_cad_region()

cad_region.locate_point((0, 0))
cad_region.sign_vector((0, 0))
cad_region.cell_complex()
cad_region.euler_characteristic()
cad_region.measure()
```

These operations reuse the existing CAD. Boolean combination also reuses a shared decomposition when all operands actually refer to the same CAD.

## `StructuredCADCell`

Structured cells expose nested algebraic bounds, including CAD-certified algebraic root branches. They are the right level for boundary descriptors, iterated integrals, curve/surface evaluation, and simplicial meshing.

Do not convert numerical mesh vertices back into exact region facts. Numerical geometry is a presentation/approximation layer over the exact CAD, not a replacement for it.
