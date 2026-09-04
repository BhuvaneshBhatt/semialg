# Region-computation architecture

The region stack is layered so symbolic semantics do not depend on numerical meshing and expensive CAD work can be reused.

```text
input formula / StandardRegion
            |
            v
 normalization + variable policy
            |
            v
   SemialgebraicRegion
      |            |
      |            +--> direct symbolic operations
      |                 membership / Boolean algebra
      |                 projection / image / preimage
      |                 relation lowering / simplification
      |
      v
 SemialgebraicContext
      |
      v
   complete CAD  <---- specialized exact backends may avoid this
      |
      v
     CADRegion
   /     |       \
  v      v        v
point   CADCell   exact cylindrical
locate  Complex   integration
         |
         +--> Euler characteristic
         +--> local dimension
         +--> connectivity/topology
         |
         v
 StructuredCADCell
         |
         +--> algebraic boundary descriptors
         +--> delineable curve/surface evaluation
         +--> numerical simplicial mesh
```

## Structural versus cached state

A `SemialgebraicRegion` is structurally identified by `(formula, variables)`. Its context, quantifier-free lowering, and CAD are lazy caches and do not affect symbolic equality/hash.

## Why point location is a CAD operation

Repeated membership at exact resolved points should not run QE repeatedly. `CADRegion.locate_point()` descends the already-built CAD tree and then evaluates the projection-tower polynomials exactly at the point.

## Why cell incidence uses the source CAD

Closed cylindrical formulas can contain reconstructed radicals. Algebraically equivalent expressions such as `Abs(x)` are not polynomial formulas and should not be fed back into CAD/QE. `CADCellComplex` therefore prefers the source CAD's recursive exact `is_cell_in_closure` relation for incidence. Formula-based subset checking remains a fallback only when source decomposition metadata is unavailable.

## Why meshes remain downstream

Triangulation samples exact algebraic boundary functions numerically. Meshes preserve CAD cell provenance and can verify sampled conformity, but they are not used to establish exact set membership, topology, or singularity facts.
