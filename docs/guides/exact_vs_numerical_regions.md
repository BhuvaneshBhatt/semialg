# Exact versus numerical region results

The region subsystem intentionally mixes exact symbolic computation with numerical geometry. The boundary between those layers is explicit.

## Exact layers

The following are exact when they return successfully:

- `SemialgebraicRegion` formulas and Boolean operations.
- Membership and relation formulas.
- Projection, image, and preimage through exact substitution/QE.
- `SemialgebraicContext` normalization and reusable problem data.
- Complete CAD cells, exact samples, and algebraic root identities.
- `CADRegion.locate_point()` and tower sign vectors for exact resolved points.
- `CADCellComplex` cell dimensions and closure incidence.
- Connected components, local dimension, and Euler characteristic from the CAD complex.
- Exact region measure/integration when symbolic evaluation succeeds.
- Singular-locus conditions under the documented algebraic-boundary definition.

Exact does not mean inexpensive. Projection and first CAD construction may be doubly exponential in the worst case.

## Numerical layers

The following are numerical geometry utilities:

- `triangulate_cad_cell()` and `triangulate_cad_region()`.
- Numerical delineable curve evaluation.
- Numerical delineable surface evaluation.
- Floating-point mesh merging and tolerance-based vertex canonicalization.

They retain provenance back to exact CAD cells and algebraic root descriptors, but their floating-point coordinates are not exact certificates.

> **Meshing guarantee.** `triangulate_cad_region(..., require_conforming=True)` verifies that adjacent sampled CAD cells induce the same sampled simplex decomposition on their common faces. It does not prove ambient isotopy between the floating-point mesh and the exact semialgebraic set.

## Singular-locus convention

`SemialgebraicRegion.singular_locus()` computes an algebraic-boundary singularity notion based on defining boundary polynomials and gradient/Jacobian degeneracy.

> A corner formed by two individually smooth boundary hypersurfaces is not automatically an algebraic singularity. Use local dimension and the CAD cell complex when the question is about stratified or manifold structure rather than algebraic hypersurface singularity.

## Symbolic points are not Boolean membership queries

```python
region.contains((a,))
```

raises if membership still depends on unresolved symbols. Use:

```python
region.membership_formula((a,))
```

or `RegionElement` when you want a symbolic condition.

## Integration results

Region integration has three useful states:

1. an exact evaluated SymPy expression;
2. an exact unevaluated iterated integral over certified CAD bounds;
3. an explicitly numerical approximation requested by the caller.

The API does not silently replace an exact request with an unlabelled floating-point estimate.
