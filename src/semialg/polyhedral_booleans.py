"""Structure-preserving exact Boolean operations for convex polyhedra."""

from __future__ import annotations


def polyhedral_intersection(left, right):
    """Return the exact canonical intersection of two full-dimensional convex polytopes.

    Returns ``None`` for an empty intersection and raises ``NotImplementedError``
    when the exact structural backend does not cover the input dimensions.
    """
    from .standard_regions import Polytope

    if not isinstance(left, Polytope) or not isinstance(right, Polytope):
        raise TypeError("polyhedral_intersection requires two Polytope objects")
    if left.ambient_dimension() != right.ambient_dimension():
        raise ValueError("ambient dimensions do not match")
    n = left.ambient_dimension()
    if left.dimension() != n or right.dimension() != n:
        raise NotImplementedError(
            "structural polytope intersection currently requires full dimension"
        )
    A = left.h_representation()
    B = right.h_representation()
    matrix = [list(A.matrix.row(i)) for i in range(A.matrix.rows)] + [
        list(B.matrix.row(i)) for i in range(B.matrix.rows)
    ]
    offsets = [A.offsets[i, 0] for i in range(A.matrix.rows)] + [
        B.offsets[i, 0] for i in range(B.matrix.rows)
    ]
    from .polyhedral import HRepresentation

    rep = HRepresentation(matrix, offsets)
    vertices = rep.vertices()
    if not vertices:
        return None
    hull_dim = Polytope(vertices).dimension()
    if hull_dim != n:
        # Lower-dimensional contact is still exact, but use the general hull representation.
        from .convex_hull import convex_hull

        return convex_hull(vertices)
    from .convex_hull import convex_hull

    return convex_hull(vertices)


def polyhedral_boolean(operation, left, right):
    """Exact structure-preserving convex-polyhedral Boolean fast path.

    Intersection is constructed geometrically. Union/difference preserve a
    canonical operand when containment makes that exact without subdivision;
    other nonconvex results deliberately fall back to ``BooleanRegion``.
    """
    from .region_structural import structural_disjoint, structural_subset
    from .standard_regions import BooleanRegion

    if operation == "intersection":
        result = polyhedral_intersection(left, right)
        return result if result is not None else BooleanRegion.intersection(left, right)
    if operation == "union":
        if structural_subset(left, right) is True:
            return right
        if structural_subset(right, left) is True:
            return left
        return BooleanRegion.union(
            left, right, assume_disjoint=structural_disjoint(left, right) is True
        )
    if operation == "difference":
        if structural_disjoint(left, right) is True:
            return left
        return BooleanRegion.difference(left, right)
    if operation == "symmetric_difference":
        return BooleanRegion.symmetric_difference(left, right)
    raise ValueError(f"unsupported polyhedral Boolean operation: {operation!r}")


__all__ = ["polyhedral_intersection", "polyhedral_boolean"]
