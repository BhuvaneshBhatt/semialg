from __future__ import annotations

import itertools
from collections.abc import Mapping
from dataclasses import dataclass

import sympy as sp

from ..cad_algorithms.cells import CylindricalSolutionCell, extract_cylindrical_solution
from ..connectivity import CADConnectedComponent, extract_cad_connectivity
from ..derived_geometry import is_compact
from ..normalization import normalize_formula, normalize_problem_variables
from ..standard_regions import (
    BoxRegion,
    ParallelepipedRegion,
    PolygonRegion,
    Polytope,
    Simplex,
    SimplexRegion,
    StandardRegion,
)


@dataclass(frozen=True)
class SimplicialComplex:
    """Finite exact simplicial complex represented by its maximal simplices."""

    vertices: tuple[tuple[sp.Expr, ...], ...]
    maximal_simplices: tuple[tuple[int, ...], ...]

    @property
    def dimension(self) -> int:
        if not self.maximal_simplices:
            return -1
        return max(len(simplex) - 1 for simplex in self.maximal_simplices)

    @property
    def simplex_count(self) -> int:
        return len(self.maximal_simplices)

    def simplices(self) -> tuple[Simplex, ...]:
        """Return the maximal simplices as canonical geometry objects."""

        return tuple(
            Simplex(tuple(self.vertices[i] for i in simplex)) for simplex in self.maximal_simplices
        )


@dataclass(frozen=True)
class SemialgebraicTriangulation:
    """Certified finite triangulation for a currently supported semialgebraic set.

    For polyhedral inputs the realization of ``complex`` is the region itself.
    The implementation does not claim a straight-simplex realization for curved
    semialgebraic sets; those require a genuine semialgebraic homeomorphism.
    """

    variables: tuple[sp.Symbol, ...]
    region_formula: sp.Expr
    complex: SimplicialComplex
    construction: str
    certified: bool = True

    @property
    def dimension(self) -> int:
        return self.complex.dimension


def _box_vertices(region: BoxRegion) -> tuple[tuple[sp.Expr, ...], ...]:
    return tuple(tuple(choice) for choice in itertools.product(*region.bounds))


def _as_polytope(region: StandardRegion) -> Polytope | None:
    if isinstance(region, Polytope):
        return region
    if isinstance(region, (SimplexRegion, PolygonRegion, ParallelepipedRegion)):
        return Polytope.from_region(region)
    if isinstance(region, BoxRegion):
        return Polytope(_box_vertices(region))
    return None


def _pulling_triangulation(polytope: Polytope) -> tuple[tuple[int, ...], ...]:
    vertices = polytope.vertices
    dimension = polytope.dimension()
    if dimension < 0:
        return ()
    if dimension == 0:
        return ((0,),)

    lattice = polytope.face_lattice()
    face_dimensions = {frozenset(face.vertex_indices): face.dimension for face in lattice.faces}
    whole = frozenset(range(len(vertices)))
    face_dimensions[whole] = dimension

    def triangulate(face_vertices: frozenset[int]) -> tuple[tuple[int, ...], ...]:
        face_dimension = face_dimensions[face_vertices]
        ordered = tuple(sorted(face_vertices))
        if len(ordered) == face_dimension + 1:
            return (ordered,)
        anchor = ordered[0]
        facets = [
            fv
            for fv, dim in face_dimensions.items()
            if dim == face_dimension - 1 and fv < face_vertices and anchor not in fv
        ]
        # Keep only maximal codimension-one subfaces of the current face.
        facets = [fv for fv in facets if not any(fv < other for other in facets)]
        simplices: list[tuple[int, ...]] = []
        for facet in sorted(facets, key=lambda item: tuple(sorted(item))):
            for simplex in triangulate(facet):
                simplices.append(tuple(sorted((anchor, *simplex))))
        return tuple(dict.fromkeys(simplices))

    return triangulate(whole)


def _polyhedral_triangulation(region: StandardRegion) -> SimplicialComplex | None:
    if isinstance(region, SimplexRegion):
        vertices = tuple(region.vertices)
        return SimplicialComplex(vertices, (tuple(range(len(vertices))),))
    if isinstance(region, PolygonRegion):
        pieces = region.triangulation()
        vertices = tuple(dict.fromkeys(vertex for piece in pieces for vertex in piece.vertices))
        positions = {vertex: i for i, vertex in enumerate(vertices)}
        simplices = tuple(tuple(positions[v] for v in piece.vertices) for piece in pieces)
        return SimplicialComplex(vertices, simplices)
    polytope = _as_polytope(region)
    if polytope is None:
        return None
    return SimplicialComplex(tuple(polytope.vertices), _pulling_triangulation(polytope))


def _real_line_triangulation(formula: sp.Expr, variable: sp.Symbol) -> SimplicialComplex:
    solution = extract_cylindrical_solution(formula, (variable,), selected_only=True)
    if solution.bounded is not True:
        raise NotImplementedError(
            "finite straight-simplex triangulation currently requires compact input"
        )
    vertices: list[tuple[sp.Expr, ...]] = []
    segments: list[tuple[int, ...]] = []

    def position(value: sp.Expr) -> int:
        point = (sp.simplify(value),)
        if point not in vertices:
            vertices.append(point)
        return vertices.index(point)

    for cell in solution.cells:
        level = cell.levels[0]
        if level.is_section:
            position(level.lower)
            continue
        if level.lower in (-sp.oo, sp.oo) or level.upper in (-sp.oo, sp.oo):
            raise NotImplementedError(
                "unbounded one-dimensional cells require a non-linear model triangulation"
            )
        left = position(level.lower)
        right = position(level.upper)
        if left != right:
            segments.append(tuple(sorted((left, right))))
    maximal: list[tuple[int, ...]] = list(dict.fromkeys(segments))
    used = {i for segment in segments for i in segment}
    maximal.extend((i,) for i in range(len(vertices)) if i not in used)
    return SimplicialComplex(tuple(vertices), tuple(maximal))


def triangulate_region(region, variables=None) -> SemialgebraicTriangulation:
    """Return an exact finite triangulation for supported compact regions.

    Convex/polyhedral canonical regions are triangulated by a pulling
    triangulation of their face lattice. Compact semialgebraic subsets of the
    real line are triangulated from exact CAD sections and sectors. General
    curved semialgebraic sets require construction of the homeomorphism in the
    semialgebraic triangulation theorem and are not approximated here.
    """

    if isinstance(region, StandardRegion):
        complex_ = _polyhedral_triangulation(region)
        if complex_ is not None:
            vars_ = tuple(sp.symbols(f"x0:{region.ambient_dimension()}", real=True))
            return SemialgebraicTriangulation(
                vars_, region.as_formula(vars_), complex_, "polyhedral-pulling"
            )
    formula = normalize_formula(region)
    vars_ = normalize_problem_variables(variables, formula)
    if not is_compact(formula, vars_):
        raise NotImplementedError(
            "triangulate_region currently requires compact non-polyhedral input"
        )
    if len(vars_) == 1:
        return SemialgebraicTriangulation(
            vars_, formula, _real_line_triangulation(formula, vars_[0]), "real-line-cad"
        )
    raise NotImplementedError(
        "general semialgebraic triangulation requires a certified semialgebraic homeomorphism"
    )


@dataclass(frozen=True)
class DimensionStratum:
    """Union of selected CAD cells having one Euclidean dimension."""

    dimension: int
    formula: sp.Expr
    cells: tuple[CylindricalSolutionCell, ...]


@dataclass(frozen=True)
class DimensionDecomposition:
    """Exact CAD-derived decomposition of a region by cell dimension."""

    variables: tuple[sp.Symbol, ...]
    formula: sp.Expr
    dimension: int
    strata: tuple[DimensionStratum, ...]


def dimension_strata(region, variables=None) -> DimensionDecomposition:
    """Decompose a semialgebraic set into unions of CAD cells by dimension."""

    formula = normalize_formula(region)
    vars_ = normalize_problem_variables(variables, formula)
    solution = extract_cylindrical_solution(formula, vars_, selected_only=True)
    by_dimension: dict[int, list[CylindricalSolutionCell]] = {}
    for cell in solution.cells:
        by_dimension.setdefault(cell.dimension, []).append(cell)
    strata = []
    for dimension in sorted(by_dimension):
        cells = tuple(by_dimension[dimension])
        pieces = tuple(cell.as_formula() for cell in cells)
        stratum_formula = sp.Or(*pieces) if len(pieces) > 1 else pieces[0]
        strata.append(DimensionStratum(dimension, stratum_formula, cells))
    exact_dimension = max(by_dimension, default=-1)
    return DimensionDecomposition(vars_, formula, exact_dimension, tuple(strata))


@dataclass(frozen=True)
class ConnectedComponentDecomposition:
    """Exact semialgebraic connected components induced by CAD adjacency."""

    variables: tuple[sp.Symbol, ...]
    formula: sp.Expr
    components: tuple[CADConnectedComponent, ...]

    @property
    def count(self) -> int:
        return len(self.components)

    @property
    def sample_points(self) -> tuple[Mapping[sp.Symbol, sp.Expr], ...]:
        return tuple(component.sample_point() for component in self.components)


def component_decomposition(region, variables=None) -> ConnectedComponentDecomposition:
    """Return exact connected-component formulas, dimensions, and samples."""

    formula = normalize_formula(region)
    vars_ = normalize_problem_variables(variables, formula)
    graph = extract_cad_connectivity(formula, vars_)
    return ConnectedComponentDecomposition(vars_, formula, graph.components)


@dataclass(frozen=True)
class HardtFiberPiece:
    """One graph or band in the one-dimensional fiber of a CAD stratum."""

    kind: str
    lower: sp.Expr
    upper: sp.Expr
    sample: sp.Expr
    stack_position: int

    @property
    def dimension(self) -> int:
        return 0 if self.kind == "section" else 1

    def model_coordinate(self, value: sp.Expr) -> sp.Expr:
        """Map a fiber value to a fixed semialgebraic model coordinate."""

        value = sp.sympify(value)
        if self.kind == "section":
            return sp.Integer(0)
        if self.lower == -sp.oo and self.upper == sp.oo:
            return value
        if self.lower == -sp.oo:
            return sp.simplify(value - self.upper)
        if self.upper == sp.oo:
            return sp.simplify(value - self.lower)
        return sp.simplify((value - self.lower) / (self.upper - self.lower))


@dataclass(frozen=True)
class HardtStratum:
    """One certified base cell over which a coordinate projection is trivial."""

    index: tuple[int, ...]
    parameters: tuple[sp.Symbol, ...]
    fiber_variable: sp.Symbol
    condition: sp.Expr
    sample: Mapping[sp.Symbol, sp.Expr]
    base_dimension: int
    fiber_pieces: tuple[HardtFiberPiece, ...]

    @property
    def fiber_dimension(self) -> int:
        return max((piece.dimension for piece in self.fiber_pieces), default=-1)

    @property
    def fiber_piece_count(self) -> int:
        return len(self.fiber_pieces)

    @property
    def fiber_component_count(self) -> int:
        if not self.fiber_pieces:
            return 0
        positions = sorted(piece.stack_position for piece in self.fiber_pieces)
        return 1 + sum(
            right - left > 1 for left, right in zip(positions, positions[1:], strict=False)
        )


@dataclass(frozen=True)
class HardtTrivialization:
    """Finite certified Hardt partition for projection with one fiber variable.

    Each nonempty stratum is a CAD base cell. Its selected final stack consists
    of delineable sections and bands, yielding the graph/band product
    trivialization used in Coste's proof for coordinate projections.
    """

    formula: sp.Expr
    parameters: tuple[sp.Symbol, ...]
    fiber_variable: sp.Symbol
    strata: tuple[HardtStratum, ...]
    image_condition: sp.Expr
    empty_fiber_condition: sp.Expr
    certified: bool = True
    method: str = "cylindrical-graph-band"

    def fiber_dimensions(self) -> tuple[int, ...]:
        return tuple(stratum.fiber_dimension for stratum in self.strata)


def hardt_trivialization(region, parameters, fiber_variable) -> HardtTrivialization:
    """Construct Hardt strata for a coordinate projection with 1D fibers.

    The variable order is ``parameters`` followed by ``fiber_variable``. An
    adapted CAD makes each selected stack over a base cell a fixed ordered list
    of sections and bands. Those pieces are semialgebraically trivial over the
    base cell, exactly as in the coordinate-projection construction preceding
    Coste's statement of Hardt's theorem.
    """

    formula = normalize_formula(region)
    params = tuple(sp.sympify(p) for p in parameters)
    fiber = sp.sympify(fiber_variable)
    if not isinstance(fiber, sp.Symbol):
        raise TypeError("fiber_variable must be a Symbol")
    if not params or not all(isinstance(p, sp.Symbol) for p in params):
        raise ValueError("parameters must contain at least one Symbol")
    if fiber in params:
        raise ValueError("fiber_variable must be distinct from the parameters")
    variables = (*params, fiber)
    solution = extract_cylindrical_solution(formula, variables, selected_only=True)
    grouped: dict[tuple[int, ...], list[CylindricalSolutionCell]] = {}
    for cell in solution.cells:
        grouped.setdefault(cell.index[:-1], []).append(cell)

    strata: list[HardtStratum] = []
    for prefix in sorted(grouped):
        cells = sorted(grouped[prefix], key=lambda cell: cell.index)
        first = cells[0]
        base_levels = first.levels[:-1]
        condition_parts = tuple(level.as_formula() for level in base_levels)
        condition = sp.And(*condition_parts) if condition_parts else sp.true
        sample = {p: first.sample[p] for p in params}
        pieces = tuple(
            HardtFiberPiece(
                cell.levels[-1].kind,
                cell.levels[-1].lower,
                cell.levels[-1].upper,
                cell.levels[-1].sample,
                cell.index[-1],
            )
            for cell in cells
        )
        strata.append(
            HardtStratum(
                prefix,
                params,
                fiber,
                condition,
                sample,
                sum(level.dimension for level in base_levels),
                pieces,
            )
        )
    conditions = tuple(stratum.condition for stratum in strata)
    image_condition = sp.Or(*conditions) if conditions else sp.false
    empty_fiber_condition = sp.Not(image_condition) if image_condition is not sp.false else sp.true
    return HardtTrivialization(
        formula, params, fiber, tuple(strata), image_condition, empty_fiber_condition
    )


__all__ = [
    "SimplicialComplex",
    "SemialgebraicTriangulation",
    "triangulate_region",
    "simplicial_betti_numbers",
    "triangulation_betti_numbers",
    "DimensionStratum",
    "DimensionDecomposition",
    "dimension_strata",
    "ConnectedComponentDecomposition",
    "component_decomposition",
    "HardtFiberPiece",
    "HardtStratum",
    "HardtTrivialization",
    "hardt_trivialization",
]


def simplicial_betti_numbers(complex_: SimplicialComplex) -> tuple[int, ...]:
    """Return exact Betti numbers of a finite simplicial complex over ``QQ``.

    All faces of the maximal simplices are generated canonically.  Boundary
    ranks are computed exactly, so this routine introduces no numerical
    tolerance into topology calculations.
    """
    dimension = complex_.dimension
    if dimension < 0:
        return ()
    faces: list[tuple[tuple[int, ...], ...]] = []
    for degree in range(dimension + 1):
        size = degree + 1
        level = {
            tuple(face)
            for simplex in complex_.maximal_simplices
            for face in itertools.combinations(sorted(simplex), size)
        }
        faces.append(tuple(sorted(level)))

    ranks = [0] * (dimension + 2)
    for degree in range(1, dimension + 1):
        rows = {face: i for i, face in enumerate(faces[degree - 1])}
        matrix = sp.zeros(len(rows), len(faces[degree]))
        for column, simplex in enumerate(faces[degree]):
            for omitted in range(len(simplex)):
                face = simplex[:omitted] + simplex[omitted + 1 :]
                matrix[rows[face], column] = -1 if omitted % 2 else 1
        ranks[degree] = int(matrix.rank())

    return tuple(
        len(faces[degree]) - ranks[degree] - ranks[degree + 1] for degree in range(dimension + 1)
    )


def triangulation_betti_numbers(region, variables=None) -> tuple[int, ...]:
    """Return all Betti numbers when :func:`triangulate_region` is supported."""
    return simplicial_betti_numbers(triangulate_region(region, variables).complex)
