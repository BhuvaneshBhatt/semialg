"""Local algebraic-geometry utilities for real polynomial varieties."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from itertools import combinations

import sympy as sp

from ._zero_testing import certified_zero
from .algebraic import (
    FractionFieldZeroDimensionalCertificate,
    FractionFieldZeroDimensionalResult,
    GTZContractedComponentCertificate,
    GTZExecutionPlan,
    GTZNodeCertificate,
    GTZPrimaryDecompositionCertificate,
    GTZPrimaryDecompositionResult,
    HilbertData,
    IndependentLocalizationCertificate,
    IndependentLocalizationResult,
    LocalizationContractionCertificate,
    LocalizationContractionResult,
    LocalPrimaryFactor,
    SaturationStabilizationCertificate,
    SaturationStabilizationResult,
    ZeroDimensionalPrimaryCertificate,
    ZeroDimensionalPrimaryComponent,
    ZeroDimensionalPrimaryResult,
    certified_independent_localization,
    clear_gtz_caches,
    contract_localized_ideal,
    gtz_cache_info,
    gtz_primary_decomposition,
    hilbert_function,
    hilbert_polynomial,
    ideal_degree,
    ideal_hilbert_data,
    plan_gtz_primary_decomposition,
    saturation_stabilization,
    verify_fraction_field_zero_dimensional_certificate,
    verify_gtz_primary_decomposition_certificate,
    verify_independent_localization_certificate,
    verify_localization_contraction_certificate,
    verify_saturation_stabilization_certificate,
    verify_zero_dimensional_primary_certificate,
    zero_dimensional_primary_decomposition,
)
from .algebraic_decomposition import (
    AlgebraicLocusPiece,
    AssociatedPrimesResult,
    DecompositionCertificate,
    DecompositionPieceCertificate,
    EquidimensionalDecomposition,
    MinimalPrimeDecompositionCertificate,
    PrimaryComponent,
    PrimaryDecompositionCertificate,
    PrimaryDecompositionResult,
    RadicalComponent,
    RadicalIdealCertificate,
    RadicalIdealResult,
    RadicalMinimalPrimeDecomposition,
    RegularChainComponent,
    RegularChainDecomposition,
    RegularChainDecompositionCertificate,
    RegularChainSplitCertificate,
    TriangularPrimalityCertificate,
    TriangularPrimalityStage,
    associated_primes,
    certified_radical_minimal_prime_decomposition,
    certify_triangular_primality,
    equidimensional_decomposition,
    primary_decomposition,
    radical_ideal,
    recursive_regular_chain_decomposition,
    verify_decomposition_certificate,
    verify_minimal_prime_decomposition_certificate,
    verify_primary_decomposition_certificate,
    verify_radical_ideal_certificate,
    verify_regular_chain_decomposition_certificate,
    verify_triangular_primality_certificate,
)
from .internal_symbols import collision_free_real_symbols, fresh_dummy
from .normalization import normalize_point, normalize_problem_variables, normalize_symbol_sequence
from .optimization_geometry import polynomial_locus_dimension


def _equations(equations) -> tuple[sp.Expr, ...]:
    if isinstance(equations, (sp.Expr, sp.Equality)):
        equations = (equations,)
    out = []
    for eq in equations:
        expr = sp.sympify(eq)
        if isinstance(expr, sp.Equality):
            expr = sp.expand(expr.lhs - expr.rhs)
        out.append(sp.expand(expr))
    return tuple(out)


@dataclass(frozen=True)
class PolynomialMapImplicitizationResult:
    """Exact elimination result for the Zariski closure of a polynomial map image."""

    parameter_variables: tuple[sp.Symbol, ...]
    image_variables: tuple[sp.Symbol, ...]
    mapping: tuple[sp.Expr, ...]
    domain_equations: tuple[sp.Expr, ...]
    equations: tuple[sp.Expr, ...]
    groebner_basis: tuple[sp.Expr, ...]

    @property
    def formula(self) -> sp.Expr:
        """Return the algebraic image closure as a conjunction of equations."""
        if not self.equations:
            return sp.true
        return sp.And(*(sp.Eq(equation, 0) for equation in self.equations))


@dataclass(frozen=True)
class ZariskiClosureResult:
    """Exact equations and ambient coordinates of a polynomial image closure."""

    variables: tuple[sp.Symbol, ...]
    equations: tuple[sp.Expr, ...]
    formula: sp.Expr
    implicitization: PolynomialMapImplicitizationResult


def _polynomial_map_data(mapping, parameters):
    maps = (
        (sp.sympify(mapping),) if isinstance(mapping, sp.Expr) else tuple(map(sp.sympify, mapping))
    )
    if not maps:
        raise ValueError("mapping must contain at least one coordinate polynomial")
    params = normalize_symbol_sequence(parameters)
    if not params:
        raise ValueError("parameters must contain at least one symbol")
    if len(set(params)) != len(params):
        raise ValueError("parameters must be distinct")
    parameter_set = set(params)
    for coordinate in maps:
        if coordinate.free_symbols - parameter_set:
            raise ValueError("mapping contains symbols that are not declared parameters")
        try:
            sp.Poly(coordinate, *params)
        except sp.PolynomialError as exc:
            raise ValueError("mapping coordinates must be polynomial in the parameters") from exc
    return maps, params


def implicitize_polynomial_map(
    mapping,
    parameters,
    *,
    image_variables=None,
    domain_equations=(),
    return_result: bool = False,
):
    """Implicitize a polynomial map by exact Groebner elimination.

    ``domain_equations`` may restrict the parameter space to an algebraic
    subvariety.  Inequalities are intentionally not accepted here: an arbitrary
    semialgebraic restriction can have the same Zariski closure as its ambient
    parameter variety, and deciding density is a separate problem.

    By default the reduced elimination generators in the image coordinates are
    returned.  ``return_result=True`` also exposes the graph Groebner basis and
    the exact closure formula.
    """

    maps, params = _polynomial_map_data(mapping, parameters)
    domain = _equations(domain_equations) if domain_equations else ()
    for equation in domain:
        if equation.free_symbols - set(params):
            raise ValueError("domain_equations contain symbols that are not declared parameters")
        try:
            sp.Poly(equation, *params)
        except sp.PolynomialError as exc:
            raise ValueError("domain_equations must be polynomial in the parameters") from exc

    if image_variables is None:
        targets = collision_free_real_symbols(
            "image", len(maps), sp.Tuple(*maps), sp.Tuple(*domain), tuple(params)
        )
    else:
        targets = normalize_symbol_sequence(image_variables)
        if len(targets) != len(maps):
            raise ValueError("image_variables must have the same length as mapping")
        if len(set(targets)) != len(targets):
            raise ValueError("image_variables must be distinct")
        if set(targets) & set(params):
            raise ValueError("image_variables must be distinct from parameters")

    graph_equations = (
        *domain,
        *(target - coordinate for target, coordinate in zip(targets, maps, strict=True)),
    )
    all_variables = (*params, *targets)
    basis = _groebner_exact(graph_equations, all_variables)
    basis_expressions = tuple(sp.expand(poly.as_expr()) for poly in basis.polys)
    eliminated = tuple(
        expression
        for expression in basis_expressions
        if not (expression.free_symbols & set(params))
    )
    if eliminated:
        image_basis = _groebner_exact(eliminated, targets)
        equations = tuple(sp.expand(poly.as_expr()) for poly in image_basis.polys)
        if equations == (sp.Integer(1),):
            equations = (sp.Integer(1),)
    else:
        equations = ()

    result = PolynomialMapImplicitizationResult(
        params, targets, maps, domain, equations, basis_expressions
    )
    return result if return_result else result.equations


def zariski_closure(
    mapping,
    parameters,
    *,
    image_variables=None,
    domain_equations=(),
    return_result: bool = False,
):
    """Return the Zariski closure of a polynomial-map image.

    The closure is represented by the elimination ideal of the graph of the
    map (and optional algebraic parameter-domain equations).  The default
    result is a Boolean conjunction of exact polynomial equalities in the image
    coordinates.
    """

    implicit = implicitize_polynomial_map(
        mapping,
        parameters,
        image_variables=image_variables,
        domain_equations=domain_equations,
        return_result=True,
    )
    result = ZariskiClosureResult(
        implicit.image_variables, implicit.equations, implicit.formula, implicit
    )
    return result if return_result else result.formula


@dataclass(frozen=True)
class ReducedAlgebraicVariety:
    """Certified reduced presentation of an affine algebraic variety."""

    variables: tuple[sp.Symbol, ...]
    source_equations: tuple[sp.Expr, ...]
    equations: tuple[sp.Expr, ...]
    radical_certificate: RadicalIdealCertificate

    @property
    def formula(self) -> sp.Expr:
        if not self.equations:
            return sp.true
        return sp.And(*(sp.Eq(equation, 0) for equation in self.equations))


@dataclass(frozen=True)
class IrreducibleAlgebraicComponent:
    """One certified minimal-prime component of a reduced variety."""

    variables: tuple[sp.Symbol, ...]
    equations: tuple[sp.Expr, ...]
    dimension: int
    degree: int
    prime_method: str

    @property
    def formula(self) -> sp.Expr:
        if not self.equations:
            return sp.true
        return sp.And(*(sp.Eq(equation, 0) for equation in self.equations))


@dataclass(frozen=True)
class ReducedComponentSingularLocus:
    """Singular locus attached to one certified irreducible component."""

    component: IrreducibleAlgebraicComponent
    formula: sp.Expr


@dataclass(frozen=True)
class MinimalPrimeIntersection:
    """Exact intersection of two or more certified minimal-prime components."""

    component_indices: tuple[int, ...]
    equations: tuple[sp.Expr, ...]
    dimension: int | None
    formula: sp.Expr


@dataclass(frozen=True)
class LocalDimensionStratum:
    """Constructible locus on which the algebraic local dimension is constant."""

    dimension: int
    component_indices: tuple[int, ...]
    formula: sp.Expr


@dataclass(frozen=True)
class SingularGeometryStratum:
    """A disjoint constructible singular-geometry stratum.

    ``component_indices`` records exactly the minimal-prime components containing
    the stratum.  ``intrinsically_singular`` distinguishes singularities of one
    or more reduced components from singularities caused solely by component
    intersection.
    """

    component_indices: tuple[int, ...]
    local_dimension: int
    intrinsically_singular: bool
    component_intersection: bool
    formula: sp.Expr


@dataclass(frozen=True)
class BranchTangentGeometry:
    """Exact local geometry of one irreducible branch at a point."""

    component_index: int
    component: IrreducibleAlgebraicComponent
    tangent_cone: TangentConeResult
    tangent_space: TangentSpaceResult
    tangent_dimension_excess: int
    multiplicity: int
    tangent_cone_degree: int


@dataclass(frozen=True)
class ComponentIntersectionGeometry:
    """Local geometry of the intersection of the branches through a point."""

    component_indices: tuple[int, ...]
    equations: tuple[sp.Expr, ...]
    dimension: int
    tangent_cone: TangentConeResult
    tangent_space: TangentSpaceResult
    tangent_dimension_excess: int
    transverse: bool
    multiplicity: int
    tangent_cone_degree: int


@dataclass(frozen=True)
class LocalBranchGeometry:
    """Certified branch, tangent, transversality, and multiplicity data at a point."""

    variables: tuple[sp.Symbol, ...]
    point: Mapping[sp.Symbol, sp.Expr]
    component_indices: tuple[int, ...]
    branches: tuple[BranchTangentGeometry, ...]
    intersection: ComponentIntersectionGeometry | None
    local_dimension: int
    singular: bool


@dataclass(frozen=True)
class StratifiedSingularGeometry:
    """Certified component/intersection/local-dimension description of a variety."""

    variables: tuple[sp.Symbol, ...]
    components: tuple[IrreducibleAlgebraicComponent, ...]
    component_singular_loci: tuple[ReducedComponentSingularLocus, ...]
    minimal_prime_intersections: tuple[MinimalPrimeIntersection, ...]
    local_dimension_strata: tuple[LocalDimensionStratum, ...]
    singular_strata: tuple[SingularGeometryStratum, ...]


def certified_radicalization(
    equations, variables=None, *, max_pieces: int | None = None
) -> ReducedAlgebraicVariety:
    """Return a certified reduced presentation ``sqrt(<equations>)``.

    The function refuses to return uncertified generators.  This makes the
    reduced ideal, rather than a possibly non-radical input presentation, the
    boundary for geometric operations such as Jacobian singular-locus tests.
    """
    eqs = _equations(equations)
    vars_ = normalize_problem_variables(variables, sp.Tuple(*eqs))
    result = radical_ideal(eqs, vars_, max_pieces=max_pieces)
    if not result.complete or result.certificate is None:
        raise NotImplementedError("could not certify the radical ideal")
    if not verify_radical_ideal_certificate(result.certificate):
        raise RuntimeError("radical-ideal certificate failed replay verification")
    return ReducedAlgebraicVariety(vars_, eqs, result.generators, result.certificate)


def irreducible_components(
    equations, variables=None, *, max_pieces: int | None = None
) -> tuple[IrreducibleAlgebraicComponent, ...]:
    """Return all certified irreducible components of an affine variety.

    Components are returned only when semialg certifies a complete,
    irredundant minimal-prime decomposition of the radical ideal.
    """
    eqs = _equations(equations)
    vars_ = normalize_problem_variables(variables, sp.Tuple(*eqs))
    decomposition = certified_radical_minimal_prime_decomposition(eqs, vars_, max_pieces=max_pieces)
    if decomposition.certificate is None or not verify_minimal_prime_decomposition_certificate(
        decomposition.certificate
    ):
        raise RuntimeError("minimal-prime certificate failed replay verification")
    if not decomposition.radical_complete or not decomposition.minimal_primes_complete:
        raise NotImplementedError("could not certify a complete minimal-prime decomposition")
    return tuple(
        IrreducibleAlgebraicComponent(
            vars_,
            component.equations,
            component.dimension,
            component.degree,
            component.prime_method,
        )
        for component in decomposition.components
        if component.prime and component.prime_method is not None
    )


def reduced_component_singular_loci(
    equations, variables=None, *, max_pieces: int | None = None
) -> tuple[ReducedComponentSingularLocus, ...]:
    """Return singular loci computed separately on certified reduced components."""
    return tuple(
        ReducedComponentSingularLocus(
            component,
            singular_locus(
                component.equations,
                component.variables,
                codimension=len(component.variables) - component.dimension,
                reduce=False,
            ),
        )
        for component in irreducible_components(equations, variables, max_pieces=max_pieces)
    )


def _formula_union(formulas: Sequence[sp.Expr]) -> sp.Expr:
    formulas = tuple(formulas)
    return sp.Or(*formulas) if formulas else sp.false


def _component_intersection_equations(
    components: Sequence[IrreducibleAlgebraicComponent],
    indices: Sequence[int],
) -> tuple[sp.Expr, ...]:
    variables = components[0].variables
    generators = tuple(equation for index in indices for equation in components[index].equations)
    if not generators:
        return ()
    basis = _groebner_exact(generators, variables)
    return tuple(sp.expand(poly.as_expr()) for poly in basis.polys)


def minimal_prime_intersections(
    equations,
    variables=None,
    *,
    max_pieces: int | None = None,
    max_order: int | None = None,
) -> tuple[MinimalPrimeIntersection, ...]:
    """Return exact intersections among certified minimal-prime components.

    The intersection of varieties is represented by the sum of their prime
    ideals.  Reduced Groebner generators make the result deterministic.
    ``max_order`` can cap the combinatorial expansion for varieties with many
    irreducible components.
    """
    components = irreducible_components(equations, variables, max_pieces=max_pieces)
    if len(components) < 2:
        return ()
    stop = len(components) if max_order is None else min(max_order, len(components))
    if stop < 2:
        return ()
    results = []
    for order in range(2, stop + 1):
        for indices in combinations(range(len(components)), order):
            generators = _component_intersection_equations(components, indices)
            impossible = any(
                sp.Poly(g, *components[0].variables).total_degree() == 0 and g != 0
                for g in generators
            )
            if impossible:
                dimension = -1
                formula = sp.false
            else:
                dimension = polynomial_locus_dimension(generators, components[0].variables)
                formula = sp.And(*(sp.Eq(g, 0) for g in generators)) if generators else sp.true
            results.append(MinimalPrimeIntersection(indices, generators, dimension, formula))
    return tuple(results)


def local_dimension_strata(
    equations, variables=None, *, max_pieces: int | None = None
) -> tuple[LocalDimensionStratum, ...]:
    """Partition a reduced variety into constructible constant-local-dimension loci.

    At a point of a reduced affine variety, local dimension is the maximum of
    the dimensions of the irreducible components through that point.  The
    returned formulas implement that statement directly, so intersections of a
    lower-dimensional component with a higher-dimensional one are assigned the
    higher local dimension rather than being reported as artificial jumps.
    """
    components = irreducible_components(equations, variables, max_pieces=max_pieces)
    dimensions = sorted({component.dimension for component in components}, reverse=True)
    strata = []
    for dimension in dimensions:
        same = tuple(i for i, c in enumerate(components) if c.dimension == dimension)
        higher = tuple(i for i, c in enumerate(components) if c.dimension > dimension)
        formula = sp.And(
            _formula_union(components[i].formula for i in same),
            sp.Not(_formula_union(components[i].formula for i in higher)),
        )
        strata.append(LocalDimensionStratum(dimension, same, formula))
    return tuple(strata)


def stratified_singular_geometry(
    equations,
    variables=None,
    *,
    max_pieces: int | None = None,
    max_intersection_order: int | None = None,
) -> StratifiedSingularGeometry:
    """Return a certified stratified description of singular algebraic geometry.

    The result separates four notions that a global Jacobian singular locus
    conflates: certified minimal-prime components, their mutual intersections,
    each component's intrinsic singular locus, and local-dimension strata.
    ``singular_strata`` is a disjoint constructible refinement by exact component
    membership and intrinsic singular status.  Membership in two or more
    components is recorded as a component-intersection singularity even when all
    participating components are individually smooth.
    """
    components = irreducible_components(equations, variables, max_pieces=max_pieces)
    if not components:
        vars_ = normalize_problem_variables(variables, sp.Tuple(*_equations(equations)))
        return StratifiedSingularGeometry(vars_, (), (), (), (), ())
    vars_ = components[0].variables
    component_singularities = reduced_component_singular_loci(
        equations, vars_, max_pieces=max_pieces
    )
    intersections = minimal_prime_intersections(
        equations, vars_, max_pieces=max_pieces, max_order=max_intersection_order
    )
    dimension_strata = local_dimension_strata(equations, vars_, max_pieces=max_pieces)

    strata = []
    component_count = len(components)
    for membership_size in range(1, component_count + 1):
        for indices in combinations(range(component_count), membership_size):
            included = _formula_union(components[i].formula for i in indices)
            # Exact membership requires every selected component and excludes all others.
            included = sp.And(*(components[i].formula for i in indices))
            excluded = _formula_union(
                components[i].formula for i in range(component_count) if i not in indices
            )
            membership_formula = sp.And(included, sp.Not(excluded))
            intrinsic_formula = _formula_union(component_singularities[i].formula for i in indices)
            local_dimension = max(components[i].dimension for i in indices)
            crossing = membership_size >= 2
            # Split each exact-membership locus into intrinsic-singular and
            # componentwise-regular pieces.  Empty formulas are harmless and keep
            # the construction purely exact without invoking real feasibility.
            intrinsic_stratum = sp.And(membership_formula, intrinsic_formula)
            if intrinsic_stratum is not sp.false:
                strata.append(
                    SingularGeometryStratum(
                        indices,
                        local_dimension,
                        True,
                        crossing,
                        intrinsic_stratum,
                    )
                )
            if crossing:
                strata.append(
                    SingularGeometryStratum(
                        indices,
                        local_dimension,
                        False,
                        True,
                        sp.And(membership_formula, sp.Not(intrinsic_formula)),
                    )
                )
    return StratifiedSingularGeometry(
        vars_,
        components,
        component_singularities,
        intersections,
        dimension_strata,
        tuple(strata),
    )


def singular_locus(
    equations,
    variables=None,
    *,
    codimension: int | None = None,
    reduce: bool = True,
    max_pieces: int | None = None,
) -> sp.Expr:
    """Return equations defining the singular locus of an algebraic variety.

    The Jacobian criterion is applied to a certified radical presentation by
    default, so nilpotent multiplicities such as ``x**2 = 0`` cannot create
    spurious geometric singularities.  Set ``reduce=False`` only when the
    caller already holds a certified reduced component.
    """

    if isinstance(equations, PolynomialMapImplicitizationResult):
        if variables is not None:
            supplied = normalize_symbol_sequence(variables)
            if supplied != equations.image_variables:
                raise ValueError("variables do not match implicitized image variables")
        vars_ = equations.image_variables
        eqs = equations.equations
    else:
        eqs = _equations(equations)
        vars_ = normalize_problem_variables(variables, sp.Tuple(*eqs))
    if reduce:
        reduced = certified_radicalization(eqs, vars_, max_pieces=max_pieces)
        eqs = reduced.equations
    if codimension is None:
        dim = polynomial_locus_dimension(eqs, vars_)
        if dim is None:
            raise NotImplementedError(
                "could not infer variety dimension; pass codimension explicitly"
            )
        codimension = max(0, len(vars_) - dim)
    if codimension < 0 or codimension > len(vars_):
        raise ValueError("invalid codimension")
    base = [sp.Eq(eq, 0) for eq in eqs]
    if codimension == 0:
        return sp.false
    jac = sp.Matrix([[sp.diff(eq, v) for v in vars_] for eq in eqs])
    if codimension > min(jac.rows, jac.cols):
        return sp.And(*base)
    minors = []
    for rows in combinations(range(jac.rows), codimension):
        for cols in combinations(range(jac.cols), codimension):
            minors.append(sp.expand(jac.extract(rows, cols).det()))
    if not minors:
        return sp.And(*base)
    return sp.And(*base, *(sp.Eq(m, 0) for m in minors))


def _point_satisfies(equations: Sequence[sp.Expr], point: Mapping[sp.Symbol, sp.Expr]) -> bool:
    return all(certified_zero(equation.subs(point)) is True for equation in equations)


def _tangent_cone_degree(cone: TangentConeResult) -> int:
    if not cone.initial_forms:
        return 1
    return ideal_degree(cone.initial_forms, cone.direction_variables)


def local_branch_geometry(
    equations,
    point,
    variables=None,
    *,
    max_pieces: int | None = None,
) -> LocalBranchGeometry:
    """Return exact local branch geometry at a point of a reduced variety.

    Each minimal-prime branch through the point is analyzed independently.
    ``tangent_dimension_excess`` is the excess of Zariski tangent dimension over
    the branch (or intersection) dimension.  Local multiplicity is computed as
    the exact degree of the ideal-theoretic tangent cone.

    For two or more branches, ``intersection.transverse`` uses the standard
    smooth-branch criterion: all participating branches must be smooth at the
    point and the tangent space of their scheme-theoretic intersection must
    have dimension ``ambient_dimension - sum(codimensions)``.  If the summed
    codimension exceeds the ambient dimension, a common point cannot be a
    transverse intersection.
    """
    eqs = _equations(equations)
    vars_ = normalize_problem_variables(variables, sp.Tuple(*eqs))
    point_map = normalize_point(point, vars_, context=(sp.Tuple(*eqs),))
    reduced = certified_radicalization(eqs, vars_, max_pieces=max_pieces)
    if not _point_satisfies(reduced.equations, point_map):
        raise ValueError("point does not lie on the algebraic variety")

    components = irreducible_components(reduced.equations, vars_, max_pieces=max_pieces)
    indices = tuple(
        i
        for i, component in enumerate(components)
        if _point_satisfies(component.equations, point_map)
    )
    branches = []
    for index in indices:
        component = components[index]
        space = tangent_space(component.equations, point_map, vars_)
        cone = tangent_cone(component.equations, point_map, vars_)
        degree = _tangent_cone_degree(cone)
        branches.append(
            BranchTangentGeometry(
                index,
                component,
                cone,
                space,
                space.dimension - component.dimension,
                degree,
                degree,
            )
        )

    intersection = None
    if len(indices) >= 2:
        intersection_equations = _component_intersection_equations(components, indices)
        intersection_dimension = polynomial_locus_dimension(intersection_equations, vars_)
        intersection_space = tangent_space(intersection_equations, point_map, vars_)
        intersection_cone = tangent_cone(intersection_equations, point_map, vars_)
        intersection_degree = _tangent_cone_degree(intersection_cone)
        codimension_sum = sum(len(vars_) - components[index].dimension for index in indices)
        expected_dimension = len(vars_) - codimension_sum
        all_smooth = all(branch.tangent_dimension_excess == 0 for branch in branches)
        transverse = (
            all_smooth
            and expected_dimension >= 0
            and intersection_space.dimension == expected_dimension
        )
        intersection = ComponentIntersectionGeometry(
            indices,
            intersection_equations,
            intersection_dimension,
            intersection_cone,
            intersection_space,
            intersection_space.dimension - intersection_dimension,
            transverse,
            intersection_degree,
            intersection_degree,
        )

    local_dimension = max(components[index].dimension for index in indices)
    singular = len(indices) >= 2 or any(branch.tangent_dimension_excess > 0 for branch in branches)
    return LocalBranchGeometry(
        vars_, point_map, indices, tuple(branches), intersection, local_dimension, singular
    )


@dataclass(frozen=True)
class TangentSpaceResult:
    variables: tuple[sp.Symbol, ...]
    point: Mapping[sp.Symbol, sp.Expr]
    jacobian: sp.Matrix
    basis: tuple[sp.Matrix, ...]
    dimension: int
    equations: tuple[sp.Expr, ...]


def tangent_space(equations, point, variables=None) -> TangentSpaceResult:
    """Return the Zariski tangent space at a point as the Jacobian nullspace."""

    eqs = _equations(equations)
    vars_ = (
        normalize_problem_variables(None, sp.Tuple(*eqs))
        if variables is None
        else normalize_symbol_sequence(variables)
    )
    point_map = normalize_point(point, vars_, context=(sp.Tuple(*eqs),))
    for eq in eqs:
        if certified_zero(eq.subs(point_map)) is not True:
            raise ValueError("point does not lie on the algebraic variety")
    jac = sp.Matrix([[sp.diff(eq, v).subs(point_map) for v in vars_] for eq in eqs])
    basis = tuple(jac.nullspace())
    directions = collision_free_real_symbols(
        "d", len(vars_), sp.Tuple(*eqs), tuple(vars_), tuple(point_map)
    )
    linear = tuple(
        sp.expand(sum(jac[i, j] * directions[j] for j in range(len(vars_))))
        for i in range(jac.rows)
    )
    return TangentSpaceResult(vars_, point_map, jac, basis, len(basis), linear)


@dataclass(frozen=True)
class TangentConeResult:
    variables: tuple[sp.Symbol, ...]
    point: Mapping[sp.Symbol, sp.Expr]
    direction_variables: tuple[sp.Symbol, ...]
    initial_forms: tuple[sp.Expr, ...]
    formula: sp.Expr
    certified: bool
    method: str

    @property
    def ideal_generators(self) -> tuple[sp.Expr, ...]:
        """Reduced Groebner-basis generators of the exact tangent-cone ideal."""
        return self.initial_forms


def _lowest_total_degree(poly: sp.Expr, variables: Sequence[sp.Symbol]) -> int:
    p = sp.Poly(sp.expand(poly), *variables)
    nonzero = [(monom, coeff) for monom, coeff in p.terms() if coeff != 0]
    if not nonzero:
        raise ValueError("zero polynomial has no order of vanishing")
    return min(sum(monom) for monom, _ in nonzero)


def _tangent_deformation(
    poly: sp.Expr,
    variables: Sequence[sp.Symbol],
    t: sp.Symbol,
) -> sp.Expr:
    """Return ``t**(-ord(poly)) * poly(t*x)`` as an exact polynomial."""

    order = _lowest_total_degree(poly, variables)
    scaled = sp.expand(poly.subs({v: t * v for v in variables}))
    quotient = sp.cancel(scaled / t**order)
    result = sp.expand(quotient)
    # This division is exact by construction.  Check it instead of allowing a
    # rational expression to leak into the Groebner computation.
    sp.Poly(result, t, *variables)
    return result


def _groebner_exact(polys: Sequence[sp.Expr], variables: Sequence[sp.Symbol]):
    """Compute an exact Groebner basis, allowing algebraic coefficients."""

    try:
        return sp.groebner(tuple(polys), *variables, order="lex", extension=True)
    except (sp.polys.polyerrors.CoercionFailed, sp.polys.polyerrors.GeneratorsError):
        # ``EX`` remains exact for symbolic coefficient expressions when SymPy
        # cannot place all exact coordinates in one algebraic extension.
        return sp.groebner(tuple(polys), *variables, order="lex", domain=sp.EX)


def _saturate_by_parameter(
    generators: Sequence[sp.Expr],
    parameter: sp.Symbol,
    variables: Sequence[sp.Symbol],
) -> tuple[sp.Expr, ...]:
    """Return generators of ``<generators> : parameter**infinity``.

    The standard Rabinowitsch elimination identity is used::

        I : t^infinity = (I + <1 - u*t>) intersect k[t, x].

    With lexicographic order and ``u`` first, the Groebner-basis elements not
    containing ``u`` generate the elimination ideal exactly.
    """

    u = fresh_dummy("tangent_cone_saturation")
    all_variables = (u, parameter, *variables)
    basis = _groebner_exact((*generators, 1 - u * parameter), all_variables)
    eliminated = tuple(
        sp.expand(poly.as_expr()) for poly in basis.polys if u not in poly.as_expr().free_symbols
    )
    if not eliminated:
        return (sp.Integer(0),)
    return eliminated


def _exact_tangent_cone_ideal(
    translated_generators: Sequence[sp.Expr],
    directions: Sequence[sp.Symbol],
) -> tuple[sp.Expr, ...]:
    """Compute ``in_m(I)`` by saturated m-adic deformation."""

    expanded = tuple(sp.expand(f) for f in translated_generators)
    nonzero = tuple(f for f in expanded if f != 0)
    if not nonzero:
        return ()

    t = fresh_dummy("tangent_cone_parameter")
    deformation = tuple(_tangent_deformation(f, directions, t) for f in nonzero)
    saturated = _saturate_by_parameter(deformation, t, directions)

    # The special fibre t=0 is the associated-graded (tangent-cone) ideal.
    fibre = tuple(sp.expand(g.subs(t, 0)) for g in saturated)
    fibre = tuple(g for g in fibre if g != 0)
    if not fibre:
        return ()

    basis = _groebner_exact(fibre, tuple(directions))
    if basis.polys == [sp.Poly(1, *directions)]:
        return (sp.Integer(1),)
    return tuple(sp.expand(poly.as_expr()) for poly in basis.polys)


def tangent_cone(equations, point, variables=None) -> TangentConeResult:
    """Return the exact ideal-theoretic Zariski tangent cone at ``point``.

    If ``I`` is the translated defining ideal and ``m`` is the maximal ideal at
    the origin, the tangent-cone ideal is the associated-graded initial ideal
    ``in_m(I)``.  It is computed exactly by the flat ``m``-adic deformation

    ``t**(-ord(f)) * f(t*d)``

    followed by saturation with respect to ``t`` and specialization at
    ``t = 0``.  Saturation is essential: it captures initial forms arising from
    cancellations between the supplied generators, so the result depends on
    the ideal rather than on a particular generating set.

    ``initial_forms`` holds a
    reduced Groebner basis for the exact tangent-cone ideal.  ``certified`` is
    always ``True``.
    """

    eqs = _equations(equations)
    vars_ = (
        normalize_problem_variables(None, sp.Tuple(*eqs))
        if variables is None
        else normalize_symbol_sequence(variables)
    )
    point_map = normalize_point(point, vars_, context=(sp.Tuple(*eqs),))
    for eq in eqs:
        if certified_zero(eq.subs(point_map)) is not True:
            raise ValueError("point does not lie on the algebraic variety")

    dirs = collision_free_real_symbols(
        "d", len(vars_), sp.Tuple(*eqs), tuple(vars_), tuple(point_map)
    )
    shift = {v: point_map[v] + d for v, d in zip(vars_, dirs, strict=True)}
    translated = tuple(sp.expand(eq.subs(shift)) for eq in eqs)
    cone_ideal = _exact_tangent_cone_ideal(translated, dirs)
    formula = sp.And(*(sp.Eq(f, 0) for f in cone_ideal)) if cone_ideal else sp.true
    return TangentConeResult(
        vars_,
        point_map,
        dirs,
        cone_ideal,
        formula,
        True,
        "saturated_m_adic_deformation",
    )


def is_singular(equations, point, variables=None, *, codimension: int | None = None) -> bool:
    """Return whether ``point`` is singular on the polynomial variety."""
    eqs = _equations(equations)
    vars_ = (
        normalize_problem_variables(None, sp.Tuple(*eqs))
        if variables is None
        else normalize_symbol_sequence(variables)
    )
    point_map = normalize_point(point, vars_, context=(sp.Tuple(*eqs),))
    for eq in eqs:
        if certified_zero(eq.subs(point_map)) is not True:
            raise ValueError("point does not lie on the algebraic variety")
    locus = singular_locus(eqs, vars_, codimension=codimension)
    value = sp.simplify(locus.subs(point_map))
    if value in (sp.true, True):
        return True
    if value in (sp.false, False):
        return False
    raise ValueError("singularity test did not reduce to an exact Boolean value")


def is_smooth(equations, variables=None, *, codimension: int | None = None) -> bool:
    """Return whether the real polynomial variety has empty singular locus."""
    from .decision import is_satisfiable

    eqs = _equations(equations)
    vars_ = normalize_problem_variables(variables, sp.Tuple(*eqs))
    locus = singular_locus(eqs, vars_, codimension=codimension)
    return not is_satisfiable(locus, vars_)


def tangent_dimension(equations, point, variables=None) -> int:
    """Return the exact Zariski tangent-space dimension at ``point``."""
    return tangent_space(equations, point, variables).dimension


__all__ = [
    "PolynomialMapImplicitizationResult",
    "ZariskiClosureResult",
    "implicitize_polynomial_map",
    "zariski_closure",
    "ReducedAlgebraicVariety",
    "IrreducibleAlgebraicComponent",
    "ReducedComponentSingularLocus",
    "local_branch_geometry",
    "LocalBranchGeometry",
    "BranchTangentGeometry",
    "ComponentIntersectionGeometry",
    "stratified_singular_geometry",
    "local_dimension_strata",
    "minimal_prime_intersections",
    "StratifiedSingularGeometry",
    "SingularGeometryStratum",
    "LocalDimensionStratum",
    "MinimalPrimeIntersection",
    "certified_radicalization",
    "irreducible_components",
    "reduced_component_singular_loci",
    "singular_locus",
    "TangentSpaceResult",
    "tangent_space",
    "TangentConeResult",
    "tangent_cone",
    "is_singular",
    "is_smooth",
    "tangent_dimension",
    "AlgebraicLocusPiece",
    "DecompositionPieceCertificate",
    "DecompositionCertificate",
    "EquidimensionalDecomposition",
    "RadicalComponent",
    "RadicalMinimalPrimeDecomposition",
    "RadicalIdealResult",
    "RegularChainComponent",
    "RegularChainDecomposition",
    "RegularChainSplitCertificate",
    "RegularChainDecompositionCertificate",
    "RadicalIdealCertificate",
    "MinimalPrimeDecompositionCertificate",
    "TriangularPrimalityStage",
    "TriangularPrimalityCertificate",
    "recursive_regular_chain_decomposition",
    "certify_triangular_primality",
    "radical_ideal",
    "certified_radical_minimal_prime_decomposition",
    "equidimensional_decomposition",
    "verify_regular_chain_decomposition_certificate",
    "verify_radical_ideal_certificate",
    "verify_triangular_primality_certificate",
    "verify_minimal_prime_decomposition_certificate",
    "PrimaryComponent",
    "PrimaryDecompositionCertificate",
    "PrimaryDecompositionResult",
    "AssociatedPrimesResult",
    "primary_decomposition",
    "associated_primes",
    "verify_primary_decomposition_certificate",
    "verify_decomposition_certificate",
    "HilbertData",
    "hilbert_function",
    "hilbert_polynomial",
    "ideal_degree",
    "ideal_hilbert_data",
    "IndependentLocalizationCertificate",
    "IndependentLocalizationResult",
    "LocalizationContractionCertificate",
    "LocalizationContractionResult",
    "SaturationStabilizationCertificate",
    "SaturationStabilizationResult",
    "certified_independent_localization",
    "contract_localized_ideal",
    "saturation_stabilization",
    "verify_independent_localization_certificate",
    "verify_localization_contraction_certificate",
    "verify_saturation_stabilization_certificate",
    "ZeroDimensionalPrimaryCertificate",
    "ZeroDimensionalPrimaryComponent",
    "ZeroDimensionalPrimaryResult",
    "zero_dimensional_primary_decomposition",
    "verify_zero_dimensional_primary_certificate",
    "FractionFieldZeroDimensionalCertificate",
    "FractionFieldZeroDimensionalResult",
    "GTZContractedComponentCertificate",
    "GTZNodeCertificate",
    "GTZPrimaryDecompositionCertificate",
    "GTZExecutionPlan",
    "GTZPrimaryDecompositionResult",
    "LocalPrimaryFactor",
    "clear_gtz_caches",
    "gtz_cache_info",
    "gtz_primary_decomposition",
    "plan_gtz_primary_decomposition",
    "verify_fraction_field_zero_dimensional_certificate",
    "verify_gtz_primary_decomposition_certificate",
]
