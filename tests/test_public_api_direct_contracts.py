import ast
import inspect
from pathlib import Path

import sympy as sp

import semialg
from semialg.formula import parse_quant_form_text

PRIMARY_CONTRACT_NAMES = (
    "PolytopeDecomposition",
    "MixedCellTetrahedralization",
    "MixedMeshTetrahedralization",
    "canonicalize_polygon",
    "canonicalize_polyhedron",
    "canonicalize_region",
    "triangulate_polytope",
    "decompose_polytope",
    "tetrahedralize_cell",
    "tetrahedralize_cells",
    "AffineBoxClip",
    "AffineHalfSpace",
    "AffineSpace",
    "Ball",
    "BooleanRegion",
    "Box",
    "CADRegion",
    "Capsule",
    "Circle",
    "Cone",
    "Cone",
    "CriticalValueImage",
    "Cube",
    "Cylinder",
    "Cylinder",
    "Dodecahedron",
    "Ellipsoid",
    "EllipsoidBoundary",
    "Exists",
    "FilledTorus",
    "ForAll",
    "Geometry",
    "HRepresentation",
    "HalfSpace",
    "Hexahedron",
    "Hyperplane",
    "Icosahedron",
    "Interval",
    "IrreducibleAlgebraicComponent",
    "Line",
    "LocalDimensionStratum",
    "MinimalPrimeIntersection",
    "Octahedron",
    "Parallelepiped",
    "Parallelepiped",
    "Parallelogram",
    "ParametricRegion",
    "Point",
    "FinitePointSet",
    "Polygon",
    "PolygonalComponent",
    "PolygonalSet",
    "PolyhedralComponent",
    "PolyhedralCone",
    "PolyhedralShell",
    "Polyhedron",
    "TetrahedralComplex",
    "PolynomialMapImplicitizationResult",
    "Polytope",
    "Prism",
    "Prism",
    "Pyramid",
    "Pyramid",
    "Ray",
    "ReducedAlgebraicVariety",
    "ReducedComponentSingularLocus",
    "RegionElement",
    "RegionNotElement",
    "RegularPolygon",
    "ResourceLimitError",
    "SemialgError",
    "SemialgStrategyFailure",
    "SemialgebraicRegion",
    "Simplex",
    "Simplex",
    "SingularGeometryStratum",
    "Sphere",
    "SphericalShell",
    "Stadium",
    "StandardRegion",
    "StratifiedSingularGeometry",
    "LocalBranchGeometry",
    "BranchTangentGeometry",
    "ComponentIntersectionGeometry",
    "Tetrahedron",
    "Simplex",
    "Torus",
    "TransformedRegion",
    "Triangle",
    "UnsupportedFragmentError",
    "ZariskiClosureResult",
    "affine_image",
    "affine_preimage",
    "affine_relative_interior_formula",
    "analyze_affine_map",
    "apply_quantifiers",
    "argmax_set",
    "argmin_set",
    "as_cad_region",
    "as_semialgebraic_region",
    "betti_number",
    "bounding_box",
    "cad",
    "centroid",
    "certified_radicalization",
    "classify_real_roots",
    "clip_affine_subspace_to_box",
    "closest_points",
    "closure_of_interior",
    "connected_component_count",
    "connected_component_samples",
    "connected_components",
    "contains_point",
    "convexity_certificate",
    "coordinate_range",
    "covariance_matrix",
    "critical_value_image",
    "critical_values",
    "deduplicate_indexed_vertices",
    "diameter",
    "discretize_region_geometry",
    "discretize_solution",
    "distance_between_regions",
    "distance_set",
    "distance_to_region",
    "equivalent",
    "euler_characteristic",
    "extrema_set",
    "fiber",
    "find_instance",
    "find_negative_point",
    "find_negative_witness_fast",
    "function_convex_partition",
    "function_convexity",
    "function_domain",
    "function_mapping_properties",
    "function_monotonic_partition",
    "function_monotonicity",
    "function_range",
    "function_sign",
    "function_sign_partition",
    "function_smoothness",
    "has_empty_interior",
    "implicitize_polynomial_map",
    "implies",
    "inertia_tensor",
    "inner_polygons",
    "inner_polyhedra",
    "integrate_over_region",
    "interior_of_closure",
    "intersects",
    "irreducible_components",
    "is_bijective",
    "is_bounded",
    "is_closed",
    "is_compact",
    "is_connected",
    "is_convex",
    "is_dense_in",
    "is_disjoint",
    "is_empty",
    "is_equal",
    "is_full_dimensional",
    "is_function_continuous",
    "is_function_smooth",
    "is_injective",
    "is_interior_disjoint",
    "is_open",
    "is_path_connected",
    "is_real_valued",
    "is_regular_closed_region",
    "is_regular_open_region",
    "is_satisfiable",
    "is_singular",
    "is_smooth",
    "is_subset",
    "is_surjective",
    "is_tautology",
    "is_zero_dimensional",
    "level_set",
    "linear_image",
    "local_dimension",
    "local_dimension_strata",
    "matrix_definiteness",
    "matrix_pd_on",
    "matrix_psd_on",
    "matrix_rank_on",
    "matrix_rank_stratification",
    "minimal_prime_intersections",
    "minkowski_sum",
    "moment_matrix",
    "nearest_point",
    "outer_polygons",
    "outer_polyhedra",
    "path_between",
    "plot_region_geometry",
    "plot_solution",
    "polygon_vertices",
    "polyhedron_face_indices",
    "polyhedron_vertices",
    "polynomial_nonnegative",
    "prove_negative",
    "prove_nonnegative",
    "prove_nonpositive",
    "prove_nonzero",
    "prove_positive",
    "prove_zero",
    "random_point",
    "random_points",
    "real_algebraic_feasibility",
    "reduce_formula",
    "reduce_region_integral",
    "reduced_component_singular_loci",
    "region_active_boundary_strata",
    "region_boundary",
    "region_boundary_result",
    "region_closure",
    "region_complement",
    "region_difference",
    "region_dimension",
    "region_image",
    "region_interior",
    "region_intersection",
    "region_measure",
    "region_moment",
    "region_nonsmooth_locus",
    "region_preimage",
    "region_product",
    "region_regular_locus",
    "region_singular_locus",
    "region_singular_locus_result",
    "region_symmetric_difference",
    "region_union",
    "region_variables",
    "replay_certificate",
    "resolve_formula",
    "sample_point",
    "sample_points",
    "scale",
    "semialgebraic_maximize",
    "semialgebraic_measure",
    "semialgebraic_minimize",
    "semialgebraic_projection",
    "sign_at",
    "sign_vector",
    "simplify_boole",
    "simplify_piecewise",
    "simplify_region",
    "simplify_system",
    "simplify_under_assumptions",
    "singular_locus",
    "solvability_conditions",
    "solve_real_algebraic_set",
    "solve_semialgebraic",
    "squared_distance_range",
    "stratified_singular_geometry",
    "local_branch_geometry",
    "polynomial_constraints",
    "active_constraints",
    "relative_interior",
    "relative_boundary",
    "semialgebraic_tangent_cone",
    "parameterization_geometry",
    "parameterization_critical_locus",
    "parameterization_critical_values",
    "parametric_cad",
    "implied_polynomial_inequality",
    "redundant_polynomial_inequalities",
    "nonnegative_combination_certificate",
    "verify_nonnegative_combination_certificate",
    "component_constraint_descriptions",
    "NonnegativeCombinationCertificate",
    "ImpliedPolynomialInequality",
    "RedundantPolynomialInequality",
    "ComponentConstraintDescription",
    "PolynomialConstraint",
    "PolynomialConstraintClause",
    "PolynomialConstraintSystem",
    "ActiveConstraintResult",
    "ParameterizationGeometry",
    "strict_feasible",
    "sublevel_set",
    "superlevel_set",
    "support_function",
    "tangent_cone",
    "tangent_dimension",
    "tangent_space",
    "topology_summary",
    "translate",
    "width",
    "zariski_closure",
    "zeng_negative_point",
    "Box",
    "FinitePointSet",
    "Interval",
    "Parallelogram",
    "TetrahedralComplex",
    "Zonotope",
    "RegionConversion",
    "convert_region",
    "polygonal_region_from_paths",
    "convex_hull",
    "random_polygon",
    "random_polytope",
    "subdivide_triangular_faces",
    "geodesic_refinement",
    "polyhedral_intersection",
    "polyhedral_boolean",
    "QuantifierEliminationResult",
    "quantifier_eliminate",
    "project_region",
)


def test_primary_api_identity_signature_and_docstrings():
    """Every root-level callable/type has an explicit import contract."""

    assert set(PRIMARY_CONTRACT_NAMES) == set(semialg.__all__) - {"__version__"}
    for name in PRIMARY_CONTRACT_NAMES:
        obj = getattr(semialg, name)
        assert callable(obj), name
        assert inspect.getdoc(obj), name
        if inspect.isfunction(obj):
            inspect.signature(obj)
        else:
            assert inspect.isclass(obj), name
            canonical_aliases = {
                "FinitePointSet": "FinitePointSet",
                "Interval": "Interval",
                "Box": "Box",
                "Parallelogram": "Parallelogram",
                "TetrahedralComplex": "TetrahedralComplex",
            }
            assert obj.__name__ == canonical_aliases.get(name, name)


def test_region_analysis_primary_function_contracts():
    x = sp.Symbol("x", real=True)
    region = semialg.SemialgebraicRegion(sp.And(x >= 0, x <= 1), (x,))

    assert semialg.local_dimension(region, {x: sp.Rational(1, 2)}) == 1
    assert semialg.local_dimension(region, {x: 2}) == -1
    assert semialg.region_variables(region) == (x,)
    singular_result = semialg.region_singular_locus_result(region)
    assert singular_result.complete
    assert singular_result.formula is sp.false
    assert semialg.region_singular_locus(region) is sp.false
    assert semialg.region_nonsmooth_locus(region) is sp.false
    assert semialg.region_active_boundary_strata(region)
    assert semialg.region_regular_locus(region) == sp.And(x >= 0, x <= 1)

    interior_of_closure = semialg.interior_of_closure(region)
    assert interior_of_closure.equals_region(sp.And(x > 0, x < 1))
    closure_of_interior = semialg.closure_of_interior(region)
    assert closure_of_interior.equals_region(sp.And(x >= 0, x <= 1))


def test_reduce_and_resolve_formula_root_contracts():
    satisfiable = parse_quant_form_text("exists x. x^2 = 1")
    impossible = parse_quant_form_text("exists x. x^2 = -1")
    tautology = parse_quant_form_text("forall x. x^2 >= 0")

    assert semialg.reduce_formula(satisfiable) is sp.true
    assert semialg.reduce_formula(impossible) is sp.false
    assert semialg.resolve_formula(satisfiable) is True
    assert semialg.resolve_formula(impossible) is False
    assert semialg.resolve_formula(tautology) is True


def test_plot_region_geometry_root_contract():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()
    try:
        plotted = semialg.plot_region_geometry(
            semialg.Box(((0, 1), (0, 1))),
            ax=ax,
        )
        assert plotted is ax
        assert ax.patches or ax.collections or ax.lines
    finally:
        plt.close(fig)


def test_every_public_function_is_called_by_the_test_suite():
    """Require a behavioral test reference for every exported function."""

    public = {name for name in semialg.__all__ if inspect.isfunction(getattr(semialg, name))}
    called: set[str] = set()
    test_dir = Path(__file__).parent
    for path in test_dir.rglob("test_*.py"):
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call):
                continue
            if isinstance(node.func, ast.Name) and node.func.id in public:
                called.add(node.func.id)
            elif (
                isinstance(node.func, ast.Attribute)
                and isinstance(node.func.value, ast.Name)
                and node.func.value.id == "semialg"
                and node.func.attr in public
            ):
                called.add(node.func.attr)
    missing = sorted(public - called)
    assert not missing, f"public functions without a direct test call: {missing}"


def test_indirect_public_type_contracts():
    """Public base/result/error types are covered through their production paths."""

    x = sp.Symbol("x", real=True)
    cad_region = semialg.as_cad_region((x >= 0) & (x <= 1), (x,))
    clip = semialg.clip_affine_subspace_to_box(
        (0, 0),
        ((1, 0),),
        ((-1, 1), (-1, 1)),
    )
    standard = semialg.Interval(0, 1)

    assert isinstance(cad_region, semialg.CADRegion)
    assert isinstance(clip, semialg.AffineBoxClip)
    assert isinstance(standard, semialg.StandardRegion)
    assert isinstance(standard, semialg.Geometry)
    assert issubclass(semialg.ResourceLimitError, semialg.SemialgError)
    assert issubclass(semialg.SemialgStrategyFailure, semialg.SemialgError)


def test_new_region_boolean_and_topological_relations_are_directly_callable():
    x = sp.Symbol("x", real=True)
    left = sp.And(x >= 0, x <= 1)
    right = sp.And(x >= 1, x <= 2)

    symmetric = semialg.region_symmetric_difference(left, right)
    assert semialg.contains_point(symmetric, (sp.Rational(1, 2),), (x,))
    assert not semialg.contains_point(symmetric, (1,), (x,))
    assert semialg.is_interior_disjoint(left, right, (x,))
