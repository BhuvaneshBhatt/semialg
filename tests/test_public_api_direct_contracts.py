import ast
import inspect
from pathlib import Path

import sympy as sp

import semialg
from semialg.formula import parse_quant_form_text

PRIMARY_CONTRACT_NAMES = (
    "is_singular",
    "is_smooth",
    "singular_locus",
    "tangent_cone",
    "tangent_dimension",
    "tangent_space",
    "CADRegion",
    "as_cad_region",
    "replay_certificate",
    "analyze_affine_map",
    "convexity_certificate",
    "function_convexity",
    "function_convex_partition",
    "function_monotonicity",
    "function_monotonic_partition",
    "function_sign_partition",
    "function_smoothness",
    "function_mapping_properties",
    "is_injective",
    "is_surjective",
    "is_bijective",
    "is_function_continuous",
    "is_function_smooth",
    "matrix_definiteness",
    "matrix_pd_on",
    "matrix_psd_on",
    "matrix_rank_on",
    "matrix_rank_stratification",
    "is_convex",
    "equivalent",
    "implies",
    "is_satisfiable",
    "real_algebraic_feasibility",
    "solve_real_algebraic_set",
    "polynomial_nonnegative_decision",
    "find_negative_point",
    "polynomial_nonnegative",
    "zeng_negative_point",
    "find_negative_witness_fast",
    "is_tautology",
    "solve_semialgebraic",
    "cad",
    "affine_transform",
    "argmax_set",
    "argmin_set",
    "centroid",
    "closest_points",
    "connected_components",
    "contains_point",
    "RegionElement",
    "RegionNotElement",
    "coordinate_range",
    "covariance_matrix",
    "diameter",
    "distance_set",
    "extrema_set",
    "has_empty_interior",
    "inertia_tensor",
    "intersects",
    "is_bounded",
    "is_closed",
    "is_compact",
    "is_connected",
    "is_dense_in",
    "is_disjoint",
    "is_empty",
    "is_equal",
    "is_full_dimensional",
    "is_open",
    "is_subset",
    "level_set",
    "linear_image",
    "minkowski_sum",
    "moment_matrix",
    "nearest_point",
    "scale",
    "squared_distance_range",
    "sublevel_set",
    "superlevel_set",
    "support_function",
    "translate",
    "width",
    "function_domain",
    "is_real_valued",
    "ResourceLimitError",
    "SemialgError",
    "SemialgStrategyFailure",
    "UnsupportedFragmentError",
    "bounding_box",
    "critical_values",
    "critical_value_image",
    "CriticalValueImage",
    "distance_between_regions",
    "distance_to_region",
    "euler_characteristic",
    "fiber",
    "is_path_connected",
    "path_between",
    "semialgebraic_image",
    "semialgebraic_preimage",
    "semialgebraic_projection",
    "semialgebraic_measure",
    "region_measure",
    "region_centroid",
    "region_covariance",
    "region_moment",
    "function_range",
    "semialgebraic_maximize",
    "semialgebraic_minimize",
    "connected_component_count",
    "connected_component_samples",
    "topology_summary",
    "betti_number",
    "solvability_conditions",
    "Exists",
    "AffineBoxClip",
    "clip_affine_subspace_to_box",
    "ForAll",
    "apply_quantifiers",
    "prove_negative",
    "prove_nonnegative",
    "prove_nonpositive",
    "prove_positive",
    "prove_zero",
    "prove_nonzero",
    "function_sign",
    "simplify_system",
    "simplify_under_assumptions",
    "local_dimension",
    "region_boundary_result",
    "region_active_boundary_strata",
    "region_nonsmooth_locus",
    "region_regular_locus",
    "region_singular_locus",
    "region_singular_locus_result",
    "integrate_over_region",
    "reduce_region_integral",
    "region_boundary",
    "region_closure",
    "region_complement",
    "region_components",
    "region_difference",
    "region_dimension",
    "region_interior",
    "region_intersection",
    "region_product",
    "region_union",
    "classify_real_roots",
    "sample_point",
    "random_point",
    "random_points",
    "affine_relative_interior_formula",
    "strict_feasible",
    "sample_points",
    "sign_at",
    "sign_vector",
    "discretize_region_geometry",
    "discretize_solution",
    "plot_region_geometry",
    "plot_solution",
    "find_instance",
    "is_zero_dimensional",
    "reduce_formula",
    "resolve_formula",
    "affine_image",
    "affine_preimage",
    "region_image",
    "region_preimage",
    "HRepresentation",
    "Point",
    "AffineSpace",
    "Hyperplane",
    "HalfSpace",
    "AffineHalfSpace",
    "Line",
    "Ray",
    "Polytope",
    "Simplex",
    "Triangle",
    "Polygon",
    "RegularPolygon",
    "ConicRegion",
    "Parallelepiped",
    "Cube",
    "Tetrahedron",
    "Octahedron",
    "Icosahedron",
    "Dodecahedron",
    "Prism",
    "Pyramid",
    "Hexahedron",
    "Ball",
    "Sphere",
    "Circle",
    "Ellipsoid",
    "EllipsoidBoundary",
    "Cylinder",
    "Cone",
    "Torus",
    "FilledTorus",
    "BallRegion",
    "BooleanRegion",
    "BoxRegion",
    "CapsuleRegion",
    "ConeRegion",
    "CylinderRegion",
    "IntervalRegion",
    "ParallelepipedRegion",
    "ParallelogramRegion",
    "ParametricRegion",
    "PointRegion",
    "PolygonRegion",
    "PolyhedronRegion",
    "PrismRegion",
    "PyramidRegion",
    "RegionDifference",
    "RegionIntersection",
    "RegionSymmetricDifference",
    "RegionUnion",
    "SimplexRegion",
    "SphereRegion",
    "SphericalShellRegion",
    "StadiumRegion",
    "Geometry",
    "StandardRegion",
    "TetrahedronRegion",
    "TransformedRegion",
    "SemialgebraicRegion",
    "as_semialgebraic_region",
    "region_closure_interior",
    "region_interior_closure",
    "region_variables",
    "is_regular_closed_region",
    "is_regular_open_region",
    "simplify_region",
    "simplify_boole",
    "simplify_piecewise",
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
            assert obj.__name__ == name


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

    closure_interior = semialg.region_closure_interior(region)
    assert closure_interior.equals_region(sp.And(x > 0, x < 1))
    interior_closure = semialg.region_interior_closure(region)
    assert interior_closure.equals_region(sp.And(x >= 0, x <= 1))


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
            semialg.BoxRegion(((0, 1), (0, 1))),
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
    standard = semialg.IntervalRegion(0, 1)

    assert isinstance(cad_region, semialg.CADRegion)
    assert isinstance(clip, semialg.AffineBoxClip)
    assert isinstance(standard, semialg.StandardRegion)
    assert isinstance(standard, semialg.Geometry)
    assert issubclass(semialg.ResourceLimitError, semialg.SemialgError)
    assert issubclass(semialg.SemialgStrategyFailure, semialg.SemialgError)
