"""Generate the multidimensional adequacy registry for root API functions."""

from __future__ import annotations

import ast
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEEP_FILES = (
    ROOT / "tests/test_root_api_multidimensional_algebraic_solving.py",
    ROOT / "tests/test_root_api_multidimensional_geometry.py",
)

TARGETS = (
    "Dodecahedron",
    "active_constraints",
    "affine_relative_interior_formula",
    "bounding_box",
    "canonicalize_polygon",
    "canonicalize_region",
    "certified_radicalization",
    "closest_points",
    "closure_of_interior",
    "component_constraint_descriptions",
    "critical_value_image",
    "decompose_polytope",
    "deduplicate_indexed_vertices",
    "discretize_region_geometry",
    "distance_between_regions",
    "distance_set",
    "euler_characteristic",
    "extrema_set",
    "fiber",
    "function_smoothness",
    "geodesic_refinement",
    "has_empty_interior",
    "implied_polynomial_inequality",
    "inertia_tensor",
    "inner_polygons",
    "inner_polyhedra",
    "interior_of_closure",
    "intersects",
    "irreducible_components",
    "is_closed",
    "is_dense_in",
    "is_full_dimensional",
    "is_open",
    "is_path_connected",
    "is_regular_closed_region",
    "is_regular_open_region",
    "is_singular",
    "is_smooth",
    "is_zero_dimensional",
    "level_set",
    "linear_image",
    "local_dimension_strata",
    "matrix_definiteness",
    "matrix_rank_stratification",
    "minimal_prime_intersections",
    "moment_matrix",
    "nearest_point",
    "outer_polygons",
    "outer_polyhedra",
    "parameterization_critical_locus",
    "parameterization_critical_values",
    "plot_region_geometry",
    "plot_solution",
    "polygon_vertices",
    "polyhedral_boolean",
    "polyhedral_intersection",
    "polyhedron_face_indices",
    "polyhedron_vertices",
    "polynomial_constraints",
    "random_point",
    "random_points",
    "random_polygon",
    "random_polytope",
    "reduced_component_singular_loci",
    "redundant_polynomial_inequalities",
    "region_boundary_result",
    "region_regular_locus",
    "region_singular_locus_result",
    "region_symmetric_difference",
    "region_variables",
    "relative_boundary",
    "relative_interior",
    "semialgebraic_tangent_cone",
    "squared_distance_range",
    "subdivide_triangular_faces",
    "sublevel_set",
    "superlevel_set",
    "tangent_dimension",
    "tetrahedralize_cells",
)


def minimum_rows():
    coverage = tomllib.loads((ROOT / "tests/public_api_coverage.toml").read_text())["apis"]
    return {
        name: {"owner": coverage[name]["owner"], "risk": coverage[name]["risk"]} for name in TARGETS
    }


def classify(test):
    n = test.lower()
    if any(
        k in n
        for k in (
            "invariant",
            "commute",
            "idempotent",
            "covariant",
            "agree",
            "reproducible",
            "same_optimum",
            "partition",
        )
    ):
        return "metamorphic"
    if any(
        k in n
        for k in (
            "boundary",
            "singular",
            "critical",
            "drop",
            "endpoint",
            "degenerate",
            "multiplicity",
        )
    ):
        return "boundary-degenerate"
    if any(
        k in n
        for k in (
            "exact",
            "contains",
            "membership",
            "realizes",
            "values",
            "vertices",
            "dimension",
            "rank",
        )
    ):
        return "independent-oracle"
    if any(k in n for k in ("presentation", "plot", "discret")):
        return "round-trip-presentation"
    return "independent-semantic"


def evidence():
    targets = set(minimum_rows())
    out = {n: [] for n in targets}
    for path in DEEP_FILES:
        tree = ast.parse(path.read_text())
        for fn in [
            n
            for n in tree.body
            if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef)) and n.name.startswith("test_")
        ]:
            calls = set()
            for node in ast.walk(fn):
                if not isinstance(node, ast.Call):
                    continue
                if isinstance(node.func, ast.Name):
                    calls.add(node.func.id)
                elif (
                    isinstance(node.func, ast.Attribute)
                    and isinstance(node.func.value, ast.Name)
                    and node.func.value.id == "sa"
                ):
                    calls.add(node.func.attr)
            for name in calls & targets:
                out[name].append(
                    (path.relative_to(ROOT).as_posix() + "::" + fn.name, classify(fn.name))
                )
    return out


def main():
    rows = minimum_rows()
    ev = evidence()
    missing = [n for n, v in ev.items() if not v]
    if missing:
        raise SystemExit(f"missing deepening evidence: {missing}")
    lines = [
        "# Generated by scripts/generate_root_api_adequacy_dimensions.py; do not edit manually.",
        "schema_version = 1",
        "",
    ]
    for name in sorted(rows):
        dimensions = ["nominal"] + sorted({d for _, d in ev[name]})
        lines += [
            f'[apis."{name}"]',
            f'owner = "{rows[name]["owner"]}"',
            f'risk = "{rows[name]["risk"]}"',
            "dimensions = [" + ", ".join(repr(x) for x in dimensions) + "]",
            "evidence = [" + ", ".join(repr(x) for x, _ in ev[name]) + "]",
            "",
        ]
    (ROOT / "tests/root_api_adequacy_dimensions.toml").write_text(
        "\n".join(lines).replace("'", '"') + "\n"
    )


if __name__ == "__main__":
    main()
