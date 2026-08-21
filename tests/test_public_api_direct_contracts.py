import ast
import inspect
from pathlib import Path

import sympy as sp

import semialg
from semialg.formula import parse_quant_form_text


def test_previously_uncovered_public_function_contracts():
    x, y = sp.symbols("x y", real=True)

    a = semialg.IntervalRegion(0, 2)
    b = semialg.IntervalRegion(1, 3)
    assert semialg.RegionIntersection(a, b).op == "intersection"
    assert semialg.RegionDifference(a, b).op == "difference"

    u = sp.Symbol("u", real=True)
    pregion = semialg.ParametricRegion((u,), ((u, 0, 1),), (u,))
    transformed, limits, jac = semialg.reduce_parametric_integral(x, [x], pregion)
    assert transformed == u and limits == ((u, 0, 1),) and jac == 1

    formula = sp.And(x >= 0, x <= 1, y >= 0, y <= 1)
    decomp = semialg.extract_structured_cad_cells(formula, [x, y])
    cyl = semialg.cylindrical_solution_from_structured(decomp)
    assert cyl.cells
    vertical = semialg.structured_cad_cells_to_vertical_bounds_2d(decomp.cells)
    assert vertical
    graph = semialg.extract_cad_connectivity(cyl)
    assert graph.cells

    comps = semialg.component_instances_text(
        "(x >= 0) & (x <= 1)", variables=["x"], return_result=True
    )
    assert comps is not None

    parsed = parse_quant_form_text("exists x. x^2 = 1")
    witness = semialg.find_instance_formula(parsed, return_result=True)
    assert witness is not None
    reduced = semialg.reduce_formula(parsed)
    assert reduced in (sp.true, True)
    resolved = semialg.resolve_formula(parsed)
    assert resolved is True

    ax = semialg.plot_region_geometry(semialg.BoxRegion(((0, 1), (0, 1))))
    assert ax is not None


def test_every_public_function_is_called_by_the_test_suite():
    """Require a behavioral test reference for every exported function."""

    public = {name for name in semialg.__all__ if inspect.isfunction(getattr(semialg, name))}
    called: set[str] = set()
    test_dir = Path(__file__).parent
    for path in test_dir.glob("test_*.py"):
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
