import sympy as sp

from semialg import Box, CriticalValueImage, betti_number, critical_value_image
from semialg.function_graph import semialgebraic_function_graph
from semialg.topology.semialgebraic import (
    simplicial_betti_numbers,
    triangulate_region,
    triangulation_betti_numbers,
)


def test_polyhedral_triangulation_drives_exact_higher_betti_numbers():
    box = Box(((0, 1), (0, 1), (0, 1)))
    triangulation = triangulate_region(box)
    assert simplicial_betti_numbers(triangulation.complex) == (1, 0, 0, 0)
    assert triangulation_betti_numbers(box) == (1, 0, 0, 0)
    assert betti_number(box, 2) == 0
    assert betti_number(box, 3) == 0


def test_real_wrapper_graphs_reuse_the_real_semialgebraic_argument_graph():
    x = sp.Symbol("x")
    t = sp.Symbol("t", real=True)
    re_graph = semialgebraic_function_graph(sp.re(x), t)
    im_graph = semialgebraic_function_graph(sp.im(x), t)
    re_subs = {x: 3, t: 3, re_graph.auxiliary_variables[0]: 3}
    im_subs = {x: 3, t: 0, im_graph.auxiliary_variables[0]: 3}
    assert sp.simplify(re_graph.formula.subs(re_subs)) is sp.true
    assert sp.simplify(im_graph.formula.subs(im_subs)) is sp.true


def test_critical_value_image_retains_positive_dimensional_singular_image():
    x, y = sp.symbols("x y", real=True)
    t = sp.Symbol("t", real=True)
    result = critical_value_image(x, sp.Eq(y**2, 0), (x, y), value_symbol=t)
    assert isinstance(result, CriticalValueImage)
    assert result.value_symbol == t
    assert result.component_images
    assert sp.simplify(result.formula) is sp.true


def test_cheap_factored_box_sign_certificate_handles_boundary_zeros():
    from semialg.proof import polynomial_sign_on_box

    z = sp.Symbol("z", real=True)
    assert polynomial_sign_on_box(z * (z - 1), {z: (sp.Integer(0), sp.Integer(1))}) == -1
    x = sp.Symbol("x", real=True)
    assert polynomial_sign_on_box(x * z + 1, {x: (0, 1), z: (0, 1)}) is None


def test_parameter_strata_use_certified_nonzero_representatives():
    from semialg.decomposition import parametric_cad

    x, p = sp.symbols("x p", real=True)
    result = parametric_cad(
        sp.And(p > 1, sp.Eq(x - p, 0)), (x,), parameters=(p,), specialize_fibers=False
    )
    assert result.strata
    assert all(
        bool(sp.simplify(stratum.condition.subs(stratum.sample))) for stratum in result.strata
    )
    assert all(stratum.sample[p] != 0 for stratum in result.strata if stratum.has_solution)


def test_intrinsic_integration_sums_disjoint_regular_cad_charts():
    from semialg import integrate_over_region

    x, y = sp.symbols("x y", real=True)
    region = sp.Or(
        sp.And(sp.Eq(y, 0), x >= 0, x <= 1),
        sp.And(sp.Eq(y, 0), x >= 2, x <= 3),
    )
    assert integrate_over_region(1, region, (x, y), measure_dimension="intrinsic") == 2
