import sympy as sp

from semialg import Box, critical_value_image
from semialg.topology.semialgebraic import (
    simplicial_betti_numbers,
    triangulate_region,
    triangulation_betti_numbers,
)


def test_topology_api_pipeline_is_compositional():
    box = Box(((0, 1), (0, 1)))
    triangulation = triangulate_region(box)
    assert simplicial_betti_numbers(triangulation.complex) == (1, 0, 0)
    assert triangulation_betti_numbers(box) == (1, 0, 0)


def test_critical_value_image_public_pipeline_handles_discrete_case():
    x = sp.Symbol("x", real=True)
    result = critical_value_image(x**2, sp.true, (x,))
    assert result.isolated_values == (0,)
