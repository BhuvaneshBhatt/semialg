import sympy as sp

import semialg

t = sp.symbols("t", real=True)
u, v = sp.symbols("u v", real=True)
y0, y1 = sp.symbols("y0 y1", real=True)


def test_cusp_parameterization_critical_locus():
    result = semialg.parameterization_geometry(
        (t**2, t**3), sp.true, (t,), image_variables=(y0, y1)
    )
    assert result.generic_rank == 1
    assert result.domain_dimension == 1
    assert result.image_dimension == 1
    assert result.generic_fiber_dimension == 0
    assert semialg.equivalent(result.critical_locus, sp.Eq(t, 0), (t,))
    assert semialg.equivalent(
        result.critical_value_set, sp.And(sp.Eq(y0, 0), sp.Eq(y1, 0)), (y0, y1)
    )


def test_surface_parameterization_generic_rank():
    result = semialg.parameterization_geometry((u, v, u * v), sp.true, (u, v))
    assert result.generic_rank == 2
    assert result.generic_fiber_dimension == 0
    assert result.critical_locus is sp.false


def test_lower_dimensional_domain_uses_intrinsic_tangent_rank():
    result = semialg.parameterization_geometry((u,), sp.Eq(v, 0), (u, v))
    assert result.domain_dimension == 1
    assert result.image_dimension == 1
    assert result.generic_rank == 1
    assert result.critical_locus is sp.false


def test_simplex_parameter_domain_uses_equality_tangent_space():
    result = semialg.parameterization_geometry(
        (u,), sp.And(sp.Eq(u + v, 1), u >= 0, v >= 0), (u, v)
    )
    assert result.generic_rank == 1
    assert result.critical_locus is sp.false


def test_focused_critical_queries():
    locus = semialg.parameterization_critical_locus((t**2, t**3), sp.true, (t,))
    values = semialg.parameterization_critical_values(
        (t**2, t**3), sp.true, (t,), image_variables=(y0, y1)
    )
    assert semialg.equivalent(locus, sp.Eq(t, 0), (t,))
    assert semialg.equivalent(values, sp.And(sp.Eq(y0, 0), sp.Eq(y1, 0)), (y0, y1))
