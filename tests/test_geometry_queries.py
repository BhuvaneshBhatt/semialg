import sympy as sp

from semialg import (
    bounding_box,
    critical_values,
    distance_between_regions,
    distance_to_region,
    euler_characteristic,
    fiber,
    is_convex,
    is_path_connected,
    path_between,
    semialgebraic_image,
    semialgebraic_preimage,
    semialgebraic_projection,
)


def test_semialgebraic_projection_eliminates_existential_variable():
    x, y = sp.symbols("x y", real=True)
    result = semialgebraic_projection(sp.And(x >= 0, x <= y, y <= 2), [x], [x, y])
    assert sp.simplify_logic(sp.Xor(result, sp.And(y >= 0, y <= 2))) is sp.false


def test_semialgebraic_image_identity_interval():
    x, u = sp.symbols("x u", real=True)
    result = semialgebraic_image(x, sp.And(x >= -1, x <= 2), [x], image_variables=[u])
    assert sp.simplify_logic(sp.Xor(result, sp.And(u >= -1, u <= 2))) is sp.false


def test_semialgebraic_preimage_and_fiber():
    x, y, u, a = sp.symbols("x y u a", real=True)
    pre = semialgebraic_preimage(x**2 + y**2, u <= 1, [x, y], target_variables=[u])
    assert pre == (x**2 + y**2 <= 1)
    assert fiber(sp.And(x >= 0, x <= a), {a: 2}) == sp.And(x >= 0, x <= 2)


def test_bounding_box_exact():
    x, y = sp.symbols("x y", real=True)
    result = bounding_box(sp.And(x >= -1, x <= 2, y >= 0, y <= 3), [x, y])
    assert result == {x: (-1, 2), y: (0, 3)}


def test_distances_are_exact():
    x = sp.symbols("x", real=True)
    assert distance_to_region((2,), sp.And(x >= 0, x <= 1), [x]) == 1
    assert distance_between_regions(x <= 0, x >= 2, [x]) == 2


def test_critical_values_include_stationary_and_boundary_values():
    x = sp.symbols("x", real=True)
    assert critical_values(x**2, sp.And(x >= -1, x <= 2), [x]) == (0, 1, 4)


def test_is_convex_fast_paths():
    x, y = sp.symbols("x y", real=True)
    assert is_convex(x**2 + y**2 <= 1, [x, y]) is True
    disconnected = sp.Or(sp.And(x >= 0, x <= 1), sp.And(x >= 2, x <= 3))
    assert is_convex(disconnected, [x]) is False


def test_path_connectivity_and_path_chain_1d():
    x = sp.symbols("x", real=True)
    region = sp.And(x >= 0, x <= 1)
    assert is_path_connected(region, [x]) is True
    path = path_between(region, (0,), (1,), [x])
    assert path.connected is True
    assert path.explicit is True
    assert path.waypoints[0][x] == 0
    assert path.waypoints[-1][x] == 1


def test_compactly_supported_euler_characteristic_cells():
    x = sp.symbols("x", real=True)
    assert euler_characteristic(sp.And(x >= 0, x <= 1), [x]) == 1
    assert euler_characteristic(sp.And(x > 0, x < 1), [x]) == -1
    assert euler_characteristic(sp.Or(sp.Eq(x, 0), sp.Eq(x, 1)), [x]) == 2


def test_critical_values_include_constant_positive_dimensional_kkt_component():
    x, y = sp.symbols("x y", real=True)
    # grad(x**2*y) vanishes on the entire line x=0, where the objective is
    # constantly zero. The objective is otherwise unbounded above and below,
    # so zero cannot enter through attained global extrema.
    assert critical_values(x**2 * y, sp.true, [x, y]) == (0,)
