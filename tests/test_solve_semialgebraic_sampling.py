import pytest
import sympy as sp

from semialg import solve_semialgebraic

pytestmark = pytest.mark.slow


def test_per_component_samples_for_one_dimensional_components():
    x = sp.Symbol("x", real=True)

    sol = solve_semialgebraic([x**2 > 1], [x], count=0, samples="per_component")

    assert sol.satisfiable is True
    assert len(sol.components) == 2
    assert sol.samples == ({x: sp.Integer(-2)}, {x: sp.Integer(2)})
    for sample in sol.samples:
        assert bool((x**2 > 1).subs(sample)) is True


def test_output_samples_can_return_per_component_samples():
    x = sp.Symbol("x", real=True)

    samples = solve_semialgebraic(
        [x**2 <= 1, sp.Ne(x, 0)], [x], count=0, samples="per_component", output="samples"
    )

    assert samples == ({x: sp.Rational(-1, 2)}, {x: sp.Rational(1, 2)})


def test_default_count_uses_component_samples_when_available():
    x = sp.Symbol("x", real=True)

    sol = solve_semialgebraic([x**2 > 1], [x], count=1)

    assert sol.samples == ({x: sp.Integer(-2)},)
    assert sol.diagnostics["capabilities"]["structural_sampling"] is True
    assert sol.diagnostics["selected_sample_mode"] == "auto"


def test_per_cell_samples_for_supported_2d_cells():
    x, y = sp.symbols("x y", real=True)
    left_square = sp.And(x >= 0, x <= 1, y >= 0, y <= 1)
    right_square = sp.And(x >= 2, x <= 3, y >= 0, y <= 1)

    sol = solve_semialgebraic(sp.Or(left_square, right_square), [x, y], count=0, samples="per_cell")

    assert len(sol.cells) == 2
    assert sol.samples == (
        {x: sp.Rational(1, 2), y: sp.Rational(1, 2)},
        {x: sp.Rational(5, 2), y: sp.Rational(1, 2)},
    )


def test_sample_mode_alias_for_cells():
    x, y = sp.symbols("x y", real=True)

    samples = solve_semialgebraic(
        [x >= 0, x <= 1, y >= x, y <= 1],
        [x, y],
        count=0,
        sample_mode="cells",
        output="samples",
    )

    assert samples == ({x: sp.Rational(1, 2), y: sp.Rational(3, 4)},)
