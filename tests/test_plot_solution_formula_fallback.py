from __future__ import annotations

from types import SimpleNamespace

import sympy as sp


def test_plot_solution_falls_back_to_sampled_formula_region_for_2d_formula():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from semialg import plot_solution

    x, y = sp.symbols("x y", real=True)
    lens = (x**2 + y**2 <= 1) & ((x - sp.Rational(1, 2)) ** 2 + y**2 <= 1)
    solution = SimpleNamespace(variables=(x, y), formula=lens, cells=(), samples=())

    fig, ax = plt.subplots()
    try:
        plotted = plot_solution(
            solution,
            bounds=[
                (-sp.Rational(5, 4), sp.Rational(7, 4)),
                (-sp.Rational(5, 4), sp.Rational(5, 4)),
            ],
            raster_resolution=48,
            ax=ax,
        )
        assert plotted is ax
        assert ax.collections
    finally:
        plt.close(fig)
