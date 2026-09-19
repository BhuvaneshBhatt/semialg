from __future__ import annotations

import sympy as sp

from semialg import (
    betti_number,
    connected_component_count,
    connected_component_samples,
    topology_summary,
)
from semialg.roadmaps import roadmap


def test_zero_dimensional_components_match_real_roots():
    x = sp.symbols("x", real=True)
    region = sp.Eq((x - 1) * (x + 2), 0)

    assert connected_component_count(region, (x,)) == 2
    samples = connected_component_samples(region, (x,))
    assert {sample[x] for sample in samples} == {-2, 1}
    assert betti_number(region, 0, (x,)) == 2


def test_circle_is_its_own_exact_roadmap():
    x, y = sp.symbols("x y", real=True)
    circle = sp.Eq(x**2 + y**2, 1)

    result = roadmap(circle, (x, y))

    assert result.dimension == 1
    assert result.rm1_certified is True
    assert result.rm2_certified is True
    assert result.roadmap_formula == circle
    assert result.component_count == 1


def test_compact_convex_set_has_projection_section_roadmap():
    x, y = sp.symbols("x y", real=True)
    disk = x**2 + y**2 <= 1

    result = roadmap(disk, (x, y))

    assert result.dimension == 1
    assert result.construction == "convex-projection-section"
    assert result.rm1_certified is True
    assert result.rm2_certified is True
    assert sp.simplify(result.roadmap_formula.subs({x: 0, y: 0})) is sp.true


def test_compact_curve_topology_from_components_and_euler():
    x, y = sp.symbols("x y", real=True)
    circle = sp.Eq(x**2 + y**2, 1)

    summary = topology_summary(circle, (x, y))

    assert summary.dimension == 1
    assert summary.connected_components == 1
    assert summary.euler_characteristic == 0
    assert summary.betti_numbers == (1, 1)
    assert betti_number(circle, 0, (x, y)) == 1
    assert betti_number(circle, 1, (x, y)) == 1


def test_compact_interval_has_no_cycle():
    x = sp.symbols("x", real=True)
    interval = sp.And(x >= 0, x <= 1)

    summary = topology_summary(interval, (x,))

    assert summary.euler_characteristic == 1
    assert summary.betti_numbers == (1, 0)


def test_compact_convex_disk_has_trivial_positive_betti_numbers():
    x, y = sp.symbols("x y", real=True)
    disk = x**2 + y**2 <= 1

    summary = topology_summary(disk, (x, y))

    assert summary.connected_components == 1
    assert summary.euler_characteristic == 1
    assert summary.betti_numbers == (1, 0, 0)
    assert betti_number(disk, 1, (x, y)) == 0
    assert betti_number(disk, 2, (x, y)) == 0
