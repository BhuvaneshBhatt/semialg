"""Full Lazard projection, native valuation and positive-dimensional curtains."""

from unittest.mock import patch

import pytest
import sympy as sp

from semialg.algebraic.roots import isolate_real_roots
from semialg.algebraic.samples import RationalSample
from semialg.cad_algorithms.decomposition import decomp_collins_complete, decomp_lazard_complete
from semialg.cad_algorithms.lifting.lazard import native_lazard_evaluate
from semialg.cad_algorithms.projection.lazard_complete import build_lazard_complete_projection
from semialg.formula import parse_formula
from semialg.qe.complete import qe_from_cad


def test_native_valuation_over_algebraic_parent():
    x, y, z = sp.symbols("x y z", real=True)
    alpha = isolate_real_roots(sp.Poly(x**2 - 2, x))[-1]
    polynomial = sp.Poly((x**2 - 2) ** 2 * (y - 1) ** 3 * (z**2 - y), x, y, z)
    with patch(
        "semialg.algebraic.samples.AlgebraicRoot.as_expr",
        side_effect=AssertionError("global reconstruction"),
    ):
        result = native_lazard_evaluate(polynomial, (x, y, z), (alpha, RationalSample(sp.S.One)))
    assert result.valuation == (2, 3)
    reduced = sp.groebner((x**2 - 2, y - 1), x, y, z).reduce(result.residue.as_expr())[1]
    assert sp.expand(reduced) == 8 * (z**2 - 1)


def test_full_operator_keeps_trailing_coefficients():
    x, y = sp.symbols("x y", real=True)
    tower = build_lazard_complete_projection((x * y**2 + y + x - 1,), (x, y))
    expressions = {p.as_expr() for p in tower.level(1).polynomials}
    assert x in expressions
    assert x - 1 in expressions
    assert any(entry.source == "trailing-coefficient" for entry in tower.level(1).entries)


@pytest.mark.parametrize("shift", [-1, 1])
def test_positive_dimensional_nullification(shift):
    x, y, t, z, w = sp.symbols("x y t z w", real=True)
    variables = (x, y, t, z, w)
    polynomial = w**2 + y * (z - shift) - x**2
    source = decomp_lazard_complete((polynomial,), variables)
    assert source.backend == "lazard-complete"
    assert source.complete
    lifts = source.tower.metadata["lazard_lifts"]
    assert any(any(v > 0 for v in item.valuation) for item in lifts)
    for value in (-2, 0, 2):
        cell = source._point_locator.locate((0, 0, value, shift, 0))
        assert cell.is_section
    result = qe_from_cad(source, variables, (("exists", w),), parse_formula(sp.Eq(polynomial, 0)))
    for xv in (-1, 0, 1):
        for yv in (-1, 0, 1):
            for zv in (shift - 1, shift, shift + 1):
                assert bool(result.subs({x: xv, y: yv, t: 2, z: zv})) is (
                    xv * xv - yv * (zv - shift) >= 0
                )


def test_ambient_backend_restarts_with_full_lazard():
    x, y, t, z, w = sp.symbols("x y t z w", real=True)
    result = decomp_collins_complete((w**2 + y * z - x**2,), (x, y, t, z, w))
    assert result.backend == "lazard-complete"
    assert result._point_locator.locate((0, 0, 4, -2, 0)).is_section


def test_native_valuation_over_nested_parent():
    from semialg.algebraic.samples import FiberRootContext
    from semialg.cad_algorithms.decomposition import _roots_in_parent_tower

    x, y, z = sp.symbols("x y z", real=True)
    alpha = isolate_real_roots(sp.Poly(x**2 - 2, x))[-1]
    context = FiberRootContext(sp.Poly(y**2 - x, x, y), (x, y), (alpha,), y)
    beta = _roots_in_parent_tower(context)[-1]
    polynomial = sp.Poly((x**2 - 2) ** 3 * (y**2 - x) ** 2 * (z**2 - y), x, y, z)
    with patch(
        "semialg.algebraic.samples.AlgebraicRoot.as_expr",
        side_effect=AssertionError("global reconstruction"),
    ):
        result = native_lazard_evaluate(polynomial, (x, y, z), (alpha, beta))
    assert result.valuation == (3, 2)
    reduced = sp.groebner((x**2 - 2, y**2 - x), y, x, z).reduce(result.residue.as_expr())[1]
    assert sp.expand(reduced) == 128 * (z**2 - y)


@pytest.mark.parametrize("ground", [sp.sqrt(2), sp.sqrt(2) + sp.sqrt(3)])
def test_algebraic_coefficient_field(ground):
    from semialg.cad_algorithms import verify_cad

    x, y = sp.symbols("x y", real=True)
    polynomial = y**2 - x - ground
    with patch(
        "semialg.algebraic.samples.AlgebraicRoot.as_expr",
        side_effect=AssertionError("global reconstruction"),
    ):
        result = decomp_lazard_complete((polynomial,), (x, y))
        replay = verify_cad(result, return_result=True)
    assert replay.valid, replay.failures
    matrix = parse_formula(sp.Eq(polynomial, 0))
    projected = qe_from_cad(result, (x, y), (("exists", y),), matrix)
    for value, expected in [(-ground - 1, False), (-ground, True), (1 - ground, True)]:
        assert bool(projected.subs(x, value)) is expected


def test_public_full_lazard_strategy():
    from semialg import cad

    x, y = sp.symbols("x y", real=True)
    result = cad(sp.Eq(y**2 - x, 0), (x, y), strategy="lazard-complete", return_result=True)
    assert result.diagnostics["planner"]["effective_backend"] == "lazard-complete"


def test_forced_restart_failure_propagates():
    from semialg.cache_utils import cache_disabled
    from semialg.errors import ExactComputationUnavailable

    x, y = sp.symbols("x y", real=True)
    with (
        cache_disabled(),
        patch(
            "semialg.cad_algorithms.decomposition.decomp_from_proj_tower",
            side_effect=ExactComputationUnavailable("forced lifting failure"),
        ),
    ):
        with pytest.raises(ExactComputationUnavailable, match="forced lifting failure"):
            decomp_collins_complete((y**2 - x,), (x, y))


@pytest.mark.parametrize("ground,expected", [(sp.sqrt(2), 0), (sp.sqrt(3), -1)])
def test_common_roots_across_ground_fields(ground, expected):
    from semialg.algebraic.comparison import compare_samples
    from semialg.algebraic.samples import FiberRootContext
    from semialg.algebraic.signs import sign_at_sample
    from semialg.cad_algorithms.decomposition import _roots_in_parent_tower

    x, y = sp.symbols("x y", real=True)
    prefix = (RationalSample(sp.S(2)),)
    left = FiberRootContext(sp.Poly(y * y - x, x, y), (x, y), prefix, y)
    right = FiberRootContext(sp.Poly(y - ground, x, y, extension=True), (x, y), prefix, y)
    with patch(
        "semialg.algebraic.samples.AlgebraicRoot.as_expr",
        side_effect=AssertionError("global reconstruction"),
    ):
        first = _roots_in_parent_tower(left)[-1]
        second = _roots_in_parent_tower(right)[0]
        assert compare_samples(first, second) == expected
        assert (
            sign_at_sample(sp.Poly(y - ground, x, y, extension=True), (*prefix, first)) == expected
        )


@pytest.mark.parametrize("backend", [decomp_collins_complete, decomp_lazard_complete])
def test_polynomial_stream(backend):
    x, y = sp.symbols("x y", real=True)
    polynomial = y * y - x
    source = backend((p for p in (polynomial,)), (x, y))
    projected = qe_from_cad(source, (x, y), (("exists", y),), parse_formula(sp.Eq(polynomial, 0)))
    for value in (-1, 0, 1):
        assert bool(projected.subs(x, value)) is (value >= 0)
