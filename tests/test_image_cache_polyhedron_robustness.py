from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor

import pytest
import sympy as sp
from hypothesis import given
from hypothesis import strategies as st

from semialg import simplify_boole
from semialg.cache_utils import BoundedLRU
from semialg.context import ExactComputationContext, computation_context
from semialg.parametric_integration import integrate_over_parametric_region
from semialg.region_integrate import integrate_over_region
from semialg.standard_region_integrate import integrate_over_standard_region
from semialg.standard_regions import (
    BoxRegion,
    IntervalRegion,
    ParametricRegion,
    PolyhedronRegion,
    SimplexRegion,
    TetrahedronRegion,
    TransformedRegion,
)
from semialg.symbolic_regions import as_semialgebraic_region


def test_noninjective_transformed_interval_integrates_image_not_preimages() -> None:
    x, t = sp.symbols("x t", real=True)
    region = TransformedRegion(IntervalRegion(-1, 1), (t**2,), (t,))

    assert region.dimension() == 1
    formula = as_semialgebraic_region(region, (x,)).quantifier_free_formula()
    assert simplify_boole(sp.Xor(formula, sp.And(x >= 0, x <= 1)), (x,)) is sp.false
    assert integrate_over_standard_region(1, region, (x,)) == 1


def test_rank_reducing_transformed_box_uses_image_dimension() -> None:
    x, y, u, v = sp.symbols("x y u v", real=True)
    region = TransformedRegion(BoxRegion(((0, 1), (0, 1))), (u, 0), (u, v))

    assert region.dimension() == 1
    assert integrate_over_standard_region(1, region, (x, y)) == 1


def test_degenerate_simplex_uses_intrinsic_image_measure() -> None:
    x, y = sp.symbols("x y", real=True)
    region = SimplexRegion(((0, 0), (1, 0), (2, 0)))

    assert region.dimension() == 1
    assert integrate_over_standard_region(1, region, (x, y)) == 2


def test_redundant_parametric_fiber_is_not_silently_zero() -> None:
    x, u, v = sp.symbols("x u v", real=True)
    region = ParametricRegion((u, v), ((u, 0, 1), (v, 0, 1)), (u,))

    with pytest.raises(NotImplementedError):
        integrate_over_parametric_region(1, (x,), region)


def test_fixed_parametric_direction_is_removed_before_jacobian() -> None:
    x, y, u, v = sp.symbols("x y u v", real=True)
    region = ParametricRegion((u, v), ((u, 0, 1), (v, u, u)), (u, v))

    assert region.dimension() == 1
    assert integrate_over_parametric_region(1, (x, y), region) == sp.sqrt(2)


def test_intrinsic_dimension_uses_exact_real_geometry() -> None:
    x, y = sp.symbols("x y", real=True)
    assert (
        integrate_over_region(1, sp.Eq(x**2 + y**2, 0), (x, y), measure_dimension="intrinsic") == 1
    )


def test_polyhedron_allows_face_sharing_but_rejects_volume_overlap() -> None:
    x, y, z = sp.symbols("x y z", real=True)
    upper = TetrahedronRegion(((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)))
    lower = TetrahedronRegion(((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, -1)))
    nested = TetrahedronRegion(((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, sp.Rational(1, 2))))

    assert integrate_over_standard_region(
        1, PolyhedronRegion((upper, lower)), (x, y, z)
    ) == sp.Rational(1, 3)
    with pytest.raises(ValueError, match="disjoint interiors"):
        integrate_over_standard_region(1, PolyhedronRegion((upper, nested)), (x, y, z))


def test_cache_clear_invalidates_live_context_mirror() -> None:
    cache = BoundedLRU[int](4, "test.live-clear")
    with computation_context() as context:
        cache.put("key", 1)
        assert cache.get("key") == 1
        cache.clear()
        assert cache.get("key") is None
        assert context.stats()["test.live-clear"]["size"] == 0
        cache.put("key", 2)
        assert cache.get("key") == 2


def test_explicitly_shared_context_is_thread_safe() -> None:
    cache = BoundedLRU[int](16, "test.shared-context")
    context = ExactComputationContext()

    def worker(index: int) -> int:
        with computation_context(context):
            cache.put(index, index * index)
            value = cache.get(index)
            if value is None:
                raise AssertionError("shared cache lost a value")
            return value

    with ThreadPoolExecutor(max_workers=4) as pool:
        values = tuple(pool.map(worker, range(12)))
    assert values == tuple(index * index for index in range(12))
    assert context.cache_size("test.shared-context") == 12


@given(
    a=st.integers(-5, 5),
    b=st.integers(-5, 5),
    c=st.integers(-5, 5),
    d=st.integers(-5, 5),
)
def test_interval_boolean_membership_identities(a: int, b: int, c: int, d: int) -> None:
    left = IntervalRegion(min(a, b), max(a, b))
    right = IntervalRegion(min(c, d), max(c, d))
    x = sp.Symbol("x", real=True)
    left_formula = as_semialgebraic_region(left, (x,)).formula
    right_formula = as_semialgebraic_region(right, (x,)).formula

    identities = (
        sp.Equivalent(
            sp.Not(sp.Or(left_formula, right_formula)),
            sp.And(sp.Not(left_formula), sp.Not(right_formula)),
        ),
        sp.Equivalent(sp.And(left_formula, sp.Or(left_formula, right_formula)), left_formula),
        sp.Equivalent(
            sp.Xor(left_formula, right_formula),
            sp.Or(
                sp.And(left_formula, sp.Not(right_formula)),
                sp.And(right_formula, sp.Not(left_formula)),
            ),
        ),
    )
    for identity in identities:
        assert simplify_boole(sp.Not(identity), (x,)) is sp.false


def test_disjunctive_univariate_simplification_never_becomes_false() -> None:
    x = sp.Symbol("x", real=True)
    result = simplify_boole((x < 0) | (x > 1), (x,), semantic=False)
    assert result is not sp.false
    assert simplify_boole(sp.Xor(result, (x < 0) | (x > 1)), (x,)) is sp.false
