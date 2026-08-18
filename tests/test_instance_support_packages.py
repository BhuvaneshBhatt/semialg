import sympy as sp

from semialg.benchmarks import (
    gen_seeded_bench_wits,
    gen_seeded_points,
)
from semialg.instances import (
    cartesian_power_tuples,
    complete_partial_assign,
    extend_partial_cons,
    find_random_section_wit,
    sample_free_assignments,
    sample_modular_points,
)
from semialg.univariate import (
    find_univar_decomp_wit,
    intv_wits_for_form,
    witness_from_interval,
)
from semialg.validation import (
    crosscheck_backend_pred,
    form_sat_by_assign,
    verify_candidate_witness,
)


def test_instance_support_01():
    assert cartesian_power_tuples([0, 1], 2) == [(0, 0), (0, 1), (1, 0), (1, 1)]
    pts = sample_modular_points(2, 3, 4, seed=1)
    assert len(pts) == 4
    assert all(len(pt) == 2 for pt in pts)


def test_instance_support_02():
    x, y = sp.symbols("x y", integer=True)
    assigns = sample_free_assignments((x, y), sample_count=3, seed=2)
    assert len(assigns) == 3
    completed = complete_partial_assign({x: 2}, (x, y), sample_count=2, seed=3)
    assert all(x in a and y in a for a in completed)
    ext = extend_partial_cons({x: 1}, sp.Eq(x + y, 3), (x, y), sample_count=10, seed=4)
    assert any(a[y] == 2 for a in ext.extensions)


def test_instance_support_03():
    x = sp.Symbol("x", integer=True)
    assert form_sat_by_assign(sp.Eq(x, 2), {x: 2})
    verdict = verify_candidate_witness(sp.Eq(x, 2), {x: 2})
    assert verdict.is_valid


def test_instance_support_04():
    x = sp.Symbol("x", integer=True)
    result = crosscheck_backend_pred(sp.Eq(x, 1), "toy", lambda a: a[x] == 1, [{x: 1}, {x: 2}])
    assert result.consistent


def test_instance_support_05():
    x = sp.Symbol("x", real=True)
    witness = find_univar_decomp_wit(x**2 - 1, sp.true, x)
    assert witness.witness in {-1, 1, sp.Integer(-1), sp.Integer(1)} or sp.Abs(witness.witness) == 1
    assert witness_from_interval((0, 2)) == 1
    interval_candidates = intv_wits_for_form(sp.And(x > 0, x < 2), x, [(0, 2)])
    assert interval_candidates[0].witness == 1


def test_instance_support_06():
    x, y = sp.symbols("x y", real=True)
    batch = gen_seeded_bench_wits((x, y), seed=5, sample_count=3)
    assert batch.seed == 5
    assert len(batch.assignments) == 3
    modular = gen_seeded_points(2, 5, seed=6, sample_count=4)
    assert len(modular.assignments) == 4
    section = find_random_section_wit(
        sp.And(sp.Eq(x + y, 2), sp.Eq(x - y, 0)), (x, y), seed=7, attempts=4
    )
    assert section.attempts >= 0
