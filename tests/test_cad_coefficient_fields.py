"""Change coefficient fields and root separation without changing exact semantics."""

import pytest
import sympy as sp

from semialg.algebraic import isolate_real_roots
from semialg.algebraic.signs import sign_at_sample
from semialg.cad_algorithms import decomp_lazard_complete, verify_cad
from semialg.formula import parse_formula
from semialg.qe.complete import qe_from_cad


@pytest.mark.parametrize(
    "roots",
    [
        (sp.S.One, sp.S(2)),
        (-sp.sqrt(2), sp.sqrt(2)),
        (sp.sqrt(2) + sp.sqrt(3), sp.sqrt(2) + sp.sqrt(3) + sp.Rational(1, 1024)),
    ],
    ids=("rational", "quadratic", "compositum-cluster"),
)
@pytest.mark.parametrize("scale", [-1, 1])
def test_field_and_separation(roots, scale):
    x = sp.Symbol("x", real=True)
    polynomial = sp.Poly(scale * (x - roots[0]) ** 2 * (x - roots[1]) ** 3, x, extension=True)
    samples = isolate_real_roots(polynomial)
    assert len(samples) == 2
    assert [sample.multiplicity for sample in samples] == [2, 3]
    for sample, expected in zip(samples, roots, strict=True):
        assert sign_at_sample(sp.Poly(x - expected, x, extension=True), (sample,)) == 0
        assert sign_at_sample(polynomial, (sample,)) == 0
    assert samples[0].interval.right <= samples[1].interval.left


def test_compositum_qe():
    x = sp.Symbol("x", real=True)
    polynomial = (x - sp.sqrt(2)) ** 2 + (x - sp.sqrt(3)) ** 2
    source = decomp_lazard_complete((polynomial,), (x,))
    assert verify_cad(source, polynomials=(polynomial,))
    matrix = parse_formula(sp.Eq(polynomial, 0))
    assert qe_from_cad(source, (x,), (("exists", x),), matrix) is sp.false
    assert qe_from_cad(source, (x,), (("forall", x),), matrix) is sp.false
