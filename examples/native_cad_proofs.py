"""Native CAD proofs for a curtain and a complex-to-real root birth."""

import sympy as sp

from semialg.algebraic.samples import RationalSample
from semialg.cad_algorithms import decomp_lazard_complete, verify_cad
from semialg.formula import parse_formula
from semialg.qe.complete import qe_from_cad
from semialg.topology.branch_correspondence import (
    certify_branch_correspondence,
    replay_branch_correspondence,
)

x, y, t, z, w = sp.symbols("x y t z w", real=True)
variables = (x, y, t, z, w)
polynomial = w**2 + y * z - x**2
source = decomp_lazard_complete((polynomial,), variables)
proof = verify_cad(source, polynomials=(polynomial,), return_result=True)
assert proof.valid, proof.failures
assert proof.checked_valuation_claims > 0
projected = qe_from_cad(source, variables, (("exists", w),), parse_formula(sp.Eq(polynomial, 0)))
for xv in (-1, 0, 1):
    for yv in (-1, 0, 1):
        for zv in (-1, 0, 1):
            assert bool(projected.subs({x: xv, y: yv, t: 3, z: zv})) is (xv * xv - yv * zv >= 0)

limit = certify_branch_correspondence(
    ((y + 1) ** 2 + x) * (y - 2),
    y,
    0,
    (x,),
    (RationalSample(sp.S.One),),
    (sp.S.Zero,),
)
assert limit.value == 2
assert replay_branch_correspondence(limit)
print("Native CAD replay and the surviving branch are certified.")
