"""Local germs, sign strata, and parameter-dependent invariants."""

import sympy as sp

from semialg import local_germ, local_sign_strata, parameter_strata

x, a = sp.symbols("x a", real=True)

germ = local_germ(sp.Ne(x, 0), (0,), (x,))
assert germ.point[x] == 0

signs = local_sign_strata((x,), germ.formula, germ.point, germ.variables)
assert {piece.signs for piece in signs} == {(-1,), (1,)}

strata = parameter_strata(
    sp.Eq(x**2, a),
    (x,),
    (a,),
    conditions_by_value={
        "zero root": sp.Eq(x, 0),
        "nonzero root": sp.Ne(x, 0),
    },
)
assert tuple(piece.value for piece in strata) == ("zero root", "nonzero root")
assert sp.simplify(strata[0].condition.subs(a, 0)) is sp.true
assert sp.simplify(strata[1].condition.subs(a, 1)) is sp.true
