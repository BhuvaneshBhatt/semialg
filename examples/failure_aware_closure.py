"""Root births, degree loss and the limits of multibase continuation."""

import sympy as sp

from semialg.algebraic.samples import RationalSample
from semialg.errors import ExactComputationUnavailable
from semialg.topology.branch_correspondence import (
    certify_branch_correspondence,
    replay_branch_correspondence,
)

x, y, z = sp.symbols("x y z", real=True)
source = (RationalSample(sp.S.One),)
proof = certify_branch_correspondence(
    ((y + 1) ** 2 + x) * (y - 2), y, 0, (x,), source, (sp.S.Zero,)
)
assert proof.value == 2
assert replay_branch_correspondence(proof)

escaping = certify_branch_correspondence(x * y**2 + y, y, 0, (x,), source, (sp.S.Zero,))
finite = certify_branch_correspondence(x * y**2 + y, y, 1, (x,), source, (sp.S.Zero,))
assert escaping.value == -sp.oo
assert finite.value == 0
assert replay_branch_correspondence(escaping)

# Two base coordinates need a continuation proof on the whole source cell;
# choosing its sample-to-boundary segment alone would not establish that claim.
try:
    certify_branch_correspondence(
        x * z**2 + y * z + 1,
        z,
        0,
        (x, y),
        (RationalSample(sp.S.One), RationalSample(sp.S(3))),
        (sp.S.Zero, sp.S.Zero),
    )
except ExactComputationUnavailable:
    status = "unavailable"
else:
    raise AssertionError("unsupported multibase degree loss was accepted")
assert status == "unavailable"
