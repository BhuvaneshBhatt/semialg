from __future__ import annotations

import sympy as sp

from semialg import SemialgebraicRegion


def test_projection_membership_matches_existential_completion():
    x, y = sp.symbols("x y", real=True)
    region = SemialgebraicRegion(sp.And(x**2 + y**2 <= 1, y >= 0), (x, y))
    projected = region.project((y,))

    assert projected.variables == (x,)
    for value, expected in ((0, True), (1, True), (2, False)):
        assert projected.contains((value,)) is expected


def test_invertible_affine_image_then_preimage_recovers_region():
    x, y = sp.symbols("x y", real=True)
    source = SemialgebraicRegion(sp.And(x >= -1, x <= 2), (x,))
    image = source.image(3 * x - 4, variables=(y,))
    roundtrip = image.preimage(3 * x - 4, (x,))
    assert roundtrip.equals_region(source)


def test_composed_images_agree_with_composed_map():
    x, y, z = sp.symbols("x y z", real=True)
    source = SemialgebraicRegion(sp.And(x >= 0, x <= 1), (x,))
    two_step = source.image(x + 1, variables=(y,)).image(y**2, variables=(z,))
    direct = source.image((x + 1) ** 2, variables=(z,))
    assert two_step.equals_region(direct)
