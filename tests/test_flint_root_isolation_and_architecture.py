from __future__ import annotations

from pathlib import Path

import sympy as sp

from semialg.algebraic.roots import isolate_real_roots


def test_clustered_newton_descartes_isolation_matches_exact_real_roots():
    x = sp.symbols("x")
    poly = sp.Poly(
        (x**2 - 2) * ((10000 * x - 14143) ** 2 - 2) * ((10000 * x - 14144) ** 2 - 2),
        x,
        domain=sp.QQ,
    )
    roots = isolate_real_roots(poly)
    expected = sorted(set(sp.real_roots(poly.as_expr())), key=lambda value: float(sp.N(value, 30)))
    assert len(roots) == len(expected)
    for root, exact in zip(roots, expected, strict=True):
        assert root.interval.left <= exact <= root.interval.right
        assert poly.count_roots(root.interval.left, root.interval.right) >= 1


def test_flint_is_internal_to_high_level_algebra_semantics():
    root = Path(__file__).resolve().parents[1] / "src" / "semialg"
    high_level = (
        root / "algebraic_decomposition.py",
        root / "algebraic" / "gtz.py",
        root / "algebraic" / "gtz_primary.py",
        root / "algebraic" / "gtz_zero_dim.py",
    )
    for path in high_level:
        text = path.read_text()
        assert "import flint" not in text
        assert "from flint" not in text
        assert "_flint_arithmetic" not in text


def test_flint_multivariate_squarefree_part_removes_multiplicity():
    from semialg._flint_arithmetic import flint_squarefree_part

    x, y = sp.symbols("x y")
    poly = sp.Poly((x + y) ** 3 * (x - y) ** 2, x, y, domain=sp.ZZ)
    result = flint_squarefree_part(poly)

    assert result is not None
    assert sp.expand(result.as_expr() - (x + y) * (x - y)) == 0
