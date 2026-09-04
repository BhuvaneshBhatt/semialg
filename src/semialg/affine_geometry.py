"""Exact structural analysis of affine maps."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import sympy as sp

from .normalization import normalize_variables


@dataclass(frozen=True)
class AffineMapAnalysis:
    """Algebraic properties of an affine map ``x -> A*x + b``.

    ``rank`` is the symbolic generic matrix rank when coefficients contain free
    parameters. ``invertible`` is ``True`` or ``False`` only when the
    determinant is certified nonzero or zero; otherwise it is ``None``.
    """

    matrix: sp.ImmutableMatrix
    offset: sp.ImmutableMatrix
    variables: tuple[sp.Symbol, ...]
    rank: int
    invertible: bool | None
    determinant: sp.Expr | None
    scaled_isometry: bool
    scale_factor: sp.Expr | None
    inverse_matrix: sp.ImmutableMatrix | None
    inverse_offset: sp.ImmutableMatrix | None

    @property
    def isometry(self) -> bool:
        return self.scaled_isometry and sp.simplify(self.scale_factor - 1) == 0


def _zero_matrix(matrix: sp.MatrixBase) -> bool:
    return all(sp.simplify(entry) == 0 for entry in matrix)


def analyze_affine_map(
    mapping: Sequence[object] | sp.MatrixBase,
    variables: Sequence[sp.Symbol | str] | None = None,
    *,
    offset: Sequence[object] | None = None,
) -> AffineMapAnalysis:
    """Analyze an affine expression map or an explicit matrix/offset pair.

    With ``variables`` supplied, ``mapping`` is interpreted as a vector of
    affine expressions.  Without variables it is interpreted as a matrix and
    ``offset`` supplies the translation vector.
    """

    if variables is None:
        matrix = sp.Matrix(mapping)
        vars_ = tuple(sp.Symbol(f"x{i + 1}", real=True) for i in range(matrix.cols))
        off = (
            sp.zeros(matrix.rows, 1)
            if offset is None
            else sp.Matrix(tuple(map(sp.sympify, offset)))
        )
        if off.shape != (matrix.rows, 1):
            raise ValueError("offset length must equal the matrix row count")
    else:
        vars_ = tuple(normalize_variables(variables, append_context_symbols=False))
        exprs = tuple(map(sp.sympify, mapping))
        matrix = sp.Matrix(exprs).jacobian(vars_)
        if any(entry.free_symbols & set(vars_) for entry in matrix):
            raise ValueError("mapping is not affine in the supplied variables")
        zero_subs = {var: 0 for var in vars_}
        off = sp.Matrix([sp.simplify(expr.subs(zero_subs)) for expr in exprs])
        if tuple(sp.simplify(expr) for expr in matrix * sp.Matrix(vars_) + off) != tuple(
            sp.simplify(expr) for expr in exprs
        ):
            raise ValueError("mapping is not affine in the supplied variables")
        if offset is not None:
            raise ValueError("offset is only accepted when mapping is supplied as a matrix")

    rank = int(matrix.rank())
    determinant = sp.simplify(matrix.det()) if matrix.rows == matrix.cols else None
    invertible: bool | None = False if matrix.rows != matrix.cols else None
    if determinant is not None:
        if determinant.is_zero is True or determinant == 0:
            invertible = False
        elif determinant.is_zero is False or (determinant.is_number and determinant != 0):
            invertible = True
    gram = sp.simplify(matrix.T * matrix)
    scale_sq: sp.Expr | None = None
    scaled = False
    if matrix.cols > 0 and gram.rows == gram.cols:
        candidate = sp.simplify(gram[0, 0])
        if candidate.is_positive is True and _zero_matrix(gram - candidate * sp.eye(matrix.cols)):
            scaled = True
            scale_sq = candidate
    scale = sp.sqrt(scale_sq) if scale_sq is not None else None
    inv_matrix = None
    inv_offset = None
    if invertible:
        inv = sp.simplify(matrix.inv())
        inv_matrix = sp.ImmutableMatrix(inv)
        inv_offset = sp.ImmutableMatrix(sp.simplify(-inv * off))
    return AffineMapAnalysis(
        sp.ImmutableMatrix(matrix),
        sp.ImmutableMatrix(off),
        vars_,
        rank,
        invertible,
        determinant,
        scaled,
        scale,
        inv_matrix,
        inv_offset,
    )


__all__ = ["AffineMapAnalysis", "analyze_affine_map"]
