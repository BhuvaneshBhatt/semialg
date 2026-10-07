"""Restricted, certifiable global multivariate mapping theorems."""

from __future__ import annotations

from dataclasses import dataclass

import sympy as sp

from .results import TruthValue


@dataclass(frozen=True)
class GlobalMapCertificate:
    classification: TruthValue
    kind: str
    jacobian: sp.Matrix
    source: str

    def verify(self):
        if self.source != "affine-full-rank":
            return False
        if any(any(v.free_symbols for v in row) for row in self.jacobian.tolist()):
            return False
        r = self.jacobian.rank()
        if self.kind == "injective":
            return (self.classification is TruthValue.TRUE) == (r == self.jacobian.cols)
        if self.kind == "surjective":
            return (self.classification is TruthValue.TRUE) == (r == self.jacobian.rows)
        if self.kind == "bijective":
            return (self.classification is TruthValue.TRUE) == (
                self.jacobian.rows == self.jacobian.cols == r
            )
        return False


def affine_global_mapping(expressions, variables, *, kind="bijective"):
    exprs = sp.Matrix(expressions)
    vars_ = tuple(map(sp.sympify, variables))
    J = exprs.jacobian(vars_)
    if any(any(v.free_symbols.intersection(vars_) for v in row) for row in J.tolist()):
        return GlobalMapCertificate(
            TruthValue.UNKNOWN, kind, J, "unsupported-nonlinear"
        )
    r = J.rank()
    ok = {
        "injective": r == len(vars_),
        "surjective": r == len(exprs),
        "bijective": len(exprs) == len(vars_) == r,
    }[kind]
    return GlobalMapCertificate(
        TruthValue.TRUE if ok else TruthValue.FALSE, kind, J, "affine-full-rank"
    )


@dataclass(frozen=True)
class ComposedMapCertificate:
    components: tuple[GlobalMapCertificate, ...]
    kind: str
    classification: TruthValue
    source: str = "composition-of-certified-maps"

    def verify(self):
        if self.source != "composition-of-certified-maps" or not self.components:
            return False
        if not all(c.verify() for c in self.components):
            return False
        if self.kind == "bijective":
            expected = all(
                c.kind == "bijective" and c.classification is TruthValue.TRUE
                for c in self.components
            )
        elif self.kind == "injective":
            expected = all(
                c.kind in ("injective", "bijective")
                and c.classification is TruthValue.TRUE
                for c in self.components
            )
        elif self.kind == "surjective":
            expected = all(
                c.kind in ("surjective", "bijective")
                and c.classification is TruthValue.TRUE
                for c in self.components
            )
        else:
            return False
        return (self.classification is TruthValue.TRUE) == expected


def compose_global_maps(certificates, *, kind="bijective"):
    certs = tuple(certificates)
    if not certs:
        return ComposedMapCertificate((), kind, TruthValue.UNKNOWN)
    if kind == "bijective":
        ok = all(
            c.kind == "bijective" and c.classification is TruthValue.TRUE and c.verify()
            for c in certs
        )
    elif kind == "injective":
        ok = all(
            c.kind in ("injective", "bijective")
            and c.classification is TruthValue.TRUE
            and c.verify()
            for c in certs
        )
    elif kind == "surjective":
        ok = all(
            c.kind in ("surjective", "bijective")
            and c.classification is TruthValue.TRUE
            and c.verify()
            for c in certs
        )
    else:
        raise ValueError("kind must be injective, surjective, or bijective")
    return ComposedMapCertificate(
        certs, kind, TruthValue.TRUE if ok else TruthValue.UNKNOWN
    )
