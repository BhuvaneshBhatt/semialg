"""Immutable result and certificate models for modular algebra algorithms."""

from __future__ import annotations

from dataclasses import dataclass

import sympy as sp


@dataclass(frozen=True)
class ModularFactorizationResult:
    """Exact factorization reconstructed from a modular integer algorithm."""

    unit: sp.Rational
    factors: tuple[tuple[tuple[sp.Rational, ...], int], ...]
    certified: bool
    method: str = "sympy_dup_zz_modular+exact_reconstruction"


@dataclass(frozen=True)
class ModularGroebnerCertificate:
    """Replayable evidence for a reconstructed Groebner basis over ``QQ``.

    ``membership_representations[i][j]`` is the exact multiplier of source
    generator ``j`` in reconstructed basis polynomial ``i``.  Together with
    exact reduction of every source generator by ``basis`` this proves ideal
    equality; an exact Buchberger computation on ``basis`` proves the
    Groebner property independently of the modular proposal.
    """

    variables: tuple[sp.Symbol, ...]
    order: str
    source_generators: tuple[sp.Expr, ...]
    basis: tuple[sp.Expr, ...]
    membership_representations: tuple[tuple[sp.Expr, ...], ...]
    primes: tuple[int, ...]
    modulus: int
    certified: bool = True
    method: str = "multimodular_crt+rational_reconstruction+exact_membership"


@dataclass(frozen=True)
class ModularGroebnerResult:
    """A certified rational Groebner basis reconstructed from finite fields."""

    basis: tuple[sp.Expr, ...]
    certificate: ModularGroebnerCertificate

    @property
    def certified(self) -> bool:
        return self.certificate.certified


@dataclass(frozen=True)
class ModularFractionFieldGroebnerCertificate:
    """Exact certificate for a modular Groebner basis over ``QQ(U)``."""

    variables: tuple[sp.Symbol, ...]
    parameters: tuple[sp.Symbol, ...]
    order: str
    source_generators: tuple[sp.Expr, ...]
    basis: tuple[sp.Expr, ...]
    membership_representations: tuple[tuple[sp.Expr, ...], ...]
    primes: tuple[int, ...]
    modulus: int
    certified: bool = True
    method: str = "modular_QQ(U)+crt+rational_function_reconstruction+exact_membership"


@dataclass(frozen=True)
class ModularFractionFieldGroebnerResult:
    basis: tuple[sp.Expr, ...]
    certificate: ModularFractionFieldGroebnerCertificate

    @property
    def certified(self) -> bool:
        return self.certificate.certified


@dataclass(frozen=True)
class ModularResultantCertificate:
    """Replayable certificate for a modularly reconstructed resultant.

    The exact verifier uses the fixed-degree Sylvester determinant on a
    deterministic tensor grid.  The recorded degree bounds guarantee that
    equality on that grid proves equality as a polynomial in the parameters.
    """

    variable: sp.Symbol
    parameters: tuple[sp.Symbol, ...]
    first: sp.Expr
    second: sp.Expr
    resultant: sp.Expr
    parameter_degree_bounds: tuple[int, ...]
    verification_grid: tuple[tuple[int, ...], ...]
    primes: tuple[int, ...]
    modulus: int
    certified: bool = True
    method: str = "multimodular_resultant+crt+rational_reconstruction+exact_sylvester_grid"


@dataclass(frozen=True)
class ModularResultantResult:
    """Certified resultant reconstructed from finite-field images."""

    resultant: sp.Expr
    certificate: ModularResultantCertificate

    @property
    def certified(self) -> bool:
        return self.certificate.certified


@dataclass(frozen=True)
class ModularSubresultantCertificate:
    """Replayable certificate for a reconstructed subresultant sequence.

    Finite fields only propose the sequence.  Replay recomputes the exact
    native subresultant normalization and compares every coefficient; this is
    stricter than checking a generic PRS recurrence, whose scale
    convention would not by itself identify SymPy's subresultants.
    """

    variable: sp.Symbol
    parameters: tuple[sp.Symbol, ...]
    first: sp.Expr
    second: sp.Expr
    polynomials: tuple[sp.Expr, ...]
    resultant_certificate: ModularResultantCertificate
    primes: tuple[int, ...]
    modulus: int
    certified: bool = True
    method: str = "multimodular_subresultants+crt+rational_reconstruction+exact_replay"


@dataclass(frozen=True)
class ModularSubresultantResult:
    """Certified subresultant sequence reconstructed from finite fields."""

    polynomials: tuple[sp.Expr, ...]
    resultant: sp.Expr
    certificate: ModularSubresultantCertificate

    @property
    def certified(self) -> bool:
        return self.certificate.certified
