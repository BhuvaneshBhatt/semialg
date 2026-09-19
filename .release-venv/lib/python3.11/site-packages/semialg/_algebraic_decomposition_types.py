"""Data model and deterministic certificate payloads for algebraic decomposition."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from hashlib import sha256

import sympy as sp


@dataclass(frozen=True)
class AlgebraicLocusPiece:
    """One exact polynomial-locus piece with a certified global dimension.

    ``dimension`` is the exact real semialgebraic dimension;
    ``algebraic_dimension`` records the Krull dimension of the polynomial
    ideal before real-set reduction. ``equidimensional`` is true only when the
    current implementation can prove that the represented reduced real set has
    no component of a different dimension. Multiplicity in the input ideal has
    no geometric significance.
    """

    variables: tuple[sp.Symbol, ...]
    equations: tuple[sp.Expr, ...]
    dimension: int
    algebraic_dimension: int
    equidimensional: bool

    @property
    def codimension(self) -> int:
        """Return the codimension in the recorded ambient polynomial ring."""
        return len(self.variables) - self.dimension


@dataclass(frozen=True)
class RegularChainSplitCertificate:
    """One replayable branch in a recursive regular-chain proof."""

    parent: tuple[sp.Expr, ...]
    children: tuple[tuple[sp.Expr, ...], ...]
    method: str
    splitter: sp.Expr | None = None


@dataclass(frozen=True)
class RegularChainDecompositionCertificate:
    """Search-independent proof payload for a recursive regular-chain union."""

    variables: tuple[sp.Symbol, ...]
    source_equations: tuple[sp.Expr, ...]
    splits: tuple[RegularChainSplitCertificate, ...]
    components: tuple[RegularChainComponent, ...]
    complete: bool
    reconstruction_complete: bool
    degree_complete: bool
    source_dimension: int
    source_degree: int
    reduced_degree: int | None


@dataclass(frozen=True)
class RegularChainComponent:
    """One recursively produced triangular component with exact invariants."""

    equations: tuple[sp.Expr, ...]
    chain: tuple[sp.Expr, ...]
    dimension: int
    degree: int
    squarefree_regular: bool


@dataclass(frozen=True)
class RegularChainDecomposition:
    """Exact union of recursively split triangular components.

    ``complete`` means every terminal ideal is certified to be the saturation
    of a squarefree regular chain. ``degree_complete`` is an independent
    Hilbert-degree check for equidimensional decompositions.
    """

    variables: tuple[sp.Symbol, ...]
    components: tuple[RegularChainComponent, ...]
    complete: bool
    reconstruction_complete: bool
    degree_complete: bool
    source_dimension: int
    source_degree: int
    reduced_degree: int | None
    methods: tuple[str, ...] = ()
    certificate: RegularChainDecompositionCertificate | None = None


@dataclass(frozen=True)
class RadicalIdealCertificate:
    """Replay payload proving generators equal the radical of the source ideal."""

    variables: tuple[sp.Symbol, ...]
    source_equations: tuple[sp.Expr, ...]
    generators: tuple[sp.Expr, ...]
    decomposition: RegularChainDecompositionCertificate


@dataclass(frozen=True)
class MinimalPrimeDecompositionCertificate:
    """Replay payload for radical reconstruction and minimal-prime claims."""

    variables: tuple[sp.Symbol, ...]
    source_equations: tuple[sp.Expr, ...]
    components: tuple[RadicalComponent, ...]
    primality_certificates: tuple[TriangularPrimalityCertificate | None, ...]
    radical_complete: bool
    minimal_primes_complete: bool
    degree_complete: bool


@dataclass(frozen=True)
class RadicalComponent:
    """One certified radical component in a reduced algebraic decomposition.

    ``prime`` is true only when semialg has an exact primality certificate; a
    radical component is never promoted to a minimal prime from triangular
    shape alone.
    """

    equations: tuple[sp.Expr, ...]
    dimension: int
    radical: bool
    prime: bool
    degree: int
    prime_method: str | None = None


@dataclass(frozen=True)
class RadicalMinimalPrimeDecomposition:
    """Certified decomposition of ``sqrt(I)`` into radical components.

    ``radical_complete`` certifies that the intersection of the returned
    radical component ideals is exactly ``sqrt(I)``.
    ``minimal_primes_complete`` additionally certifies that every component is
    prime and that the prime cover is irredundant; in that case the components
    are exactly the minimal primes of ``I``.
    """

    variables: tuple[sp.Symbol, ...]
    components: tuple[RadicalComponent, ...]
    radical_complete: bool
    minimal_primes_complete: bool
    degree_complete: bool
    certificate: MinimalPrimeDecompositionCertificate | None = None


@dataclass(frozen=True)
class RadicalIdealResult:
    """Exact generators for ``sqrt(I)`` when certified.

    ``complete`` is true only when ``generators`` have been independently
    reconstructed as the intersection of certified radical components and
    checked in both directions against the input radical.
    """

    variables: tuple[sp.Symbol, ...]
    generators: tuple[sp.Expr, ...]
    complete: bool
    method: str
    decomposition: RegularChainDecomposition
    certificate: RadicalIdealCertificate | None = None


@dataclass(frozen=True)
class DecompositionPieceCertificate:
    """Replay data for one reduced-real decomposition piece."""

    equations: tuple[sp.Expr, ...]
    dimension: int
    algebraic_dimension: int
    equidimensional: bool


@dataclass(frozen=True)
class DecompositionCertificate:
    """Replayable certificate for an equidimensional decomposition result.

    The certificate records the normalized input, canonical final pieces, the
    exact reduction/splitting mechanisms exercised, and a deterministic digest.
    :func:`verify_decomposition_certificate` replays the certified decomposition
    from the normalized input and checks that all recorded semantic claims match.
    """

    variables: tuple[sp.Symbol, ...]
    normalized_equations: tuple[sp.Expr, ...]
    pieces: tuple[DecompositionPieceCertificate, ...]
    complete: bool
    methods: tuple[str, ...]
    max_pieces: int
    digest: str


def _certificate_digest(
    variables: Sequence[sp.Symbol],
    normalized_equations: Sequence[sp.Expr],
    pieces: Sequence[DecompositionPieceCertificate],
    complete: bool,
    methods: Sequence[str],
    max_pieces: int,
) -> str:
    payload = repr(
        (
            tuple(map(str, variables)),
            tuple(sp.srepr(eq) for eq in normalized_equations),
            tuple(
                (
                    tuple(sp.srepr(eq) for eq in piece.equations),
                    piece.dimension,
                    piece.algebraic_dimension,
                    piece.equidimensional,
                )
                for piece in pieces
            ),
            bool(complete),
            tuple(methods),
            int(max_pieces),
        )
    ).encode("utf-8")
    return sha256(payload).hexdigest()


@dataclass(frozen=True)
class EquidimensionalDecomposition:
    """Exact finite cover of a reduced polynomial locus grouped by dimension.

    ``complete`` means every returned piece is certified equidimensional, so
    the result is suitable as the substrate for component-relative Jacobian
    regularity.  If it is false, the union is still exactly the same geometric
    zero set, but at least one piece may contain components of several local
    dimensions and callers must use a more complete method before reasoning
    componentwise.
    """

    variables: tuple[sp.Symbol, ...]
    pieces: tuple[AlgebraicLocusPiece, ...]
    complete: bool
    certificate: DecompositionCertificate | None = None

    @property
    def dimensions(self) -> tuple[int, ...]:
        return tuple(sorted({piece.dimension for piece in self.pieces}, reverse=True))

    def pieces_of_dimension(self, dimension: int) -> tuple[AlgebraicLocusPiece, ...]:
        """Return all pieces whose certified global dimension is ``dimension``."""
        return tuple(piece for piece in self.pieces if piece.dimension == dimension)


@dataclass(frozen=True)
class TriangularPrimalityStage:
    """One exact successive-extension step in a triangular primality proof."""

    leader: sp.Symbol
    degree: int
    factor_degrees: tuple[int, ...]
    reconstruction_verified: bool
    factorization_method: str


@dataclass(frozen=True)
class TriangularPrimalityCertificate:
    """Replay data for primality/reducibility of a squarefree regular chain."""

    variables: tuple[sp.Symbol, ...]
    chain: tuple[sp.Expr, ...]
    stages: tuple[TriangularPrimalityStage, ...]
    complete: bool
    prime: bool
    splitters: tuple[sp.Expr, ...] = ()
    method: str | None = None


@dataclass(frozen=True)
class _TriangularPrimalityAnalysis:
    """Internal exact analysis of a triangular chain over successive fields.

    ``prime_method`` is populated only when primality is proved. ``splitters``
    contains exact polynomial representatives of a certified factorization in
    the quotient by the preceding triangular prefix. An empty result means the
    current exact coefficient-domain backend cannot decide the next extension.
    """

    prime_method: str | None = None
    splitters: tuple[sp.Expr, ...] = ()


@dataclass(frozen=True)
class PrimaryComponent:
    """One exact primary component with its certified associated prime."""

    equations: tuple[sp.Expr, ...]
    radical: tuple[sp.Expr, ...]
    dimension: int
    degree: int
    primary: bool
    method: str


@dataclass(frozen=True)
class PrimaryDecompositionCertificate:
    """Replay payload for an exact primary decomposition."""

    variables: tuple[sp.Symbol, ...]
    source_equations: tuple[sp.Expr, ...]
    components: tuple[PrimaryComponent, ...]
    complete: bool
    irredundant: bool


@dataclass(frozen=True)
class PrimaryDecompositionResult:
    """Certified primary decomposition of a polynomial ideal."""

    variables: tuple[sp.Symbol, ...]
    components: tuple[PrimaryComponent, ...]
    complete: bool
    irredundant: bool
    method: str
    certificate: PrimaryDecompositionCertificate | None = None


@dataclass(frozen=True)
class AssociatedPrimesResult:
    """Certified associated primes obtained from a primary decomposition."""

    variables: tuple[sp.Symbol, ...]
    primes: tuple[tuple[sp.Expr, ...], ...]
    complete: bool
    method: str
    primary_decomposition: PrimaryDecompositionResult
