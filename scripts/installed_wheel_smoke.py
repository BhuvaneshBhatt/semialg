"""Smoke the built wheel from outside the source checkout.

This script intentionally uses only public APIs and examples mirrored in the
user documentation.  CI runs it with the wheel installed in a fresh virtual
environment and with the current working directory outside the repository.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import sympy as sp

import semialg
from semialg.algebraic_decomposition import associated_primes, primary_decomposition
from semialg.sos_certificates import SOSCertificate, verify_sos_certificate


def _assert_installed_import(source_root: Path | None) -> None:
    module_path = Path(semialg.__file__).resolve()
    if source_root is not None:
        source_package = (source_root / "src" / "semialg").resolve()
        if module_path.is_relative_to(source_package):
            raise AssertionError(f"semialg imported from source checkout: {module_path}")


def _public_api_smoke() -> None:
    x, y = sp.symbols("x y", real=True)
    assert semialg.is_satisfiable((x >= -2) & (x <= 3), [x])
    assert semialg.is_equal(x**2 <= 1, (x >= -1) & (x <= 1), [x])

    decomposition = primary_decomposition(((x + y) ** 2, y * (x + y)), (x, y))
    assert decomposition.complete
    assert semialg.replay_certificate(decomposition).verified is True

    positivity = semialg.polynomial_nonnegative((x**2 + 1), (x,), return_result=True)
    assert positivity.decision is True


def _documentation_examples() -> None:
    # Primary-decomposition tutorial: embedded associated prime.
    x, y = sp.symbols("x y")
    ideal = ((x + y) ** 2, y * (x + y))
    result = primary_decomposition(ideal, (x, y))
    assert result.complete and len(result.components) == 2
    primes = associated_primes(ideal, (x, y))
    assert primes.complete and len(primes.primes) == 2

    # SOS tutorial: exact Gram verification, independent of any numerical SDP backend.
    certificate = SOSCertificate(
        polynomial=x**2 + 1,
        variables=(x,),
        monomial_basis=(sp.Integer(1), x),
        gram_matrix=sp.ImmutableMatrix.eye(2),
    )
    assert verify_sos_certificate(x**2 + 1, certificate)

    # CAD tutorial: factored and expanded source formulas are semantically equivalent.
    factored = (x - 1) * (x + 1) >= 0
    expanded = x**2 - 1 >= 0
    assert semialg.is_equal(factored, expanded, [x])


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-root", type=Path)
    args = parser.parse_args(argv)
    _assert_installed_import(args.source_root)
    _public_api_smoke()
    _documentation_examples()
    print(f"installed semialg smoke passed from {Path(semialg.__file__).resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
