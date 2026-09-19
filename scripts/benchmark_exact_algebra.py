#!/usr/bin/env python3
"""Deterministic exact-algebra performance benchmark."""

from __future__ import annotations

import json
from time import perf_counter

import sympy as sp

from semialg.algebraic import (
    modular_groebner_basis_fraction_field,
    verify_modular_fraction_field_groebner_certificate,
)
from semialg.algebraic_function_fields import (
    MonogenicFunctionField,
    RationalFunctionField,
    compress_primitive_element,
    verify_primitive_element_compression,
)


def timed(call):
    start = perf_counter()
    value = call()
    return value, perf_counter() - start


def main():
    u, x, y = sp.symbols("u x y")
    domain = sp.QQ.frac_field(u)
    systems = (
        ("quadratic", (x**2 - u, y - x)),
        ("denominators", (x / (u + 1) + y, y**2 - u / (u + 2))),
        ("mixed", ((u + 1) * x + y, x * y - u)),
    )
    groebner_rows = []
    for name, generators in systems:
        modular, modular_seconds = timed(
            lambda generators=generators: modular_groebner_basis_fraction_field(
                generators, (x, y), (u,), max_primes=8
            )
        )
        direct, direct_seconds = timed(
            lambda generators=generators: sp.groebner(
                generators, x, y, order="grevlex", domain=domain
            )
        )
        certified = modular is not None and verify_modular_fraction_field_groebner_certificate(
            modular.certificate
        )
        same = modular is not None and modular.basis == tuple(
            sp.cancel(poly.as_expr()) for poly in direct.polys
        )
        groebner_rows.append(
            {
                "case": name,
                "modular_seconds": modular_seconds,
                "direct_seconds": direct_seconds,
                "certified": certified,
                "same_basis": same,
                "primes": list(modular.certificate.primes) if modular else [],
            }
        )

    a, b = sp.symbols("a b")
    base = RationalFunctionField(tuple())
    lower = MonogenicFunctionField(base, a, (-2, 0, 1))
    top = MonogenicFunctionField(lower, b, (lower.convert(-3), lower.zero, lower.one))
    compression, compression_seconds = timed(lambda: compress_primitive_element(top))

    payload = {
        "fraction_field_groebner": groebner_rows,
        "primitive_compression": {
            "seconds": compression_seconds,
            "degree": compression.field.degree,
            "steps": len(compression.steps),
            "certified": verify_primitive_element_compression(compression),
        },
    }
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
