# Certified algebraic roots

semialg uses exact algebraic roots throughout CAD lifting, sign determination, finite algebraic solving, and root classification. The advanced interval utilities expose the same certification discipline to callers that need explicit rational isolating intervals or separators.

This guide is about **univariate polynomial roots over `QQ`**. It does not replace the existing transcendental solver; see [Transcendental and algebraic solving scope](transcendental_scope.md).

## Isolating and counting roots in a rational interval

`certify_polynomial_root_interval()` counts **distinct** real roots of the square-free part of an exact rational polynomial in a rational interval. Multiplicity belongs to root objects and multiplicity APIs; interval topology counts distinct locations.

The fast exact path uses a Descartes variation count after the interval transformation. Variation 0 proves the interval is root-free and variation 1 proves a unique interior root. Larger variation counts are ambiguous, so semialg falls back to exact Sturm counting. This makes Sturm a certification fallback rather than the default cost on every interval.

<!-- semialg-exec -->
```python
import sympy as sp

from semialg.algebraic import (
    certify_polynomial_root_interval,
    verify_polynomial_root_interval_certificate,
)

x = sp.Symbol("x")
polynomial = (x**2 - 2) * (x - 3)
certificate = certify_polynomial_root_interval(polynomial, 0, 2, var=x)

assert certificate.root_count == 1
assert certificate.unique
assert verify_polynomial_root_interval_certificate(
    certificate,
    polynomial=polynomial,
    var=x,
    left=0,
    right=2,
)
```

The verifier can replay a certificate by itself. When the caller also supplies the expected polynomial, variable, interval, or endpoint convention, the verifier additionally binds the certificate to that external problem statement. This prevents accidentally accepting a valid certificate that belongs to a different root query.

## Open and closed endpoints

Endpoint roots are represented independently from interior root counts. `include_left` and `include_right` determine whether endpoint roots contribute to `root_count`.

<!-- semialg-exec -->
```python
import sympy as sp

from semialg.algebraic import certify_polynomial_root_interval

x = sp.Symbol("x")
polynomial = x * (x - 1)
closed = certify_polynomial_root_interval(polynomial, 0, 1, var=x)
open_interval = certify_polynomial_root_interval(
    polynomial,
    0,
    1,
    var=x,
    include_left=False,
    include_right=False,
)

assert closed.root_count == 2
assert open_interval.root_count == 0
assert closed.left_is_root and closed.right_is_root
```

## Rational separators

`rational_between_algebraic_reals(left, right)` returns an exact rational strictly between two ordered exact real algebraic values. It refines certified isolating intervals until their separation itself proves that the rational lies between the inputs.

CAD sector sampling uses the same separator, so standalone root utilities and cylindrical lifting do not maintain independent refinement policies.

<!-- semialg-exec -->
```python
import sympy as sp

from semialg.algebraic import isolate_real_roots, rational_between_algebraic_reals

x = sp.Symbol("x")
left, right = isolate_real_roots(sp.Poly(x**2 - 2, x, domain=sp.QQ))
separator = rational_between_algebraic_reals(left, right)

assert separator.is_Rational
assert sp.simplify(left.as_expr() < separator) is sp.true
assert sp.simplify(separator < right.as_expr()) is sp.true
```

## Sign-stable neighborhoods

`certified_sign_stable_root_neighborhood(root, companions)` refines the distinguished root interval until it contains exactly that root and no root of any companion polynomial. Every nonzero companion therefore has a fixed exact sign throughout the returned interval.

A companion that vanishes at the distinguished root is rejected: a nonzero sign cannot be stable on a neighborhood containing such a zero.

<!-- semialg-exec -->
```python
import sympy as sp

from semialg.algebraic import (
    certified_sign_stable_root_neighborhood,
    isolate_real_roots,
)

x = sp.Symbol("x")
positive_root = isolate_real_roots(sp.Poly(x**2 - 2, x, domain=sp.QQ))[1]
neighborhood = certified_sign_stable_root_neighborhood(
    positive_root,
    (x + 1, 3 - x),
)

assert neighborhood.signs == (1, 1)
assert neighborhood.interval.left < neighborhood.interval.right
```

## Exactness boundary

These interval APIs require exact univariate polynomials over `QQ`. They do not accept arbitrary analytic functions, numerical root brackets, or floating-point approximations as proof objects. Existing transcendental solving remains a separate subsystem with its own supported fragment and certification rules.

For signatures and result fields, see the [Algebraic reference](../reference/algebraic.md). For the general proof model, see [Exactness and certification](../concepts/exactness_and_certification.md) and [Certificate model](../concepts/certificate_model.md).
