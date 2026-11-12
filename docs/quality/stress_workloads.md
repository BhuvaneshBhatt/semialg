# Large exact workloads

Use these tests to exercise algorithms near their practical size limits and to compare changes in algebraic preprocessing, quotient arithmetic, certificate replay and CAD. Each test checks an independently constructed mathematical answer. Runtime is reported, not asserted: making a case faster should keep it passing.

The regular suite includes the existing computational and performance regression cases. The workloads on this page carry the separate `stress` marker and require explicit selection:

```sh
python -m pytest -m stress tests/stress --durations=0
```

Run one case in a fresh process when measuring cold performance:

```sh
python -m pytest -m stress tests/stress/test_exact_workloads.py::test_dense_grid_rur --durations=0
```

For a bounded profiling run, use the diagnostic runner. Its profile includes instrumentation overhead, so compare timings with the same profiling mode:

```sh
python tools/profile_tests.py --output build/stress-profile --timeout 3600 -- -m stress tests/stress --durations=0
```

## Workloads and exact answers

| Case | Input | Expected answer | Main work |
|---|---|---|---|
| Dense grid RUR | Four dense quartics obtained from a Cartesian grid by an invertible linear map | 256 distinct rational points and a degree-256 squarefree defining polynomial | Native F4, quotient arithmetic, separator search and exact RUR verification |
| Multiple grid primary decomposition | Three dense sextics with root multiplicities 1, 2 and 3 | 27 primary components with total quotient length 216 | Primary splitting, nilpotents and exact certificate replay |
| Curved grid topology | A degree-48 inequality formed from 16 vertical and 16 parabolic boundaries inside a sheared box | 145 connected components | CAD lifting, adjacency and exact connectivity |
| Exceptional parameter fibers | Existence of x in [-1,1] with a T_64(x)=b, where T_64 is a Chebyshev polynomial | b² ≤ a², including the zero-polynomial fiber a=b=0 | Native parameter preprocessing, specialization guards and complete real elimination |

### Dense grid RUR

Let L_i=x_i+Σ_j x_j. The matrix of this change of coordinates is I+11ᵀ, whose inverse maps grid coordinates r to x_i=r_i−Σ_j r_j/(n+1). Impose ∏_{r=0}^{d−1}(L_i−r)=0 for each i. There are exactly dⁿ distinct rational points.

The test constructs those points independently, evaluates the returned separating form at each point and builds the expected defining polynomial as the product of the resulting linear factors. It also checks every coordinate equation modulo the defining polynomial and verifies that the separating form reconstructs the parameter.

### Multiple grid primary decomposition

Replace each factor L_i−r with (L_i−r)^(r+1). The support still has dⁿ points, but the local length at a grid point is ∏_i(r_i+1), so the total quotient length is [d(d+1)/2]ⁿ.

The test checks every component's support and degree, the complete multiplicity multiset, irredundancy and exact certificate replay. Counting components alone would miss incorrect multiplicities.

### Curved grid topology

Put v=y−x² and p(t)=∏_{r=1}^{2m}(t−r). Inside 0<x,v<2m+1, require p(x)p(v)>0. The map (x,v)↦(x,v+x²) is a homeomorphism. There are m+1 positive intervals and m negative intervals in either coordinate, so the region has (m+1)²+m² connected components.

The test checks the component count and a one-to-one correspondence between returned components and independently chosen rational witnesses in the grid rectangles. With m=8, it expects 145 components.

### Exceptional parameter fibers

The polynomial T_64 maps [-1,1] onto [-1,1]. Therefore a T_64(x)=b has a solution in that interval exactly when b²≤a². At a=0, the original equality becomes b=0; replacing the whole problem with its generic rational-function-field fiber would lose this case.

The test checks equivalence of the complete returned parameter formula with b²≤a² and checks the zero and nonzero constant fibers explicitly. It uses the complete solvability API; a partial root-count classification is insufficient for this contract.

## Resource handling

The dense algebra cases use explicit native reconstruction budgets, with primes starting at 1,000,003, up to 256 images and an 8,192-bit modulus. Larger primes reduce accidental zero coefficients in modular supports. Exact result verification remains enabled.

These workloads can require substantial time and memory. An external timeout reports an unfinished run; it does not establish mathematical correctness or count as a pass. Keep dependency versions, process state, native build and resource limits fixed when comparing results. Use the smallest representative size while investigating a failure, then rerun the full workload.
