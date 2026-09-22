# Multidimensional root-API testing

A public function is not considered deeply tested merely because several examples call it. `semialg` distinguishes **contract count** from **semantic dimensions**: independent tests should give an API qualitatively different ways to fail.

The generated registry `tests/root_api_adequacy_dimensions.toml` records 79 depth-tracked APIs and the additional dimension used to deepen each one. `scripts/generate_root_api_adequacy_dimensions.py` regenerates the registry from executable tests, and `test_root_api_multidimensional_adequacy_contract.py` verifies that every registry entry points to a test that directly calls the API.

| Dimension | Purpose | Typical example |
|---|---|---|
| `nominal` | Establish the ordinary public contract | known exact result |
| `metamorphic` | Detect inconsistencies without duplicating the implementation | translation, scaling, renaming, idempotence, or equivalent representation |
| `boundary-degenerate` | Exercise regime changes and singular cases | repeated roots, endpoints, rank drops, active constraints |
| `independent-oracle` | Compare with an independently known mathematical answer | exact vertices, dimension, distance, or rank |
| `independent-semantic` | Exercise a distinct semantic obligation | decomposition limits, topology, component structure |
| `round-trip-presentation` | Ensure presentation helpers do not alter mathematical data | discretize/plot stability |

The registry is an adequacy floor, not a claim of exhaustive correctness. A new regression should normally strengthen the dimension that exposed it, and a new algorithmic backend should preferentially gain differential or metamorphic evidence rather than another nominal example.

## Generated exact cases

Property-based tests are most useful when the oracle is known by construction. The deepening suite therefore uses generated exact integer/rational cases where appropriate—for example, zero-dimensional systems remain zero-dimensional after translations and nonzero scalar multiplication, and seeded samples must remain members of their exact canonical region. This avoids treating floating-point agreement as a proof.

## Algebraic certification boundary

The deepening pass found a zero-dimensional radicalization defect: the regular-chain squarefreeness check had treated every zero-dimensional branch as automatically squarefree. A shifted repeated root such as `(x - 2)**2` could therefore survive `certified_radicalization` unchanged. Zero-dimensional branches now pass the same exact regular-chain/separant certification as other branches. The regression is retained as a metamorphic multiplicity test.
