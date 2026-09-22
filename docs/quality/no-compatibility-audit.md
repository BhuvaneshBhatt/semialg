# No-backward-compatibility audit

Semialg 1.2 intentionally carries no deprecated root-name compatibility registry.
The package root resolves only canonical names from `PUBLIC_EXPORTS`.

The cleanup removed the dead `COMPAT_EXPORTS`, `COMPAT_PUBLIC_EXPORTS`, and
`COMPAT_REPLACEMENTS` machinery; the `RootFunction` alias; the
`compare_alg_numbers` / `sort_algebraic_numbers` aliases (the canonical names are
`compare_samples` / `sort_samples`); and `canonicalize_qe_formula` (the canonical
operation is `simplify_qe_formula`).  Equational-constraint use is represented by
`ProjectionConfig.use_ecs`; it is not a second projection backend name.

## Forwarding-function audit

An AST audit was run over root-exported functions whose executable body is a
single returned call.  The surviving forwarders were inspected as API operations,
not deleted mechanically.  They fall into intentional categories:

- **parameterized predicates:** the six `prove_*` functions select distinct sign
  relations while sharing `_prove_sign`;
- **result/projection pairs:** e.g. `region_singular_locus` versus
  `region_singular_locus_result`;
- **named matrix predicates:** `matrix_psd_on` and `matrix_pd_on` select distinct
  definiteness relations;
- **high-level geometry vocabulary:** `is_bounded`, `is_compact`, `is_closed`,
  `is_subset`, `is_equal`, and `is_disjoint` delegate to lower-level region
  reasoning while providing the primary query vocabulary;
- **representation adapters:** symbolic-region and CAD helpers convert between
  typed representations rather than preserving an obsolete name;
- **topology adapters:** `triangulation_betti_numbers` first triangulates a region
  and then applies the simplicial invariant, so it is not synonymous with
  `simplicial_betti_numbers`.

No additional forwarding function is a pure duplicate synonym with identical
public semantics and argument contract.  Thin implementation is
therefore not treated as evidence of legacy status by itself.
