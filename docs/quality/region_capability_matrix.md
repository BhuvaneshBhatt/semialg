# Region-computation capability matrix

| Capability | Status | Exact/certified layer | Notes |
|---|---|---|---|
| Unified symbolic region | Implemented | Exact | `SemialgebraicRegion(formula, variables)` |
| Named standard shapes | Implemented | Exact | Convertible to unified regions |
| Membership expressions | Implemented | Exact | Symbolic lowering via `RegionElement` |
| Subset/equality/disjointness | Implemented | Exact | QE/CAD-backed where needed |
| Boolean operations | Implemented | Exact | Formula layer and reusable CAD layer |
| Projection | Implemented | Exact | Existential QE |
| Image | Implemented | Exact | Graph relation + elimination, with fast affine cases elsewhere |
| Preimage | Implemented | Exact | Substitution where applicable |
| Reusable CAD region | Implemented | Exact | `CADRegion` |
| CAD signatures | Implemented | Structural exact summary | Includes a stable digest of decomposition structure |
| Point location/sign vector | Implemented | Exact for resolved algebraic points | Reuses existing CAD |
| Dimension-stratified cell complex | Implemented | Exact | Codimension-one closure incidence |
| Connected components | Implemented | Exact/partial by topology path | CAD connectivity infrastructure |
| Euler characteristic | Implemented | Exact | Compact-support semialgebraic invariant by default |
| Local dimension | Implemented | Exact | From selected CAD cell closures |
| Algebraic singular locus | Implemented | Exact/partial | Boundary-polynomial/Jacobian semantics |
| Regular-open/closed operations | Implemented | Exact | Interior/closure identities |
| Exact measure/integration | Implemented | Exact/partial | Symbolic or exact unevaluated CAD integrals |
| Delineable curve/surface evaluation | Implemented | Numerical over exact branch descriptors | Root identity comes from CAD |
| CAD-cell simplicial meshing | Implemented | Numerical | Higher-dimensional Freudenthal subdivision |
| Adjacency-aware conforming mesh | Implemented | Numerical + structural conformity check | Shared sampled faces/vertices |
| Oriented cellular boundary maps | Unsupported | — | The cell complex stores unsigned incidence |
| Betti numbers/homology | Unsupported | — | Requires oriented chain-complex machinery |
| Exact roadmaps/path parameterizations | Unsupported | — | The connectivity graph is not a full algebraic roadmap |
| Certified isotopic triangulation | Unsupported | — | The mesh is not an ambient-isotopy certificate |
| General manifold/Whitney stratification | Unsupported | — | Singularity/local-dimension tools expose lower-level information |
