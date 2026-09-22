# API naming

The package root is reserved for common mathematical questions and broadly useful region and geometry operations. Specialized result types, certificate machinery, lower-level reasoning primitives, and symbolic relation constructors live in expert namespaces.

## Region relations and operations

Boolean questions use `is_subset`, `is_equal`, `is_disjoint`, `is_interior_disjoint`, and `intersects`. Lower-level decision implementations such as `region_subset`, `region_equal`, and `region_disjoint` are expert APIs. Symbolic relation constructors use `RegionSubset`, `RegionEqual`, and `RegionDisjoint` for unevaluated assertions.

Set constructors use the `region_*` prefix: `region_union`, `region_intersection`, `region_difference`, `region_symmetric_difference`, `region_complement`, and `region_product`. Topological constructors follow the same convention with `region_closure`, `region_interior`, and `region_boundary`. The composition names `closure_of_interior` and `interior_of_closure` state their operation order directly.

## Components, moments, and measure

`connected_components` is the primary connected-component operation and uses certified connectivity. `semialg.regions.explicit_region_components` is a specialist helper for syntactically explicit pieces when full connectivity certification is not required.

`centroid` and `covariance_matrix` are the normalized moment APIs. `region_measure` is abstraction-aware: canonical geometry defaults to intrinsic measure and formula regions to ambient measure. `semialgebraic_measure` is formula-oriented and also supports explicit bounds and parameter stratification.

## Polynomial nonnegativity

`polynomial_nonnegative` is the primary entry point. It returns a Boolean by default and accepts `return_result=True` for the certified strategy trace. SOS, Zeng, ARS, and CAD selection and witness-search tuning are options on the same operation.

## Images and preimages

`affine_image` and `affine_preimage` dispatch by representation and preserve canonical structure when possible. `linear_image` is the zero-offset specialization. `region_image` and `region_preimage` handle general symbolic maps. Formula-level `semialgebraic_image` and `semialgebraic_preimage` remain in `semialg.geometry_queries` for callers that need direct formula behavior.

## Boolean structured regions

`BooleanRegion` represents Boolean combinations of `StandardRegion` objects. Its structured constructors are `BooleanRegion.union`, `.intersection`, `.difference`, and `.symmetric_difference`; the lowercase `region_*` operations remain representation-independent.

## Small specialized entry points

Some thin functions encode standard mathematical vocabulary rather than duplicate semantics. `matrix_pd_on` and `matrix_psd_on` specialize `matrix_definiteness`; `linear_image` specializes `affine_image`; and `is_connected` shares the semialgebraic connectedness contract with `is_path_connected`.
