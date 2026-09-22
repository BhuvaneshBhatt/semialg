# Moments, centroids, and covariance

The moment APIs are wrappers around `integrate_over_region`, so their supported cases follow the integration engine.

## Public API

- `region_moment`
- `centroid`
- `covariance_matrix`

## Examples

```python
import sympy as sp
from semialg import region_moment, centroid, covariance_matrix

x, y = sp.symbols("x y", real=True)

region_moment(x**2 <= 1, [x], powers=[2])
# 2/3

centroid(sp.And(x >= 0, y >= 0, x + y <= 1), [x, y])
# {x: 1/3, y: 1/3}

covariance_matrix(sp.And(x >= 0, x <= 1, y >= 0, y <= 1), [x, y])
# Matrix([[1/12, 0], [0, 1/12]])
```

## Notes

Centroid and covariance require finite positive measure in the selected measure dimension. Unsupported or infinite-measure cases should fail conservatively.
