# Known issues

This page tracks concrete implementation limitations that are narrower than the mathematical scope described in [Limitations](limitations.md).

## CAD resource controls

Some CAD interfaces validate resource-control arguments such as cell or time limits more broadly than the underlying algorithms can currently enforce as hard interruption budgets. Do not rely on these options as a process-level timeout mechanism.

## Numerical meshing is not a topology certificate

Numerical triangulation and conforming-mesh helpers are useful for visualization and downstream numerical work, but their output is not an ambient-isotopy certificate. Use the exact topology APIs when the topological claim itself must be certified.
