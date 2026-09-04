# Package metadata

`semialg.__version__` exposes the installed package version. The authoritative
source-package value is the `project.version` field in `pyproject.toml`; package
builds and installed metadata should agree with that value.

Repository URLs, supported Python versions, dependencies, optional dependency
groups, and package classifiers are also defined in `pyproject.toml`.



## Primary API overview

This table is the substantive coverage target for the primary APIs assigned to this reference page. Each entry states the API's primary role; the family contract and detailed sections below explain shared algorithms, exactness guarantees, and limitations. It is maintained together with `docs/reference/primary_api_manifest.toml`, and documentation tests require every root-level primary API to map here rather than merely appearing in the generated public index.

| API | Kind | Role / return |
|---|---|---|
| `__version__` | str | Installed semialg package version string. |
