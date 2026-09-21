# Changelog

This file contains information related to the development of the project, such as the features implemented in each version.

---

## [Unreleased]

## [0.3.0a1] - 2026-09-20

### Changed
- Promoted the stabilized 0.3 development line to its first alpha release.
- Established the current package-root API, solver lifecycle, result semantics, topology validation, reporting behavior, executable examples, and optional mesh subsystem as the alpha baseline.

### Notes
- This is a pre-release intended for broader real-world testing before 0.3.0.
- Backward-incompatible API changes may still occur before the final 0.3.0 release.

## [0.3.0.dev1] - 2026-09-20

### Changed
- Enforced model membership invariants when adding elements and unified `simple_report()` across all five public model types.
- Added common pre-assembly topology validation, including explicit errors for empty models and orphan nodes.
- Modernized executable FEM examples around explicit public imports and test-friendly model builders.
- Hardened the optional mesh subsystem around explicit Gmsh execution, temporary-file cleanup, meshio-based triangle loading, and clearer errors.
- Corrected full-circle geometry generation to use four valid Gmsh quarter-circle arcs.
- Renamed the public surface-hole helper to `subtract_surfaces()` and removed the invalid unused arc-surface helper.

### Removed
- Removed unsupported legacy modules for ad-hoc I/O, JSON model loading, frozen plotting helpers, incomplete 3D experiments, obsolete report templates, and unintegrated material/section helpers.
- Removed the legacy `.nusa` sample data and materials example that depended on those modules.
- Kept `nusa.mesh` as the supported optional peripheral subsystem.

### Tests
- Added API-integrity, topology, example-smoke, and mesh regression coverage.
- Expanded the suite to 139 passing tests, with GitHub Actions also validating package builds and Sphinx documentation.

## [0.3.0.dev0] - 2026-09-19

### Changed
- Recovered active development and stabilized the finite-element solver architecture.
- Replaced legacy nested solver state with vector-based global force and displacement state.
- Introduced explicit `assemble()` / `stiffness_matrix` semantics and separated assembly invalidation from solution invalidation.
- Decoupled internal solver indexing from public node labels.
- Normalized the visualization API and clarified applied loads, generalized nodal forces, and reactions.
- Defined an explicit top-level public API instead of wildcard package exports.
- Modernized packaging around `pyproject.toml`, optional mesh dependencies, and automated CI.

### Fixed
- Added explicit singular-system errors instead of least-squares fallback.
- Corrected nonzero prescribed-displacement handling in the reduced system.
- Corrected CST orientation handling using signed area.
- Improved topology-state invalidation and persistence of loads and boundary conditions.

### Tests
- Expanded numerical and API regression coverage across Spring, Bar, Truss, Beam, and LinearTriangle models.

## [0.3.dev0] - 2020-09-02

### Added
- Created `model_reader` function to generate models from text files with JSON structure. Currently implemented for Spring and Truss element models.

## [0.3.dev0] - 2020-09-01

### Added
- Added `simple_report` for Spring-type element.
- Created and added Nusa logo.

## [0.3.dev0] - 2020-08-15

### Added
- Added `version.py`.
- Adjusted `setup.py`.

---

## [0.2.0] - 2018-11-16

### Fixed
- Minor bugs in beam examples, which were outdated with respect to the current version.

---

## [0.1.0] - 2016-05-19

### Changed
- Displacements and forces are now stored in a dictionary for better control over components and improved readability.  
  Currently implemented only in the `BeamModel`.

## [0.1.0] - 2016-01-10

### Changed
- Method names changed from "mixedUp" style to `lower_case_with_underscores` as recommended by PEP8.

## [0.1.0] - 2016-06-01

### Changed
- Updated containers from list to dictionary for Spring-type models.
- Tested with existing examples. Force implementation for elements still pending.
