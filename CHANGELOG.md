# Changelog

This file contains information related to the development of the project, such as the features implemented in each version.

---

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




