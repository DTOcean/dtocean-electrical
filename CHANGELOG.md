# Change Log

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/)
and this project adheres to [Semantic Versioning](http://semver.org/).

## [Unreleased]

## [1.0.2] - 2017-03-28

### Added

- Add change log.

### Fixed

- Fixed issue where bathymetry points would be excluded if the local gradient
  exceeded 0. This now conforms to the user set limit.
- Fixed issue with the substation location algorithm, when devices where in
  symmetric layouts.

## [1.0.1] - 2017-03-09

### Fixed

- Replaced incorrect Shapely object when calculating design limits.

## [1.0.0] - 2017-01-05

### Added

- Initial import of dtocean-electrical from SETIS.
