# Change Log

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/)
and this project adheres to [Semantic Versioning](http://semver.org/).

## [Unreleased]

### Added

- Reimplement constraints recording and plotting.
- Add change log.

### Changed

- Changed grid processing to improve methodology for joining lease area and
  cable corridor.

### Fixed

- Replaced incorrect Shapely object when calculating design limits.
- Fixed issue where bathymetry points would be excluded if the local gradient
  exceeded 0. This now conforms to the user set limit.
- Fixed issue with the substation location algorithm, when devices where in
  symmetric layouts.
- Fixed bug with checking soil types.

### Removed

- Removed incorrect catch for KeyErrors when searching for cable routes.

## [1.0.0] - 2017-01-05

### Added

- Initial import of dtocean-electrical from SETIS.
