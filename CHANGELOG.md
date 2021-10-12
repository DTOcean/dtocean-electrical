# Change Log

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/)
and this project adheres to [Semantic Versioning](http://semver.org/).

## [2.1.0] - 2021-10-12

### Added

-   Add checks for paths between nodes which contain less than one grid point.

### Changed

-   Change behaviour of devices_per_string input so that the routing algorithm
    runs with only this value when set and not all values up to it.

### Fixed

-   Fix substation locating algorithm (for single substation).

## [2.0.0] - 2019-03-06

### Added

-   Re-implement constraints recording and plotting.
-   Add change log.
-   Add example which uses pickled inputs exported from dtocean-core.
-   Added edge_buffer option to ConfigurationOptions class to ensure substation
    is located a certain distance away from the lease area edge.
-   The selected export cable voltage is now returned with the outputs.

### Changed

-   Changed grid processing to improve methodology for joining lease area and
    cable corridor.
-   Excluded grid points are now removed from Grid.grid_pd and Grid.points as
    well as Grid.graph.
-   Emit a nx.NetworkXNoPath error in the optimiser if there are no point in the
    selected tool's graph.
-   Made calls to Dijkstra's shortest path algorithm more efficient for radial
    networks.
-   Accelerated grid distance and gradient processing by converting parts of the
    grid DataFrame to numpy arrays and reorganising other DataFrame
    manipulation.
-   Changed networkx single_source_dijkstra function call to match 2.0 API.
-   Optimise by lowest cost per power transmitted, rather than just cost.
-   Allow array rated powers above the maximum scope to be calculated, but still
    raise a warning.
-   Array and export cable databases are now provided as separate inputs.

### Removed

-   Removed incorrect catch for KeyErrors when searching for cable routes.

### Fixed

-   Replaced incorrect Shapely object when calculating design limits.
-   Fixed issue where bathymetry points would be excluded if the local gradient
    exceeded 0. This now conforms to the user set limit.
-   Fixed issue with the substation location algorithm, when devices where in
    symmetric layouts.
-   Fixed bug with checking soil types.
-   Use power histogram bin centres for power calculations rather than edges.
-   Add missing soil types for setting burial protection index.
-   Ensure all Grid attributes are referenced by the point "id".
-   Fixed Network.make_cable_routes so that cable points are returned in the
    correct order.
-   Improved algorithm for determining the grid spacing.
-   Ensured that umbilical cable paths are correctly joined to the inter-arrays
    cables.
-   Fixed approximation of umbilical cable termination point when path only has
    one point.
-   Fixed umbilical cable length calculation by changing device rotation frame
    from global to local coordinates.
-   Improved speed of umbilical calculation.
-   Ensure that the burial depth can be calculated if cable paths intersect.
-   Fixed issue where umbilical cables would incorrectly reduce the static cable
    lengths.
-   Improved memory efficiency for grid processing.

## [1.0.0] - 2017-01-05

### Added

-   Initial import of dtocean-electrical from SETIS.
