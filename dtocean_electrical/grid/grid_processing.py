# -*- coding: utf-8 -*-
"""
Functions to process the bathymetry data into a usual graph for cable routing
algorithms.

author:: Adam Collin <a.collin@ed.ac.uk>

"""
# Start logging
from datetime import datetime
from shapely.geometry import Polygon
from grid import Grid, GridPoint
import networkx as nx
import numpy as np
from scipy.spatial import distance
import logging
module_logger = logging.getLogger(__name__)


def grid_processing(site_data, export_data, options, n_neighbours=8):

    '''Function to control bathymetry processing stages. This first creates
    a merged bathymetry and then creates a Grid object.

    Args:
        site_data (object) [-]: Instance of the ElectricalSiteData class.
        export_data (object) [-]: Instance of the ElectricalExportData class.
        options (object) [-]: Instance of the ConfigurationOptions class.
        n_neighbours (int) [-]: Number of neighbours to consider when building
            network graph.

    Attributes:
        exclusion_zones (list) [-]: List of Shapely Polygon objects.
        merged_grid (pd.DataFrame) [-]: Site and export data merged into one
            area.
        grid (object) [-]: Instance of the Grid class.

    Returns:
        grid

    Note:
        n_neighbours is set to eight by default. Allowed values are 4 and 8.

    '''

    exclusion_zones = get_exclusions(site_data, export_data)

    merged_grid = merge_bathymetry(site_data, export_data)

    grid = \
        make_graph_object(merged_grid, n_neighbours, exclusion_zones, options)

    grid.lease_boundary = lease_area_boundary(site_data)

    module_logger.info("Grid prepared...")

    return grid


def lease_area_boundary(site_data):

    '''Get the lease area boundary.

    Notes:
        This is to be replaced with site.lease_boundary so update interface.

    '''

    # make shapely point list of lease area
    lease_area_points_x = site_data.bathymetry.x.tolist()
    lease_area_points_y = site_data.bathymetry.y.tolist()

    # This assumes a four sided square/rectangle shape

    lease_min_x = min(lease_area_points_x)
    lease_max_x = max(lease_area_points_x)
    lease_min_y = min(lease_area_points_y)
    lease_max_y = max(lease_area_points_y)

    polygon_edge = [[lease_min_x, lease_min_y],
                    [lease_min_x, lease_max_y],
                    [lease_max_x, lease_max_y],
                    [lease_max_x, lease_min_y]]

    lease_area_polygon = Polygon(polygon_edge)

    return lease_area_polygon


def merge_bathymetry(site_data, export_data):

    '''Merge two bathymetry data sets.

    Args:
        site_data (object) [-]: Instance of the ElectricalSiteData class.
        export_data (object) [-]: Instance of the ElectricalExportData class.

    Attributes:
        export_area (set) [m]: Set of x, y coordinate tuples of export area.
        lease_area (set) [m]: Set of x, y coordinate tuples of lease area.
        shared_area (list) [m]: List of x, y coordinate tuples present in lease
            area and export area.
        shared_x (tuple) [m]: Tuple of x coordinates of points.
        shared_y (tuple) [m]: Tuple of y coordinates of points.
        site_index (pd.DataFrame) [-]: i and j indices of site from shared
            area.
        export_index (pd.DataFrame) [-]: i and j indices of export from shared
            area.
        new_site_data (pd.DataFrame) [-]: Updated data for site area.
        merged_grid (pd.DataFrame) [-]: Merge of export and site areas.

    Return:
        merged_grid

    '''

    module_logger.info("Merging bathymetry...")

    export_area = set(zip(export_data.bathymetry.x, export_data.bathymetry.y))

    lease_area = set(zip(site_data.bathymetry.x, site_data.bathymetry.y))

    shared_area = \
        merge_bathymetry_xy_intersection_test(export_area, lease_area)

    shared_x, shared_y = zip(*[(point[0], point[1]) for point in shared_area])

    merge_bathymetry_z_intersection_test(
        shared_x, shared_y, export_data.bathymetry, site_data.bathymetry)

    site_index = \
        get_index_values_of_selection(site_data.bathymetry, shared_x, shared_y)

    export_index = \
        get_index_values_of_selection(
            export_data.bathymetry, shared_x, shared_y)

    new_site_data = \
        offset_index_values(site_data.bathymetry, site_index, export_index)

    # remove common points from lease area
    new_site_data.drop(
        new_site_data.index[site_index.index], inplace=True)

    merged_grid = export_data.bathymetry.append(site_data.bathymetry)

    merged_grid['id'] = list(range(len(merged_grid)))
    merged_grid = merged_grid.reset_index(drop=True)

    return merged_grid


def make_graph_object(area_pd, n_neighbours, exclusions, options):

    '''Creates instance of the Grid class and also generates the NetworkX graph
    object.

    Args:
        area_pd (pd.DataFrame) [-]: Area to be converted.
        n_neighbours (int) [-]: Number of neighbour directions to consider.
        exclusions (list) [-]:  List of Shapely Polygon objects.
        options (object) [-]: Instance of ConfigurationOptions class.

    Attributes:
        grid (object) [-]: Instance of Grid object.
        grid_dict (dict) [-]: Grid represented as dictionary,
            Key: i, j indices; Values: index, x coord,  y coord, z coord.
        graph_ready (dict) [m]: Structured dictionary for NetworkX of edge
            information for all points.
        grads (dict) [degrees] Structured dictionary of gradients.

    returns:
        grid

    '''

    grid = Grid(area_pd)

    module_logger.info("Creating bathymetry grid points...")

    grid.add_points_to_grid()

    grid_dict = grid.grid_as_dict()

    graph_ready = get_neighbours_distance(grid_dict)

    grads = get_neighbours_gradient(grid_dict)

    module_logger.info("Building network graph...")

    grid.graph = nx.Graph(graph_ready)

    remove_exclusion_zones(exclusions, grid)

    apply_gradient_constraint(options.equipment_gradient_constraint,
                              grid,
                              grads)

    apply_equipment_constraints(options, grid,
                                preselected=options.installation_tool)

    module_logger.info("Constraints checked...")

    return grid


def merge_bathymetry_xy_intersection_test(area_one, area_two):

    '''Test for intersection in the x - y plane.

    Args:
        area_one (set) [m]: Set of x, y coordinate tuples of an area.
        area_two (set) [m]: Set of x, y coordinate tuples of an area.

    Attributes:
        common (list) [m]: List of x, y coordinate tuples present in area_one
            and area_two.

    Returns:
        common

    Notes:
        Tests have been separated to provide better error feedback to the user.

    '''

    common = list(area_one.intersection(area_two))

    if not common:

        errStr = ("No intersection was found between the lease area "
                  "and cable corridor")

        raise RuntimeError(errStr)

    return common


def merge_bathymetry_z_intersection_test(x_list, y_list, area_one,
                                         area_two):

    '''Test for intersection in the z direction.

    Args:
        x_list (tuple) [m]: Tuple of x coordinates of points.
        y_list (tuple) [m]: Tuple of y coordinates of points.
        area_one (pd.DataFrame) [-]: First area under consideration.
        area_two (pd.DataFrame) [-]: Second area under consideration.

    Attributes:
        area_one_z (pd.Series) [-]: First area z coordinates.
        area_two_z (pd.Series) [-]: Second area z coordinates.
        err_z (list) [m]: Difference between area_one_z and area_two_z.

    Notes:
        Tests have been separated to provide better error feedback to the user.

    '''

    area_one_z = get_z_dimension(area_one, x_list, y_list)

    area_two_z = get_z_dimension(area_two, x_list, y_list)

    err_z = [z_one - z_two for z_one, z_two in zip(area_one_z, area_two_z)]

    if sum(err_z) != 0:

        errStr = ("Error in depth at intersection between the lease "
                  "area and the cable corridor.")

        raise RuntimeError(errStr)

    return


def get_z_dimension(area, x_list, y_list):

    '''Get z coordinate of area at points defined by x_list and y_list. The z
    coordinate is defined in column 'layer 1 start'.

    Args:
        area (pd.DataFrame) [-]: Area under consideration.
        x_list (tuple) [m]: Tuple of x coordinates of points.
        y_list (tuple) [m]: Tuple of y coordinates of points.

    Attributes:
        z (pd.Series) [m]: z coordinates.

    Returns:
        z

    '''

    z = (area[(area['x'].isin(x_list)) &
              (area['y'].isin(y_list))]['layer 1 start'])

    return z


def get_index_values_of_selection(area, x_list, y_list):

    '''Get index values of bathymetry pd.DataFrame.

    Args:
        area (pd.DataFrame) [-]: Area representation which must include column
            headers x, y, i and j.
        x_list (tuple) [m]: Tuple of x coordinates of points.
        y_list (tuple) [m]: Tuple of y coordinates of points.

    Attributes:
        index_values (pd.DataFrame) [-]: Returned i and j indices in selection
            area.

    Returns:
        index_values

    '''

    index_values = \
        area[(area['x'].isin(x_list)) & (area['y'].isin(y_list))][['i', 'j']]

    return index_values


def offset_index_values(area, area_indices, other_indices):

    '''Offset the i and j index values in an area with respect to index values
    in a neighbour area.

    Args:
        area (pd.DataFrame) [-]: Area under consideration.
        area_indices (pd.DataFrame) [-]: i and j indices from shared area in
            area.
        other_indices (pd.DataFrame) [-]: i and j indices from shared area in
            other area.

    Attributes:
        area_i (list) [-]: List of i indices from shared area in area.
        area_j (list) [-]: List of i indices from shared area in area.
        other_i (list) [-]: List of i indices from shared area in neighbour
            area.
        other_j (list) [-]: List of i indices from shared area in neighbour
            area.
        diff_i (int) [-]: Difference in i indices; acts as offset.
        diff_j (int) [-]: Difference in j indices; acts as offset.
        new_i (list) [-]: Updated list of i indices from shared area in area.
        new_j (list) [-]: Updated list of i indices from shared area in area.

    Returns:

    '''

    area_i = area_indices.i.tolist()
    area_j = area_indices.j.tolist()
    other_i = other_indices.i.tolist()
    other_j = other_indices.j.tolist()

    # find difference in i and j
    diff_i = [(export - lease) for (export, lease) in zip(other_i, area_i)][0]
    diff_j = [(export - lease) for (export, lease) in zip(other_j, area_j)][0]

    # add difference to area
    new_i = [point + diff_i for point in area.i]
    new_j = [point + diff_j for point in area.j]

    # update area pd.DataFrame
    area.drop(['i', 'j'], 1, inplace=True)
    area['i'] = new_i
    area['j'] = new_j

    return area


def get_neighbours_distance(grid_dict):

    '''Get distance between neighbouring points for all points defined in
    grid_dict.

    Args:
        grid_dict(dict) [-]: Grid represented as dictionary,
            Key: i, j indices; Values: index, x coord,  y coord, z coord.

    Attributes:
        edge_weights (dict) [-]: Structured dictionary for NetworkX.

    Return:
        edge_weights.

    '''

    module_logger.info("Checking distance to neighbours...")

    edge_weights = {point[0]: get_weights(key, grid_dict)
                    for key, point in grid_dict.iteritems()}

    return edge_weights


def get_neighbours_gradient(grid_dict):

    '''Get gradient between neighbouring points for all points defined in
    grid_dict.

     Args:
        grid_dict(dict) [-]: Grid represented as dictionary,
            Key: i, j indices; Values: index, x coord,  y coord, z coord.

    Attributes:
        gradients (dict) [-]: Structured dictionary of gradients.

    Return:
        gradients.

    '''

    module_logger.info("Checking gradient to neighbours...")

    gradients = {point[0]: get_gradients(key, grid_dict)
                 for key, point in grid_dict.iteritems()}

    return gradients


def get_exclusions(site_data, export_data):

    '''Collect exclusion zones from site and export area together.

    Args:
        site_data (object) [-]: Instance of ElectricalSiteData object.
        export_data (object) [-]: Instance of ElectricalExportData object.

    Attributes:
        all_exclusions (list) [-]:  List of Shapely Polygon objects.

    Returns:
        all_exclusions

    Notes:
        Site and export exclusion_zones attribute are lists of Polygon objects.

    '''

    all_exclusions = []

    if site_data.exclusion_zones is not None:

        all_exclusions.extend(site_data.exclusion_zones)

    if export_data.exclusion_zones is not None:

        all_exclusions.extend(export_data.exclusion_zones)

    return all_exclusions


def remove_exclusion_zones(all_exclusions, grid):

    '''Remove exclusion zones from the area. Call the remove_exclusion_zones
    method from the Grid class.

    Args:
        all_exclusions (list) [-]:  List of Shapely Polygon objects.
        grid (object) [-]: Instance of Grid object.

    Attributes:
        excluded_points (list) [-]: Index of points to be removed.

    '''

    module_logger.info("Checking for exclusion zones...")

    if all_exclusions:

        excluded_points = grid.get_exclusion_zone_points(all_exclusions)

        grid.graph = grid.remove_points_from_graph(excluded_points, grid.graph)

        msg = ("Number of points removed in exclusion zones: {}".format(
               len(excluded_points)))

        module_logger.info(msg)

    return


def apply_gradient_constraint(maximum_gradient, grid, grads):

    '''Remove points which breach gradient constraint. Call gradient_constraint
    method from Grid class.

    Args:
        maximum_gradient (float) [deg]: Gradient constraint to be applied.
        grid (object) [-]: Instance of Grid object.
        grads (dict) [degrees] Structured dictionary of gradients.

    '''

    module_logger.info("Checking gradient constraints...")

    grid.gradient_constraint(maximum_gradient, grads)

    return


def apply_equipment_constraints(options,
                                grid,
                                constraint_type=2,
                                preselected=None):

    '''Find areas compatible with installation equipment. Considers two
    approaches to defining area. Type # 1 finds installers which are valid at
    all points. Type # 2 creates graphs for each installer individually.

    Args:
        options (object) [-]: Instance of ConfigurationOptions object.
        grid (object) [-]: Instance of Grid object.
        constraint_type (int) [-]: Defines the type of seabed filtering to be
            applied.

    Attributes:
        compatability_matrix (pd.DataFrame) [-]: Equipment soil compatability
            matrix.
        compatibility_matrix_dict (dict) [-]: Equipment soil compatability
            matrix in dictionary structure.
        valid_installers (list) [-]: List of valid installers for the grid.

    Returns:
        valid_installers
        
    Notes:
        Allow the user to enter only one tool for now.

    '''

    module_logger.info("Checking equipment soil compatibility...")
    
    compatability_matrix = options.equipment_soil_compatibility
    
    compatibility_matrix_dict = options.equipment_soil_compatibility_dict


    if preselected:

        valid_installers = preselected

        soils = compatibility_matrix_dict[preselected]

        grid.check_equipment_soil_compatibility(preselected, soils)

    elif constraint_type == 1:

        # compatability_matrix = options.equipment_soil_compatibility

        valid_installers = (
            grid.check_equipment_soil_compatibility_site(compatability_matrix))

    elif constraint_type == 2:

        valid_installers = []

        # compatibility_matrix_dict = options.equipment_soil_compatibility_dict

        for technique, soils in compatibility_matrix_dict.iteritems():

            grid.check_equipment_soil_compatibility(technique, soils)

            valid_installers.append(technique)

    return valid_installers


def get_gradients(indices, points):

    '''Get gradients between point defined by indices and all neighbours.

    Args:
        indices (tuple) [-]: i and j indices of point under consideration.
        points (dict) [-]: Grid represented as dictionary,
            Key: i, j indices; Values: index, x coord,  y coord, z coord.

    Attributes:
        p1 (tuple) [-]: Coordinates of point under consideration, x, y, z.
        local_list (list) [-]: List of neighbour points x, y, z coords.
        markers (list) [-]: List of neighbour points index values.
        gradients (list) [m]: Edge gradients.
        gradients_dict (dict) [-]: Structured dictionary of gradients.

    Return:
        gradients_dict

    '''

    p1 = (points[indices][1], points[indices][2], points[indices][3])

    local_list, markers = get_adjacent_points(indices, points)

    gradients = [gradient(p1, p2) for p2 in local_list]

    gradients_dict = make_weight_dict(gradients, markers)

    return gradients_dict


def get_weights(indices, points):

    '''Get distance between point defined by indices and all neighbours.

    Args:
        indices (tuple) [-]: i and j indices of point under consideration.
        points (dict) [-]: Grid represented as dictionary,
            Key: i, j indices; Values: index, x coord,  y coord, z coord.

    Attributes:
        p1 (tuple) [-]: Coordinates of point under consideration, x, y, z.
        local_list (list) [-]: List of neighbour points x, y, z coords.
        markers (list) [-]: List of neighbour points index values.
        edges (list) [m]: Edge lengths.
        edges_dict (dict) [-]: Structured dictionary for NetworkX.

    Return:
        edges_dict

    '''

    p1 = (points[indices][1], points[indices][2], points[indices][3])

    local_list, markers = get_adjacent_points(indices, points)

    edges = [edge_length(p1, p2) for p2 in local_list]

    edges_dict = make_weight_dict(edges, markers)

    return edges_dict


def get_adjacent_points(indices, points):

    '''Find neighbours in group points of point defined by indices.

    Args:
        indices (tuple) [-]: i and j indices of point under consideration.
        points (dict) [-]: Grid represented as dictionary,
            Key: i, j indices; Values: index, x coord,  y coord, z coord.

    Attributes:
        i_index (int) [-]: i index of point under consideration.
        j_index (int) [-]: j index of point under consideration.
        local_list (list) [-]: List of neighbour points x, y, z coords.
        markers (list) [-]: List of neighbour points index values.

    Return:
        local_list
        markers

    Note:
        based on http://stackoverflow.com/questions/2373306/

    '''

    i_index = indices[0]
    j_index = indices[1]

    local_list = []
    append = local_list.append
    markers = []
    append_markers = markers.append

    for x, y in [(i_index+i, j_index+j)
                 for i in (-1, 0, 1) for j in (-1, 0, 1) if i != 0 or j != 0]:

        if (x, y) in points:

            append((points[(x, y)][1:]))
            append_markers(points[(x, y)][0])

    return local_list, markers


def make_weight_dict(edges, markers):

    '''Convert edge data into format required for networkx graph object.

    Args:
        edges (list) [m]: Edge lengths.
        markers (list) [-]: Index of points upon which the edges terminate.

    Attributes:
        weight_dict (dict) [-]: Structured dictionary.

    Return:
        weight_dict.

    '''

    weight_dict = \
        {marker: {'weight': edge} for marker, edge in zip(markers, edges)}

    return weight_dict


def edge_length(p1, p2):

    '''Find the distance between two points.

    Args:
        p1 (tuple) [m]: Point one x, y, z coordinates.
        p1 (tuple) [m]: Point two x, y, z coordinates.

    Attributes:
        edge_length (float) [m]: Distance.

    Return:
        edge_length.

    '''

    edge_length = distance.euclidean((p1[0], p1[1], p1[2]),
                                     (p2[0], p2[1], p2[2]))

    return edge_length


def gradient(p1, p2):

    '''Find the gradient between two points.

    Args:
        p1 (tuple) [m]: Point one x, y, z coordinates.
        p1 (tuple) [m]: Point two x, y, z coordinates.

    Attributes:
        delta_x (float) [m]: Distance in x direction.
        delta_y (float) [m]: Distance in y direction.
        delta_z (float) [m]: Distance in z direction.
        horizontal_distance (float) [m]: Distance in x-y plane.
        angle (float) [rad]: Gradient expressed in rads.

    Return:
        angle (float) [degrees]: Gradient expressed in degrees.

    '''

    delta_x = p2[0] - p1[0]
    delta_y = p2[1] - p1[1]
    delta_z = p2[2] - p1[2]

    horizontal_distance = np.sqrt(np.square(delta_x)+np.square(delta_y))
    angle = np.arctan(delta_z/horizontal_distance)

    return abs(np.degrees(angle))
