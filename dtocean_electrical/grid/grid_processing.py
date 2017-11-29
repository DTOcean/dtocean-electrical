# -*- coding: utf-8 -*-
"""
Functions to process the bathymetry data into a usual graph for cable routing
algorithms.

author:: Adam Collin <a.collin@ed.ac.uk>,
         Mathew Topper <dataonlygreater@gmail.com>

"""
import pickle
import logging

import numpy as np
import pandas as pd
import networkx as nx
from scipy.spatial import distance
from shapely.geometry import MultiPoint, Point

from .grid import Grid

# Start logging
module_logger = logging.getLogger(__name__)


def grid_processing(site_data, export_data, options):

    '''Function to control bathymetry processing stages. This first creates
    a merged bathymetry and then creates a Grid object.

    Args:
        site_data (object) [-]: Instance of the ElectricalSiteData class.
        export_data (object) [-]: Instance of the ElectricalExportData class.
        options (object) [-]: Instance of the ConfigurationOptions class.

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

    clipped_export, lease_polygon = clip_grid(site_data.bathymetry,
                                              export_data.bathymetry)
    
    network_graph = make_graph(site_data.bathymetry, clipped_export)
    gradient_dict = make_gradient_dict(site_data.bathymetry, clipped_export)
    
    #    # Dump these for hot start
    #    nx.write_gpickle(network_graph, "network_graph.pkl")
    #    
    #    with open("gradient_dict.pkl", "wb") as fstream:
    #        pickle.dump(gradient_dict, fstream)
    #        
    #    # Load these for hot start
    #    network_graph = nx.read_gpickle("network_graph.pkl")
    #    
    #    with open("gradient_dict.pkl", "rb") as fstream:
    #        gradient_dict = pickle.load(fstream)

    grid = make_grid(site_data.bathymetry,
                     clipped_export,
                     network_graph,
                     lease_polygon)
    
    grid.remove_exclusion_zones(exclusion_zones)
    constrained_lines = grid.gradient_constraint(
                                        options.equipment_gradient_constraint,
                                        gradient_dict)
    
    # Free memory
    del(gradient_dict)

    apply_equipment_constraints(grid,
                                options, 
                                preselected=options.installation_tool)

    return grid, exclusion_zones, constrained_lines


def clip_grid(grid_df_static, grid_df_to_clip): 

    """Remove grid_df_static from grid_df_to_clip, leaving some overlapping
    points. Return a copy of grid_df_to_clip with points removed."""
    
    # Get a convex hull grid_df_static
    static_coords = [(x, y) for x, y in zip(grid_df_static.x,
                                            grid_df_static.y)]
    static_multipoint = MultiPoint(static_coords)
    static_poly = static_multipoint.convex_hull
    
    while True:
        
        # Remove all points *inside* grid_df_static from grid_df_to_clip
        inside = map(lambda x, y: static_poly.contains(Point(x, y)),
                     grid_df_to_clip.x,
                     grid_df_to_clip.y)
        
        grid_df_clipped = grid_df_to_clip[~np.array(inside)]
        clipped_coords = [(x, y) for x, y in zip(grid_df_clipped.x,
                                                 grid_df_clipped.y)]
        
        shared_coords = list(set(clipped_coords).intersection(
                                                        set(static_coords)))
        
        # If we have some values then leave the loop
        if shared_coords: break
        
        if static_poly.area <= 0.:
            
            errStr = ("No overlapping nodes were found between lease area "
                      "and cable corridor.")
            raise RuntimeError(errStr)
        
        # Get maximum grid_df_to_clip spacing
        x_unique = sorted(set(grid_df_to_clip.x))
        y_unique = sorted(set(grid_df_to_clip.y))
        
        dx = x_unique[1] - x_unique[0]
        dy = y_unique[1] - y_unique[0]
        
        maxd = max(dx, dy)
        
        # Reduce the static_grid polygon
        static_poly = static_poly.buffer(-maxd / 2.)
        
    # Offset the ids of the clipped_grid to avoid clash with the grid_df_static
    offset_index = max(grid_df_static.id) + 1
    grid_df_clipped.loc[:, 'id'] = grid_df_clipped['id'].apply(
                                                    lambda x: x + offset_index)
        
    return grid_df_clipped, static_poly


def make_graph(grid_df_static, grid_df_clipped):

    '''Creates instance NetworkX graph object connecting the two grids.
    '''

    module_logger.info("Calculating lease area distances...")

    id_indexed_grid_df = grid_df_static.set_index("id")
    id_array, x_array, y_array, z_array = make_grid_arrays(grid_df_static)
    graph_dict_static = get_neighbours_distance(grid_df_static,
                                                id_indexed_grid_df,
                                                id_array,
                                                x_array,
                                                y_array,
                                                z_array)
    
    module_logger.info("Calculating cable corridor distances...")

    id_indexed_grid_df = grid_df_clipped.set_index("id")
    id_array, x_array, y_array, z_array = make_grid_arrays(grid_df_clipped)
    graph_dict_clipped = get_neighbours_distance(grid_df_clipped,
                                                 id_indexed_grid_df,
                                                 id_array,
                                                 x_array,
                                                 y_array,
                                                 z_array)
    
    assert (set(graph_dict_static.keys()) &
            set(graph_dict_clipped.keys())) == set()
    
    # Get dataframe of overlapping points
    overlap_df = pd.merge(grid_df_static,
                          grid_df_clipped,
                          on=["x", "y"])
    
    # Create zero length edges in the lease area dict
    add_overlap_edges(graph_dict_static, overlap_df)
    
    assert (set(graph_dict_static.keys()) &
            set(graph_dict_clipped.keys())) == set()
    
    # Merge the two dicts
    graph_dict_final = dict(graph_dict_static, **graph_dict_clipped)
    
    module_logger.info("Building network graph...")
    graph = nx.Graph(graph_dict_final)
    
    return graph
    
    
def add_overlap_edges(graph_dict_static, overlap_df):
    
    for _, point in overlap_df.iterrows():
        graph_dict_static[point.id_x][point.id_y] = {'weight': 0.}
        
    return


def make_gradient_dict(grid_df_static, grid_df_clipped):

    '''Creates NetworkX graph input dictionary for gradients over the two
    grids.
    '''
    
    id_indexed_grid_df = grid_df_static.set_index("id")
    id_array, x_array, y_array, z_array = make_grid_arrays(grid_df_static)

    module_logger.info("Calculating lease area gradients...")
    graph_dict_static = get_neighbours_gradient(grid_df_static,
                                                id_indexed_grid_df,
                                                id_array,
                                                x_array,
                                                y_array,
                                                z_array)
    
    id_indexed_grid_df = grid_df_clipped.set_index("id")
    id_array, x_array, y_array, z_array = make_grid_arrays(grid_df_clipped)

    module_logger.info("Calculating cable corridor gradients...")
    graph_dict_clipped = get_neighbours_gradient(grid_df_clipped,
                                                 id_indexed_grid_df,
                                                 id_array,
                                                 x_array,
                                                 y_array,
                                                 z_array)
    
    assert (set(graph_dict_static.keys()) &
            set(graph_dict_clipped.keys())) == set()
    
    # Merge the two dicts
    graph_dict_final = dict(graph_dict_static, **graph_dict_clipped)
    
    return graph_dict_final


def make_grid(grid_df_static,
              grid_df_clipped,
              network_graph,
              lease_polygon):

    '''Creates instance of the Grid class.
    '''
    
    module_logger.info("Preparing grid...")
    
    # Drop the i,j columns and join the two grids
    grid_df_static = grid_df_static.drop(["i", "j"], axis=1)
    grid_df_clipped = grid_df_clipped.drop(["i", "j"], axis=1)
    grid_df_joined = grid_df_static.append(grid_df_clipped, ignore_index=True)
    
    # Initalise the Grid object
    grid = Grid(grid_df_joined,
                network_graph,
                lease_polygon)
    
    return grid


def get_neighbours_distance(grid_df,
                            id_indexed_grid_df,
                            id_array,
                            x_array,
                            y_array,
                            z_array):

    '''Get distance between neighbouring points for all points defined in
    grid_df.

    Return:
        edge_weights (dict) [-]: Structured dictionary for NetworkX.

    '''

    edge_weights = {point.id: get_metric_weights(point.id,
                                                 id_indexed_grid_df,
                                                 id_array,
                                                 x_array,
                                                 y_array,
                                                 z_array,
                                                 edge_length)
                                            for _, point in grid_df.iterrows()}

    return edge_weights


def get_neighbours_gradient(grid_df,
                            id_indexed_grid_df,
                            id_array,
                            x_array,
                            y_array,
                            z_array):

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

    edge_gradients = {point.id: get_metric_weights(point.id,
                                                   id_indexed_grid_df,
                                                   id_array,
                                                   x_array,
                                                   y_array,
                                                   z_array,
                                                   gradient)
                                            for _, point in grid_df.iterrows()}

    return edge_gradients


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


def apply_equipment_constraints(grid,
                                options,
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


def get_metric_weights(index,
                       id_indexed_grid_df,
                       id_array,
                       x_array,
                       y_array,
                       z_array,
                       metric):
    
    point_ids, metrics = get_metric_edges(index,
                                          id_indexed_grid_df,
                                          id_array,
                                          x_array,
                                          y_array,
                                          z_array,
                                          metric)
    metrics_weights = weights_from_metrics(point_ids, metrics)

    return metrics_weights


def get_metric_edges(index,
                     id_indexed_grid_df,
                     id_array,
                     x_array,
                     y_array,
                     z_array,
                     metric):

    '''Find metric along edges to neighbours of point defined by index.

    Args:
        index (tuple) [-]: id of point under consideration.
        id_indexed_grid_df (pd.DataFrame) [-]: Grid indexed by id
        ij_indexed_grid_df (pd.DataFrame) [-]: Grid indexed by i, j indices
        metric (function): metric to apply between points

    '''

    id_point = id_indexed_grid_df.loc[index]
    id_point_coords = (id_point.x, id_point.y, id_point['layer 1 start'])
    i_index = int(id_point.i)
    j_index = int(id_point.j)
    
    search_ij = [(i_index + i, j_index + j) for i in (-1, 0, 1)
                                            for j in (-1, 0, 1)
                                                        if i != 0 or j != 0]
    
    metrics = []
    point_ids = []
    
    for i, j in search_ij:
        
        if i < 0 or j < 0: continue
        
        point_id = id_array[i, j]
        
        if np.isnan(point_id): continue
            
        ij_point_coords = (x_array[i, j],
                           y_array[i, j],
                           z_array[i, j])
        
        ij_metric = metric(id_point_coords, ij_point_coords)
        
        point_ids.append(int(point_id))
        metrics.append(ij_metric)
            
    return point_ids, metrics


def weights_from_metrics(indices, metrics):

    '''Convert metric data into format required for networkx graph object.

    Args:
        indices (list) [-]: Index of points upon which the edges terminate.
        metrics (list) [-]: Metric values.

    Return:
        weight_dict (dict) [-]: Structured dictionary.

    '''

    weight_dict = {index: {'weight': metric}
                                for index, metric in zip(indices, metrics)}

    return weight_dict


def edge_length(p1, p2):

    '''Find the distance between two points.

    Args:
        p1 (tuple) [m]: Point one x, y, z coordinates.
        p1 (tuple) [m]: Point two x, y, z coordinates.

    Return:
        edge_length (float) [m]

    '''

    edge_length = distance.euclidean((p1[0], p1[1], p1[2]),
                                     (p2[0], p2[1], p2[2]))

    return edge_length


def gradient(p1, p2):

    '''Find the gradient between two points.

    Args:
        p1 (tuple) [m]: Point one x, y, z coordinates.
        p1 (tuple) [m]: Point two x, y, z coordinates.

    Return:
        angle (float) [degrees]: Gradient expressed in degrees.

    '''

    delta_x = p2[0] - p1[0]
    delta_y = p2[1] - p1[1]
    delta_z = p2[2] - p1[2]

    horizontal_distance = np.sqrt(np.square(delta_x) + np.square(delta_y))
    angle = np.arctan(delta_z / horizontal_distance)

    return abs(np.degrees(angle))


def make_grid_arrays(grid_df):
    
    grid_df = grid_df.rename(columns={'layer 1 start': 'depth'})
    
    i_dim = grid_df.i.max() + 2
    j_dim = grid_df.j.max() + 2
    
    id_array = np.ones((i_dim, j_dim)) * np.nan
    x_array = np.ones((i_dim, j_dim)) * np.nan
    y_array = np.ones((i_dim, j_dim)) * np.nan
    z_array = np.ones((i_dim, j_dim)) * np.nan
    
    for row in grid_df.itertuples():
        
        idx = row.i
        jdx = row.j
        
        id_array[idx, jdx] = row.id
        x_array[idx, jdx] = row.x
        y_array[idx, jdx] = row.y
        z_array[idx, jdx] = row.depth
        
    return id_array, x_array, y_array, z_array
