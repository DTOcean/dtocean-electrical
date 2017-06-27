# -*- coding: utf-8 -*-
"""This defines a Grid object composed of a number of Point objects.

Author: Adam Collin
Affiliation: The University of Edinburgh
Email: a.collin@ed.ac.uk

"""

import logging
module_logger = logging.getLogger(__name__)

import pandas as pd
from shapely.geometry import Point, LineString


class Grid(object):

    '''Data structure for the grid. This is composed of a number of GridPoint
    objects.

    Args:
        all_points (pd.DataFrame): Pandas DataFrame of all points. columns: id,
            i, j, x, y, layer1depth, layer1type,..,..,layerNdepth, layerNtype,
            where: i and j are local indices; x, y are grid coordinates and
            layerndepth, laterntype define the start depth and the type of
            layer n (from n=1:N).

    Attributes:
        n_points (int): The number of Point objects in the grid.
        points (list): List of the Point objects.
        grid_pd (pd.DataFrame): Pandas DataFrame of all points.
        graph_ready (dict): Grid data saved as networkx compatible dictionary.
        all_ids (list): List of all point ids.
        all_i (list): List of all point i indices.
        all_j (list): List of all point j indices.
        all_x (list): List of all point x coordinates.
        all_y (list): List of all point y coordinates.
        graph (object): Networkx graph object of the grid.
        soil_types (list): All soil types in the grid area.
        jetting_graph (dict): Dictionary of Networkx graph object of the grid
            filtered for the jetting installation tool for each burial
            protection index.
        ploughing_graph (dict): Dictionary of Networkx graph object of the grid
            filtered for the ploughing installation tool for each burial
            protection index.
        cutting_graph (dict): Dictionary of Networkx graph object of the grid
            filtered for the cutting installation tool for each burial
            protection index.
        dredging_graph (dict): Dictionary of Networkx graph object of the grid
            filtered for the dredging installation tool for each burial
            protection index.

#############
# For removal
#        jetting_graph (object): Networkx graph object of the grid filtered for
#            the jetting installation tool.
#        ploughing_graph (object): Networkx graph object of the grid filtered
#            for the ploughing installation tool.
#        cutting_graph (object): Networkx graph object of the grid filtered for
#            the cutting installation tool.
#        dredging_graph (object):Networkx graph object of the grid filtered for
#            the dredging installation tool.
#############

    Returns:
        none

    '''

    def __init__(self, all_points,
                       graph,
                       lease_boundary):

        self.grid_pd = all_points
        self.graph = graph
        self.lease_boundary = lease_boundary
        self.n_points = len(all_points)
        self.points = {}
        self.all_ids = []
        self.all_x = []
        self.all_y = []
        self.soil_types = self.get_soil_types()
        self.soil_coverage = self.get_soil_coverage()
        self.jetting_graph = None
        self.ploughing_graph = None
        self.cutting_graph = None
        self.dredging_graph = None
        
        self.add_points_to_grid()
        
        return

    def __str__(self):

        '''Override the print command for this object to display some useful
        information.

        '''

        return 'This grid has ' + str(self.n_points) + ' points.'

    def get_soil_types(self):

        '''Find all soil types present in the lease area.

        Returns:
            list: Soil types.

        '''
        
        all_soil_types = pd.unique(self.grid_pd['layer 1 type']).tolist()

        return all_soil_types

    def get_soil_coverage(self):
        
        '''Soil type coverages, as percentage of area.

        Returns:
            dict: {Soil_type: coverage}

        '''

        soil_coverage = {key: None for key in self.soil_types}
        
        for soil in self.soil_types:

            number_soil_points = \
                self.grid_pd['layer 1 type'].value_counts()[soil]
                
            percent_soil_points = float(number_soil_points)/self.n_points
            
            soil_coverage[soil] = percent_soil_points

        return soil_coverage

    def add_points_to_grid(self):

        '''Add Point objects to the grid.

        Args:

        Attributes:
            all_points (list): List of the Point objects.
            n_points (int): The number of Point objects in the grid.
            all_ids (list): List of all point ids.
            all_i (list): List of all point i indices.
            all_j (list): List of all point j indices.
            all_x (list): List of all point x coordinates.
            all_y (list): List of all point y coordinates.

        Returns:
            none.

        '''

        all_points = {row.id: GridPoint(row.id, row)
                                    for _, row in self.grid_pd.iterrows()}

        self.points = all_points

        id_grid_pd = self.grid_pd.set_index("id", drop=False)
        self.all_ids = id_grid_pd.id
        self.all_x = id_grid_pd.x
        self.all_y = id_grid_pd.y

        return

    def get_exclusion_zone_points(self, exclusion_zones):

        '''Create list of points in the exclusion zone and remove these from
        the networkx graph object.

        Args:
            exclusion_zones (list): List of exclusions zones as Shapely
                objects.

        Attributes:
            all_exclusions (list): List of Point ids to be removed from the
                networkx graph object.
            local_exclusions (list): List of Point ids to be removed from the
                networkx graph object for the given exclusion zone.

        Returns:
            all_exclusions

        Note:
            This utilises Shapely functions to determine association. Could be
            improved using matplotlib contains_points?

        '''

        all_exclusions = []

        for exclusion_zone in exclusion_zones:

            local_exclusions = \
                [point.index for point in self.points.values()
                 if exclusion_zone.contains(point.shapely_point)]

            all_exclusions.extend(local_exclusions)

            local_exclusions = \
                [point.index for point in self.points.values()
                 if exclusion_zone.intersects(point.shapely_point)]

            all_exclusions.extend(local_exclusions)

        return list(set(all_exclusions))

    def remove_exclusion_zones(self, all_exclusions):
    
        '''Remove exclusion zones from the area.
    
        Args:
            all_exclusions (list) [-]:  List of Shapely Polygon objects.
            grid (object) [-]: Instance of Grid object.
    
        Attributes:
            excluded_points (list) [-]: Index of points to be removed.
    
        '''
    
        module_logger.info("Checking for exclusion zones...")
    
        if all_exclusions:
    
            excluded_points = self.get_exclusion_zone_points(all_exclusions)
            
            # Remove from graph
            self.graph.remove_nodes_from(excluded_points)
            
            # Remove from dataframe
            self.grid_pd = self.grid_pd[self.grid_pd.id != excluded_points]
            
            # Remove from points dict
            for delp in excluded_points:
                self.points.pop(delp, None)
                
            msg = ("Number of points removed in exclusion zones: {}".format(
                   len(excluded_points)))
            module_logger.info(msg)
                
        return

    def gradient_constraint(self, gradient_limit, all_grads):

        '''Apply a gradient constraint to remove edges to neighbours which
        exceed this threshold.

        Args:
            gradient_limit (float) [deg]: The gradient limit.
            all_grads (dict) [deg]: Gradient between point and all neighbours.

        Attributes:
            constrained_edges (list): List of edges which exceed the gradient
                limit.

        Returns:
            none.

        '''
        
        module_logger.info("Checking gradient constraints...")

        constrained_edges = \
            [(key, sub_key) for key, val in all_grads.iteritems()
                 for sub_key, sub_val in val.iteritems()
                     if sub_val['weight'] > gradient_limit]
            
        n_edges = len(constrained_edges)
        
        logMsg = "Removing {} constrained edges".format(n_edges)
        module_logger.info(logMsg)

        self.graph.remove_edges_from(constrained_edges)
        constrained_lines = self._make_lines(constrained_edges)

        return constrained_lines

    def check_equipment_soil_compatibility_site(self, install_matrix):

        '''Check for suitable installation equipment in the installation area.
        This finds installation equipment which can install at all points.

        Args:
            install_matrix (pd.DataFrame): Pandas DataFrame of equipment-
            installation compatibility matrix.

        Returns:
            list: List of installation techniques which can install in the
                whole area.

        '''

        return (install_matrix[install_matrix[self.soil_types].sum(axis=1) == +
                len(self.soil_types)].index.tolist()
                )

    def check_equipment_soil_compatibility(self, technique, soil_list):

        '''This creates a graph for each installation equipment. This is
        repeated for each level of the burial protection index in a two stage
        filtering process.

        Args:
            technique (str): The installation technique under consideration.
            soil_list (list): The list of soils compatible with the technique.

        Attributes:
            tool_points (list): List of Point ids where the technique can
                install.
            updated_graph (object): Networkx graph object with points removed.

        Returns:
            none.

        Note:
            This calls graph_filter and burial_protection_index.

        '''

        points = self.grid_pd[
                    -self.grid_pd['layer 1 type'].isin(soil_list)].ix[:, 'id']

#        tool_points = self.graph_filter(soil_list)
#
#        # default is to look at all three burial protection levels
#        final_graphs = self.burial_protection_index(tool_points)

        updated_graph = self.graph.copy()
        updated_graph.remove_nodes_from(points)
        
        if technique == 'Jetting':

            self.jetting_graph = updated_graph

        elif technique == 'Ploughing':

            self.ploughing_graph = updated_graph

        elif technique == 'Cutting':

            self.cutting_graph = updated_graph

        elif technique == 'Dredging':

            self.dredging_graph = updated_graph

        return

    def burial_protection_index(self, tool_points):

        '''Check the soil type against the burial protection index and find
        intersection with tool compatibility. Then update the graph by removing
        points which are not valid for the given tool-bpi combination.

        Args:
            graph (object):

        Attrbitues:
            bpi (dict): The burial protection index filter. The key is the bpi
                level, value is the list of compatible soil type(s).
            valid_points (set): Points to install at.
            points_to_remove (list): Points to be removed from the Networkx
                graph object.
            updated_graph (object): Networkx graph object updated with points
                removed.

        Returns:
            filtered_grid (dict): Networkx graph object(s) for the given
                installation tool and all bpi levels.

        Note:
            The BPI is hard coded here and not exposed to the user. This is the
            second stage of a two stage filter process.

            This calls find_points_to_remove and remove_points_from_grid.

        '''

        bpi = {'one': ['loose sand', 'medium sand', 'dense sand',
                       'very soft clay', 'soft clay', 'firm clay',
                       'stiff clay', 'hard glacial till'
                       ],
               'two': ['loose sand', 'medium sand', 'dense sand', 'firm clay',
                       'stiff clay', 'hard glacial till'
                       ],
               'three': ['dense sand']}

        filtered_grid = {}

        for bpi_level, soil_list in bpi.iteritems():

            valid_points = set(tool_points).intersection(
                set(self.graph_filter(soil_list)))

            points_to_remove = self.find_points_to_remove(valid_points)

            updated_graph = self.graph.copy()
            updated_graph.remove_nodes_from(points_to_remove)

            filtered_grid[bpi_level] = updated_graph

        return filtered_grid

    def find_points_to_remove(self, valid_points):

        '''This function determines the points to be removed by comparing the
        valid points against the whole set.

        Args:
            valid_points (set): Points to install at.

        Returns:
            points_to_remove (list): List of points to be removed.

        '''

        points_to_remove = set(self.all_ids).symmetric_difference(valid_points)

        return list(points_to_remove)

    def graph_filter(self, soil_list):

        '''This returns points which are compatible, i.e. it must be inverted
        to be removed from the graph.

        Args:
            soil_list (list): The list of soils compatible with the technique.
            points (list): List of Point ids where the technique can install.

        Returns:
            points

        '''

        # This is the inverse - returns points which are not in the list
        # points = (
        # self.grid_pd[-self.grid_pd['layer 1 type'].isin(soil_list)].ix[:,'id'])

        points = self.grid_pd[
            self.grid_pd['layer 1 type'].isin(soil_list)].ix[:, 'id']

        return points
    
    def _make_lines(self, edge_id_list):
        
        all_lines = []
        
        for edge in edge_id_list:
                        
            start_row = self.grid_pd.loc[edge[0]]
            end_row = self.grid_pd.loc[edge[1]]
            
            line = LineString([(start_row["x"], start_row["y"]),
                               (end_row["x"], end_row["y"])])
            
            all_lines.append(line)
                    
        return all_lines


class GridPoint(object):

    '''Data structure for grid point data.

    Args:
        index (int): Unique grid point id.
        i (int): i index.
        j (int): j index.
        x (int): x coordinate corresponding to the lease area.
        y (int): y coordinate corresponding to the lease area.
        layer_depths (float): layer depth.
        layer_types (string): the layer soil type.

    Attributes:
        z (type): Description.
        neighbours (list): list of neighbours ids.
        neighbours_distance (list): list of edge distance to neighbours.
        neighbours_gradients (list): list of gradients to neighbours.
        target_burial_depth (float): target burial depth based on soil type.
        shapely_point (shapely geometry point object) []: Point stored in
            shapely format.

    Returns:
        none.

    Note:
        This is currently considering one soil layer.

    '''

    def __init__(self, index, data):

        self.index = index
        self.x = data['x']
        self.y = data['y']
        self.z = data['layer 1 start']
        self.neighbours = []
        self.neighbours_distance = []
        self.neighbours_gradient = []
        self.shapely_point = Point(self.x, self.y, self.z)

    def __str__(self):

        '''Override the print command for this object to display some useful
        information.

        '''

        s = ('Point ' + str(self.index) + ': (' + str(self.x) + ', ' +
             str(self.y) + ')')

        return s
