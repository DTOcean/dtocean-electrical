# -*- coding: utf-8 -*-
"""This defines a Grid object composed of a number of Point objects.

Author: Adam Collin
Affiliation: The University of Edinburgh
Email: a.collin@ed.ac.uk

"""

import pandas as pd
from shapely.geometry import Point
from copy import deepcopy


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

    def __init__(self, all_points):

        self.n_points = len(all_points)
        self.points = []
        self.grid_pd = all_points
        self.graph_ready = dict
        self.all_ids = []
        self.all_i = []
        self.all_j = []
        self.all_x = []
        self.all_y = []
        self.graph = None
        self.soil_types = self.get_soil_types()
        self.soil_coverage = self.get_soil_coverage()
        self.jetting_graph = None
        self.ploughing_graph = None
        self.cutting_graph = None
        self.dredging_graph = None
        self.lease_boundary = None

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

        return pd.unique(self.grid_pd['layer 1 type']).tolist()

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

        all_points = \
            [GridPoint(index, row) for index, row in self.grid_pd.iterrows()]

        self.points = all_points
#        self.n_points = len(self.points)
        self.all_ids = self.grid_pd.id
        self.all_i = self.grid_pd.i
        self.all_j = self.grid_pd.j
        self.all_x = self.grid_pd.x
        self.all_y = self.grid_pd.y

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
                [point.index for point in self.points
                 if exclusion_zone.contains(point.shapely_point)]

            all_exclusions.extend(local_exclusions)

            local_exclusions = \
                [point.index for point in self.points
                 if exclusion_zone.intersects(point.shapely_point)]

            all_exclusions.extend(local_exclusions)

        return list(set(all_exclusions))

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

        Note:
            This calls _remove_edges to remove the edge from the graph object.

        '''

        constrained_edges = \
            [(key, sub_key) for key, val in all_grads.iteritems()
             for sub_key, sub_val in val.iteritems() if sub_val['weight'] > 0]

        self._remove_edges(constrained_edges)

        return

    def _remove_edges(self, constrained_edges):

        '''Remove edge between two points in the networkx graph object.

        Args:
            constrained_edges (list): List of tuples containing u and v, where
            u and v are the Point ids defining the edge(s) to be removed.

        Returns:
            none.

        '''

        self.graph.remove_edges_from(constrained_edges)

        return

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

        points = \
           self.grid_pd[-self.grid_pd['layer 1 type'].isin(soil_list)].ix[:, 'id']

#        tool_points = self.graph_filter(soil_list)
#
#        # default is to look at all three burial protection levels
#        final_graphs = self.burial_protection_index(tool_points)

        updated_graph = (
            self.remove_points_from_graph(points, deepcopy(self.graph)))

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

            updated_graph = self.remove_points_from_graph(points_to_remove,
                                                          deepcopy(self.graph))

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

    def remove_points_from_graph(self, points_list, graph):

        '''Remove points from the networkx graph object.

        Args:
            points_list (list): List of points to be removed.
            graph (object): Networkx graph object from which the points are to
                be removed from.

        Returns:
            graph (object): Networkx graph object with points removed.

        '''

        graph.remove_nodes_from(points_list)

        return graph

    def grid_as_dict(self):

        points_as_dict = \
            {(value.i, value.j):
                (idx, value.x, value.y, value['layer 1 start'])

                for idx, value in self.grid_pd.iterrows()}

        return points_as_dict


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
        self.i = data['i']
        self.j = data['j']
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
