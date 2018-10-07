# -*- coding: utf-8 -*-

#    Copyright (C) 2017-2018 Mathew Topper
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

from copy import deepcopy

import numpy as np
import pytest
import networkx as nx

from dtocean_electrical.optim_codes.array_layout import (
                                                dijkstra,
                                                get_export,
                                                calculate_distance_dijkstra,
                                                calculate_saving_vector)


def test_dijkstra(graph):

    length, path = dijkstra(graph, 676, 694)
    
    assert np.isclose(length, 181.2425646506693)
    assert path == range(676, 695)
        

def test_dijkstra_target_missing(grid):

    with pytest.raises(nx.NetworkXNoPath):
        dijkstra(grid.graph, 676, 695, grid)


def test_dijkstra_nopath(grid):
    
    graph_copy = deepcopy(grid.graph)
    
    # Remove a rectangle of nodes
    nodes_to_remove = [658, 659, 660, 661, 662, 663, 664, 665, 666, 685, 693,
                       712, 713, 714, 715, 716, 717, 718, 719, 720]
    
    graph_copy.remove_nodes_from(nodes_to_remove)
    
    assert graph_copy.has_node(690)
    
    with pytest.raises(nx.NetworkXNoPath):
        dijkstra(graph_copy, 676, 690, grid)
        

def test_get_export(grid):
    
    cp_loc = (491820., 6502180)
    landing_loc = (491760., 6500320.)
    
    length, path = get_export(cp_loc, landing_loc, grid, grid.graph)

    assert length > 0.
    assert path[-1] == 153


def test_calculate_distance_dijkstra(grid):
    
    layout_grid = [[None, 334],
                   [None, 344],
                   [None, 604],
                   [None, 614]]
    
    substation_location = (491770., 6502090.)

    (distance_array,
     path_array) = calculate_distance_dijkstra(layout_grid,
                                               substation_location,
                                               grid,
                                               grid.graph)
    
    assert distance_array.shape == (5, 5)
    assert np.sum(distance_array) > 0.
    assert distance_array[1, 2] == distance_array[2, 1]
    
    assert path_array.shape == (5, 5)
    assert len(path_array[1, 2]) > 0
    assert path_array[1, 2] == path_array[2, 1][::-1]
        
              
def test_calculate_saving_vector(grid):
    
    layout_grid = [[None, 334],
                   [None, 344],
                   [None, 604],
                   [None, 614]]
    
    substation_location = (491770., 6502090.)

    (distance_array,
     path_array) = calculate_distance_dijkstra(layout_grid,
                                               substation_location,
                                               grid,
                                               grid.graph)
    
    saving_vector = calculate_saving_vector(distance_array, 4)
    
    for vector in saving_vector:
        
        assert vector[0] != vector[1]
        
        if vector[1] == 0:
            assert np.isclose(vector[2], 0)
        else:
            assert vector[2] > 0.
        