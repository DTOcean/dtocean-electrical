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

from shapely.geometry import LineString, Polygon


def test_get_exclusion_zone_points(grid):
    
    exclusion_zone = Polygon([(491800., 6502100.),
                              (491900., 6502100.),
                              (491900., 6502200.),
                              (491800., 6502200.)])
        
    result = grid.get_exclusion_zone_points([exclusion_zone])
    
    # Point inside the exclusion zone
    # [231, 8, 15, 491850., 6502150., -35.98, 'hard rock', 0.0],
    
    assert 231 in result


def test_remove_exclusion_zones(grid):
    
    grid_copy = deepcopy(grid)
    
    exclusion_zone = Polygon([(491800., 6502100.),
                              (491900., 6502100.),
                              (491900., 6502200.),
                              (491800., 6502200.)])
        
    grid_copy.remove_exclusion_zones([exclusion_zone])
    
    # Point inside the exclusion zone
    # [231, 8, 15, 491850., 6502150., -35.98, 'hard rock', 0.0],

    assert not grid_copy.grid_pd.id.isin([231]).any()
    assert not 231 in grid_copy.points
    assert not grid_copy.graph.has_node(231)


def test_make_lines(grid):
    
    edge_id_list = [[224, 225,
                     225, 226,
                     226, 227,
                     227, 228]]
    
    all_lines = grid._make_lines(edge_id_list)
    
    for line in all_lines:
        
        assert isinstance(line, LineString)
        assert line.length > 0
