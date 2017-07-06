
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
