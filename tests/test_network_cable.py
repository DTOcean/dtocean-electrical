
from dtocean_electrical.network.cable import get_burial_depths
    
    
def test_get_burial_depths_target(grid):
    
    burial_depth = get_burial_depths([36, 37], grid.grid_pd, 10)
        
    assert burial_depth == [10, 10]
    
    
def test_get_burial_depths(grid):
    
    burial_depth = get_burial_depths([36, 37], grid.grid_pd)
    
    assert burial_depth == [0., 0.]
