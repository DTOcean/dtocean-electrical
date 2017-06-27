
import pandas as pd

from dtocean_electrical.network.cable import ArrayCable, ExportCable
from dtocean_electrical.network.network import Network


def test_Network_make_cable_routes():
    
    network = Network(None,
                      None,
                      None,
                      None,
                      None,
                      None,
                      None,
                      None)
    
    grid_dict = {"id": [1, 2, 3,
                        4, 5, 6,
                        7, 8, 9,
                        10, 11, 12,
                        13, 14, 15,
                        16, 17, 18],
                 "x": [0, 10, 20,
                       0, 10, 20,
                       0, 10, 20,
                       0, 10, 20,
                       0, 10, 20,
                       0, 10, 20],
                 "y": [0, 0, 0,
                       10, 10, 10,
                       20, 20, 20,
                       30, 30, 30,
                       40, 40, 40,
                       50, 50, 50],
                 'layer 1 start': [0, 0, 0,
                                   10, 10, 10,
                                   20, 20, 20,
                                   30, 30, 30,
                                   40, 40, 40,
                                   50, 50, 50],
                'layer 1 type': ["hard rock", "hard rock", "hard rock",
                                 "hard rock", "hard rock", "hard rock",
                                 "hard rock", "hard rock", "hard rock",
                                 "hard rock", "hard rock", "hard rock",
                                 "hard rock", "hard rock", "hard rock",
                                 "hard rock", "hard rock", "hard rock"]
                 }
    
    grid = pd.DataFrame(grid_dict)
    all_x = grid.x
    all_y = grid.y
    
    array_cable = ArrayCable(0,
                             20.,
                             0,
                             0,
                             [14, 11, 8],
                             [0.1] * 3,
                             [True] * 3,
                             'device',
                             'collection point',
                             0,
                             0)
    
    export_cable = ExportCable(0,
                               20.,
                               0,
                               1,
                               [8, 5, 2],
                               [0.1] * 3,
                               [True] * 3,
                               'collection point',
                               0)
    
    network.array_cables = [array_cable]
    network.export_cables = [export_cable]
    
    network.make_cable_routes(grid, all_x, all_y)
    markers = network.cable_routes.marker
    
    assert len(network.cable_routes) == 6
    assert all(x <= y for x, y in zip(markers, markers[1:]))
    