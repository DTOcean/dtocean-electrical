
import numpy as np
import pandas as pd

from dtocean_electrical.network.cable import (ArrayCable,
                                              ExportCable)
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


def test_Network__map_component_types():
    
    network = Network(None,
                      None,
                      None,
                      None,
                      None,
                      None,
                      None,
                      None)
    
    types = ['export',
             'array',
             'wet-mate',
             'dry-mate',
             'substation',
             'passive hub',
             'umbilical']
    
    result = network._map_component_types(types)
    
    assert result == ['export_cable',
                      'array_cable',
                      'wet_mate_connectors',
                      'dry_mate_connectors',
                      'collection_points',
                      'collection_points',
                      'dynamic_cable']


def test_Network_calculate_annual_yield():
    
    network = Network(None,
                      None,
                      None,
                      None,
                      None,
                      None,
                      None,
                      None)
    
    network.power_histogram = [0.5, 0.5]
    network.array_power_output = [2, 0.5]
    
    annual_yield = network.calculate_annual_yield()

    assert annual_yield == 8760000000.0 + 8760000000.0 / 4


def test_Network_calculate_annual_yield_zero():
    
    network = Network(None,
                      None,
                      None,
                      None,
                      None,
                      None,
                      None,
                      None)
    
    network.power_histogram = [0.5, 0.5]
    network.array_power_output = [2, np.nan]
    
    annual_yield = network.calculate_annual_yield()

    assert annual_yield == 0.
    
    
def test_Network_calculate_lcoe():

    network = Network(None,
                      None,
                      None,
                      None,
                      None,
                      None,
                      None,
                      None)
    
    network.total_cost = 10
    network.annual_yield = 2
    
    network.calculate_lcoe()
    
    assert network.lcoe == 5e3
    
    
def test_Network_calculate_lcoe_inf():

    network = Network(None,
                      None,
                      None,
                      None,
                      None,
                      None,
                      None,
                      None)
    
    network.total_cost = 10
    network.annual_yield = 0
    
    network.calculate_lcoe()
    
    assert network.lcoe == np.inf
