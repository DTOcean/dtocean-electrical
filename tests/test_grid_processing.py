
import pytest
import numpy as np
import pandas as pd
import networkx as nx

from dtocean_electrical.grid.grid import Grid
from dtocean_electrical.grid.grid_processing import (clip_grid,
                                                     edge_length,
                                                     gradient,
                                                     get_metric_edges,
                                                     weights_from_metrics,
                                                     get_metric_weights,
                                                     get_neighbours_distance,
                                                     add_overlap_edges,
                                                     get_neighbours_gradient,
                                                     make_graph,
                                                     make_gradient_dict,
                                                     make_grid)
    

def test_clip_grid(lease, export):
    
    grid_df_clipped, static_poly = clip_grid(lease, export)
    
    assert len(grid_df_clipped) < len(export)

    assert set(export.id).intersection(set(lease.id))
    assert not set(grid_df_clipped.id).intersection(set(lease.id))
    
    lease_coords = [(x, y) for x, y in zip(lease.x, lease.y)]
    result_coords = [(x, y) for x, y in zip(grid_df_clipped.x,
                                            grid_df_clipped.y)]
        
    shared_coords = list(set(result_coords).intersection(set(lease_coords)))
    
    assert len(shared_coords) > 0
    assert static_poly.area > 0


def test_clip_grid_fail(lease, export):
    
    new_lease = lease.copy()
    new_lease['x'] = new_lease['x'].apply(lambda x: x - 100000)
    
    with pytest.raises(RuntimeError):
        clip_grid(new_lease, export)


@pytest.mark.parametrize("p2, expected", [
    ((1, 1, 1), np.sqrt(3)),
    ((1, 1, -1), np.sqrt(3))]
    )
def test_edge_length(p2, expected):
    
    p1 = (0, 0, 0)
    result = edge_length(p1, p2)
    
    assert np.isclose(result, expected)
    

@pytest.mark.parametrize("p2, expected", [
    ((0, 1, 1), 45),
    ((0, 1, -1), 45)]
    )
def test_gradient(p2, expected):
    
    p1 = (0, 0, 0)
    result = gradient(p1, p2)
    
    assert np.isclose(result, expected)


def test_get_metric_edges():
    
    grid_dict = {"id": [0, 1, 2, 3, 4, 5, 6, 7],
                 "i": [0, 1, 2, 0, 1, 2, 0, 1],
                 "j": [2, 2, 2, 1, 1, 1, 0, 0],
                 "x": [1., 2., 3., 1., 2., 3., 1., 2.],
                 "y": [6., 6., 6., 5., 5., 5., 4., 4.],
                 'layer 1 start': [0., 0., 0., 0., 0., 0., 0., 0.]}
    
    grid_df = pd.DataFrame(grid_dict)
    point_ids, metrics = get_metric_edges(4, grid_df, edge_length)
    
    assert len(point_ids) == len(metrics)
    assert set(point_ids) == set([0, 1, 2, 3, 5, 6, 7])
    assert set(metrics) == set([1.0, np.sqrt(2)])


def test_weights_from_metrics():
    
    indices = [0, 2, 6, 1, 3, 5, 7]
    metrics = [np.sqrt(2), np.sqrt(2), np.sqrt(2), 1.0, 1.0, 1.0, 1.0]

    result = weights_from_metrics(indices, metrics)
    
    assert set(result.keys()) == set(indices)
    
    for index, weight in zip(indices, metrics):
        assert np.isclose(result[index]['weight'], weight)
        
        
def test_get_metric_weights():
    
    grid_dict = {"id": [0, 1, 2, 3, 4, 5, 6, 7],
                 "i": [0, 1, 2, 0, 1, 2, 0, 1],
                 "j": [2, 2, 2, 1, 1, 1, 0, 0],
                 "x": [1., 2., 3., 1., 2., 3., 1., 2.],
                 "y": [6., 6., 6., 5., 5., 5., 4., 4.],
                 'layer 1 start': [0., 0., 0., 0., 0., 0., 0., 0.]}
    
    grid_df = pd.DataFrame(grid_dict)
    result = get_metric_weights(4, grid_df, edge_length)
    
    indices = [0, 2, 6, 1, 3, 5, 7]
    metrics = [np.sqrt(2), np.sqrt(2), np.sqrt(2), 1.0, 1.0, 1.0, 1.0]
    
    assert set(result.keys()) == set(indices)
    
    for index, weight in zip(indices, metrics):
        assert np.isclose(result[index]['weight'], weight)


def test_get_neighbours_distance():
    
    grid_dict = {"id": [0, 1, 2, 3, 4, 5, 6, 7],
                 "i": [0, 1, 2, 0, 1, 2, 0, 1],
                 "j": [2, 2, 2, 1, 1, 1, 0, 0],
                 "x": [1., 2., 3., 1., 2., 3., 1., 2.],
                 "y": [6., 6., 6., 5., 5., 5., 4., 4.],
                 'layer 1 start': [0., 0., 0., 0., 0., 0., 0., 0.]}
    
    grid_df = pd.DataFrame(grid_dict)
    result = get_neighbours_distance(grid_df)
    
    assert 4 in result
    
    indices = [0, 2, 6, 1, 3, 5, 7]
    metrics = [np.sqrt(2), np.sqrt(2), np.sqrt(2), 1.0, 1.0, 1.0, 1.0]
    
    assert set(result[4].keys()) == set(indices)
    
    for index, weight in zip(indices, metrics):
        assert np.isclose(result[4][index]['weight'], weight)
        
    assert 0 in result
    
    indices = [4, 1, 3]
    metrics = [np.sqrt(2), 1.0, 1.0]
    
    assert set(result[0].keys()) == set(indices)
    
    for index, weight in zip(indices, metrics):
        assert np.isclose(result[0][index]['weight'], weight)


def test_add_overlap_edges():

    grid_dict = {"id": [0, 1, 2, 3, 4, 5, 6, 7],
                 "i": [0, 1, 2, 0, 1, 2, 0, 1],
                 "j": [2, 2, 2, 1, 1, 1, 0, 0],
                 "x": [1., 2., 3., 1., 2., 3., 1., 2.],
                 "y": [6., 6., 6., 5., 5., 5., 4., 4.],
                 'layer 1 start': [0., 0., 0., 0., 0., 0., 0., 0.]}
    overlap_dict = {"id_x": [0],
                    "id_y": [100]}
    
    grid_df = pd.DataFrame(grid_dict)
    grid_dict = get_neighbours_distance(grid_df)
    overlap_df = pd.DataFrame(overlap_dict)
    
    pre_length = len(grid_dict[0])
    
    add_overlap_edges(grid_dict, overlap_df)
    
    assert len(grid_dict[0]) == pre_length + 1
    assert 100 in grid_dict[0].keys()


def test_get_neighbours_gradient():
    
    grid_dict = {"id": [0, 1, 2, 3, 4, 5, 6, 7],
                 "i": [0, 1, 2, 0, 1, 2, 0, 1],
                 "j": [2, 2, 2, 1, 1, 1, 0, 0],
                 "x": [1., 2., 3., 1., 2., 3., 1., 2.],
                 "y": [6., 6., 6., 5., 5., 5., 4., 4.],
                 'layer 1 start': [0., 0., 0., 0., 1., 0., 0., 0.]}
    
    grid_df = pd.DataFrame(grid_dict)
    result = get_neighbours_gradient(grid_df)
    
    angle_rad = np.arctan(1. / np.sqrt(2))
    angle_deg = np.degrees(angle_rad)

    indices = [0, 2, 6, 1, 3, 5, 7]
    metrics = [angle_deg, angle_deg, angle_deg, 45., 45., 45., 45.]
    
    assert set(result[4].keys()) == set(indices)
    
    for index, weight in zip(indices, metrics):
        assert np.isclose(result[4][index]['weight'], weight)

    
def test_make_graph(lease, export):
    
    clipped, _ = clip_grid(lease, export)
    graph = make_graph(lease, clipped)
    
    assert isinstance(graph, nx.Graph)
    
    distance = nx.dijkstra_path_length(graph, 500, 2000)

    assert distance > 0.


def test_make_gradient_dict(lease, export):
    
    clipped, _ = clip_grid(lease, export)
    graph_dict = make_gradient_dict(lease, clipped)
    
    assert set(graph_dict.keys()) == set(lease.id) | set(clipped.id)


def test_make_grid(lease, export):

    clipped, lease_polygon = clip_grid(lease, export)
    network_graph = make_graph(lease, clipped)
    
    grid = make_grid(lease,
                     clipped,
                     network_graph,
                     lease_polygon)
    
    assert isinstance(grid, Grid)
