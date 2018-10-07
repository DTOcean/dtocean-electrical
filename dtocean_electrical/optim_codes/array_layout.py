# -*- coding: utf-8 -*-

#    Copyright (C) 2016 Adam Collin
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

"""
This module defines the DTOcean electrical subsystems array routing functions.

.. module:: array_layout
   :platform: Windows
   :synopsis: Intra-array cable routing functions.
   
.. moduleauthor:: Adam Collin <adam.collin@ieee.org>
.. moduleauthor:: Mathew Topper <mathew.topper@dataonlygreater.com>
"""

import logging
import itertools

import numpy as np
import networkx as nx
import scipy.spatial
from scipy import spatial
from scipy.misc import comb
from shapely.geometry import LineString, Point, MultiPoint


module_logger = logging.getLogger(__name__)


def snap_to_grid(grid, point, lease):
    
        '''Snap a point to the grid.
        
        Args:
            grid (np.array) [m]: Array of x and y coordinates.
            point (tuple) [m]: Coordinates of point under consideration, x and
                y coordinates.
            
        Attributes:
            new_coords (list) [m]: Coordinates of nearest point, x, y and z.
            
        Returns:
            tuple

        '''

        new_coords = grid[spatial.KDTree(grid).query(
            np.array(point))[1]].tolist()
        
        # and add z coord
        z = lease[(lease.x == new_coords[0]) &
                    (lease.y == new_coords[1])]\
                    ['layer 1 start'].values[0]
    
        new_coords.append(z)
        new_coords = [float(i) for i in new_coords]
        
        return tuple(new_coords)

def set_substation_to_edge(line, lease_area_ring, lease_bathymetry, lease):
    
    # set to lease edge
    # find poi between line of intial to shore and lease area

    # can make function
    
    lease_x = lease_bathymetry.x.tolist()
    lease_y = lease_bathymetry.y.tolist()

    grid_to_search = np.array([lease_x, lease_y]).T

    if line.intersects(lease_area_ring):
        
        poi = lease_area_ring.intersection(line)
        poi = [poi.x, poi.y]

        # snap to nearest point
        interim_estimate_snapped = snap_to_grid(
            grid_to_search, poi, lease_bathymetry)

        # then shift
        
        # find cp_loc in list of points
        grid_point = lease_bathymetry[
            (lease_bathymetry.x == interim_estimate_snapped[0]) & 
            (lease_bathymetry.y == interim_estimate_snapped[1])]
                
        # get neighbours
        check_x_direction = [0, 0, -1, +1, -1, +1, -1, +1]
        check_y_direction = [-1, +1, 0, 0, -1, +1, -1, +1]
        
        neighbour_ids = [
            (grid_point.i.item() - i_shift, grid_point.j.item() - j_shift)
            for i_shift, j_shift
            in zip(check_x_direction, check_y_direction)]
        
        for neighbour in neighbour_ids:
            
            neighbour = lease_bathymetry[
                            (lease_bathymetry.i == neighbour[0]) &
                            (lease_bathymetry.j == neighbour[1])]

            neighbour_shapely = Point(neighbour.x, neighbour.y)

            if neighbour_shapely.within(lease):

                interim_estimate = (neighbour.x.item(), neighbour.y.item())

                break

    return interim_estimate

def substation_in_site(grid_df, cp_loc, lease):
    
    # find cp_loc in list of points
    grid_point = grid_df[(grid_df.x == cp_loc[0]) & (grid_df.y == cp_loc[1])]
            
    # get neighbours
    check_x_direction = [0, 0, -1, +1, -1, +1, -1, +1]
    check_y_direction = [-1, +1, 0, 0, -1, +1, -1, +1]
    
    neighbour_ids = [
        (grid_point.i.item() - i_shift, grid_point.j.item() - j_shift)
        for i_shift, j_shift
        in zip(check_x_direction, check_y_direction)]
    
    for neighbour in neighbour_ids:
        
        neighbour = grid_df[
                        (grid_df.i == neighbour[0]) &
                        (grid_df.j == neighbour[1])]
                        
        neighbour_shapely = Point(neighbour.x, neighbour.y)
        
        if neighbour_shapely.within(lease):
            
            cp_estimate = (neighbour.x.item(), neighbour.y.item())

            break
        
    return cp_estimate

def closeness_test(device_points, cp_loc, threshold):
    
    close = False
            
    for oec in device_points:

        if oec.distance(Point(cp_loc[:2])) < threshold:

            close = True
            
    return close

def offset_cp(device_loc, export, cp_loc_estimate, position, distance):

    # make points
    device_points = []

    for item in device_loc:
        device_points.append(Point(item[0], item[1]))

    point_collection = MultiPoint(device_points)

    point_x, point_y = zip(*[(point[0], point[1]) for
                        point in
                        zip(*point_collection.envelope.exterior.xy)[:-1]])

    edges = zip(*point_collection.envelope.exterior.xy)
    edge_strings = []
    
    for id_, point in enumerate(edges[:-1]):
        # make line
        edge_strings.append(LineString([point, edges[id_+1]]))
    
    for item in edge_strings:

        poi = export.intersection(item)

        if poi:

            area_centre = (item.centroid.xy[0][0], item.centroid.xy[1][0])

            if position == 'edge':

                new_cp_loc = Point(area_centre[0], area_centre[1])

            elif position == 'beyond':

                new_cp_loc = offset_cp_local(cp_loc_estimate[0],
                                             area_centre,
                                             distance)

            elif position == 'external':
                
                # shift cp beyond lease area into export cable corridor
                # not a valid solution for wp4
                new_cp_loc = offset_cp_local(cp_loc_estimate[0],
                                             area_centre,
                                             distance)

            else:

                # add warning that substation location not found
                pass

    return (new_cp_loc.x, new_cp_loc.y)
        
def offset_cp_local(cp_loc, array_edge, distance):
    
    '''Offset the collection point beyond the array boundary by value specified
    in distance.

    Args:
        cp_loc
        array_edge

    '''

    extended_line = extend_line(cp_loc, array_edge)
    # make line from edge to end of extended line
    # first_point = (array_edge[0], array_edge[1])
    last_point = (extended_line.xy[0][1], extended_line.xy[1][1])
    new_end = LineString([array_edge, last_point]).interpolate(distance)

    return new_end

def extend_line(p1, p2):
    
    '''Extend line in p1 -> p2 direction.

    http://stackoverflow.com/questions/33159833/shapely-extending-line-feature

    Args:
        p1
        p2

    '''

    ratio = 5
    a = p1
    b = (p1[0]+ratio*(p2[0]-p1[0]), p1[1]+ratio*(p2[1]-p1[1]))
    
    return LineString([a,b])


def dijkstra(graph, a, b, grid=None):
    
    def log_error():
        
        id_grid_pd = grid.grid_pd.set_index("id")
        target_missing = b not in id_grid_pd.index
        
        if target_missing:
            
            logMsg = "Target node '{}' is not in grid".format(b)
        
        else:
        
            source_node = id_grid_pd.loc[a]
            target_node = id_grid_pd.loc[b]
            logMsg = ("Point ({}, {}) not reachable from point "
                      "({}, {}), with node numbers {} and {}").format(
                                                          target_node["x"],
                                                          target_node["y"],
                                                          source_node["x"],
                                                          source_node["y"],
                                                          b,
                                                          a)
                    
        module_logger.debug(logMsg)
        
        return
    
    if not graph.has_node(a):
        raise nx.NetworkXNoPath("node {} not in graph".format(a))

    try:
        length, path = nx.single_source_dijkstra(graph, a, b)
    except nx.NetworkXNoPath:
        if grid is not None: log_error()
        raise nx.NetworkXNoPath("node {} not reachable from {}".format(b, a))

    return length, path


def get_export(cp_loc, landing_loc, grid, graph):
    
    '''Get the export cable route and length.
    '''
    
    a = get_single_location(landing_loc[:2], grid.grid_pd)
    b = get_single_location(cp_loc[:2], grid.grid_pd)
    
    length, route = dijkstra(graph, a, b, grid)
    
    return length, route


def calculate_distance(array_layout, n_oec, substation_location):
    
    '''Calculate the distance between all devices in array_layout and between
    all devices and a fixed point defined by substation_location.

    Args:
        array_layout (dict) [m]: locations of each oec in the array.
        n_oec (int) [-]: the number of oecs in the array.
        substation_location (float) [-]: substation_location as x,y,z.

    Attributes:
        distance_list (list) [m]: list of device-device and device-substation
            spacings.
        distance_array (np.array) [m]: structured array of device-device and
            device-substation spacings.
            
    Returns:
        distance_array

    ''' 

    keylist = []

    layoutlist = []

    for key in array_layout:
        keylist.append(int(key[6:]))

    for i in sorted(keylist):
        layoutlist.append(array_layout.get(str('Device00')+str(i)))

    layoutlist.insert(0,(substation_location[0],substation_location[1],0.0))

    distance_list = ([scipy.spatial.distance.euclidean(a,b)
                      for a, b
                      in itertools.product(np.asarray(layoutlist),repeat=2)])

    distance_array = np.array(distance_list).reshape(n_oec+1, n_oec+1)

    return distance_array


def calculate_distance_dijkstra(layout_grid,
                                substation_location,
                                grid,
                                graph):
                                    
    '''Calculate the distance between all devices in layout_grid and between
    all devices and a fixed point defined by substation_location. This uses
    dijkstras algorithm to calculate the seabed distance.

    Args:
        layout_grid
        substation_location (float) [-]: substation_location as x,y,z.
        grid (object) [-]: grid object of seabed.
        graph
            
    Returns:
        distance_array
        path_array
        
    '''
        
    layoutlist = [x[1] for x in layout_grid]
    
    # Insert substation location
    substation_location = get_single_location(substation_location[:2],
                                              grid.grid_pd)
    layoutlist.insert(0, substation_location)
    
    # Initialise output arrays
    nlocs = len(layoutlist)
    distance_array = np.zeros([nlocs, nlocs])
    path_array = np.empty([nlocs, nlocs], dtype=object)
    
    ncombos = int(comb(nlocs, 2))
    module_logger.debug("{} node combinations found".format(ncombos))
    
    layout_ids = range(nlocs)
    
    for i, j in itertools.combinations(layout_ids, 2):
        
        module_logger.debug("Evaluating node combination "
                            "({}, {})".format(i, j))
        
        a = layoutlist[i]
        b = layoutlist[j]
    
        ij_length, ij_path = dijkstra(graph, a, b, grid)
        
        distance_array[i, j] = ij_length
        path_array[i, j] = ij_path
                  
        # Reverse paths
        distance_array[j, i] = ij_length
        path_array[j, i] = ij_path[::-1]
    
    return distance_array, path_array


def get_single_location(point, grid_points):
    
    '''Get the location of a single component on the grid.

    '''
    
    grid = np.array(grid_points[['x','y']])
    
    new_coords = grid[spatial.KDTree(grid).query(np.array(point))[1]].tolist()
    
    # and add z coord
    idx = grid_points[(grid_points.x == new_coords[0]) & 
            (grid_points.y == new_coords[1])]['id'].values[0]

#    new_coords.append(z)
#    new_coords = [float(i) for i in new_coords]
    
    return idx

def create_new_for_analysis(layout, layout_grid):
    
    '''Map the devices to the grid.
    
    Note:
        This assumes direct overlapping. Will this always be the case?
        Unlikely but check this.

    '''

    new = []    
    
    for item in layout_grid:
        
        for key, value in layout.iteritems():
            
            if key == 'Device' + str(item[0]).zfill(3):
                
                new.append((item[0],item[1],value[:2]))
                
    return new
    
def get_device_locations(layout, site_grid):
    
    '''Map the devices to the grid.
    
    Note:
        This assumes direct overlapping. Will this always be the case?
        Unlikely but check this.

    '''
    
    device_on_grid = []

    for key, value in layout.iteritems():
        for point in site_grid.points.values():
            if np.isclose(point.x, value[0]) and np.isclose(point.y, value[1]):
                device_on_grid.append((int(key[6:]),
                                       point.index,
                                       (value[0], value[1]))
                                       )

    device_on_grid = sorted(device_on_grid, key=lambda x: x[0])
    
    return device_on_grid


def calculate_saving_vector(distance_vector, n_oec):
    
    '''Description.
    
    Args:
        args (type): Description.
    
    Attributes:
        attributes (type): Description.
    
    Returns:
        saving_vector_filtered (list): Description.

    '''
    
    saving_vector = [(i, j, distance_vector[0][i] - distance_vector[j][i])
                        for i in range(0, n_oec + 1)
                            for j in range(0, n_oec + 1) if i != j]

    saving_vector_sorted = sorted(saving_vector,
                                  key=lambda i: (float(i[2])), reverse=True)
    # tidy up savings vector by removing self connecting nodes
    saving_vector_filtered = []
    for point in saving_vector_sorted:
        if point[0] != point[1] and point[2] >= 0.0:
            saving_vector_filtered.append(point)

    return saving_vector_filtered


def check_in_route(point_to_check, route):
    
    '''Description.
    
    Args:
        args (type): Description.
    
    Attributes:
        attributes (type): Description.
    
    Returns:
        returns (type): Description.

    '''
    
    if (point_to_check[0],0) in route:
        in_set = True
    else:
        in_set = False
    return in_set
    
def check_in_path(point_to_check, path):
    
    '''Description.
    
    Args:
        args (type): Description.
    
    Attributes:
        attributes (type): Description.
    
    Returns:
        returns (type): Description.

    '''
    
    same_path = []
    for i in range(0, len(path)):
        if point_to_check[0] in path[i] and point_to_check[1] in path[i]:
            same_path.append(True)
        else:
            same_path.append(False)
    if not sum(same_path):
        same_path = False
    else:
        same_path = True
    return same_path

def check_neighbour_number(route, point_to_check):
    
    '''Description.
    
    Args:
        args (type): Description.
    
    Attributes:
        attributes (type): Description.
    
    Returns:
        returns (type): Description.

    '''
    
    counter_u = 0
    one_neighbour = []
    for i in range(0,len(route)):
        arc = route[i]

        if arc[0] == point_to_check[1]:
            counter_u+=1
        elif arc[1] == point_to_check[1]:
            counter_u+=1
    
        if counter_u > 1:
            one_neighbour.append(False)
        else:
            one_neighbour.append(True)
    
    if False in one_neighbour:
        t = False
    else:
        t = True
    return t

def check_path_capacity(path, point_to_check, cap):
    
    '''Description.
    
    Args:
        args (type): Description.
    
    Attributes:
        attributes (type): Description.
    
    Returns:
        returns (type): Description.

    '''
    
    for i in path:
        if point_to_check[0] in i:
            # get length
            path_length_k = len(i)-1
        if point_to_check[1] in i:
            # get length
            path_length_u = len(i)-1
    
    # sum path lengths and check capacity
    total = path_length_k + path_length_u
    if total > cap:
        cap_exceed = True
    else:
        cap_exceed = False
    return cap_exceed

def crossing_dijkstra(path, route, path_array, devices, site_grid):
    
    '''Controller for checking cable crossing using dijkstras algorithm
    generated paths.
    
    '''

    # make line
    line1 = path_array[path[0]][path[1]]
    line1 = make_linestring(line1, site_grid)
    
    # initiate crossing vector
    cross = []
    for edge in route:
        # make line
        line2 = path_array[edge[0]][edge[1]]
        line2 = make_linestring(line2, site_grid)
        if line1.intersection(line2).is_empty:
            cross.append(False)
        elif line1.intersection(line2) == Point(devices[path[0]-1][2][0],devices[path[0]-1][2][1]):
            cross.append(False)
        elif line1.intersection(line2) == Point(devices[path[1]-1][2][0],devices[path[1]-1][2][1]):
            cross.append(False)
        elif line1.distance(Point(devices[path[0]-1][2][0],devices[path[0]-1][2][1])) < 1e-3:
            cross.append(False)
        elif line1.distance(Point(devices[path[1]-1][2][0],devices[path[1]-1][2][1])) < 1e-3:
            cross.append(False)
        else:
            cross.append(True)
            # Need to check if intersection is at device
#            pass
#            cross.append(True)

    return sum(cross)

def make_linestring(path, site_grid):

    line_path = []    
    for point in path:
        line_path.append((site_grid.points[point].x,
                          site_grid.points[point].y))
    
    return LineString(line_path)

def crossing(path, layout, route):
    
    '''Description.
    
    Args:
        args (type): Description.
    
    Attributes:
        attributes (type): Description.
    
    Returns:
        returns (type): Description.

    '''
    
    pref = str('Device00')+str(path[0])
    p2ref = str('Device00')+str(path[1])
    p = np.array((layout[pref][0],layout[pref][1]))
    p2= np.array((layout[p2ref][0],layout[p2ref][1]))
    
    # initiate crossing vector
    cross = []
    for edge in route:
        # p and p2 are defined by point
        # q and q2 are defined by the edge in R
        qref = str('Device00')+str(edge[0])
        q2ref = str('Device00')+str(edge[1])
        q = np.array((layout[qref][0],layout[qref][1]))
        q2= np.array((layout[q2ref][0],layout[q2ref][1]))
        cross_check = _check_cable_crossing(p, p2, q, q2)
        cross.append(cross_check)

    return sum(cross)

def _check_cable_crossing(p, p2, q, q2):
    
    '''Check if line1 - defined by point p and point p2 - and line2 - defined
    by point q and point q2 intersect. 
    
    Args:
        p (numpy.ndarray): point defining the start of line1.
        p2 (numpy.ndarray): point defining the end of line1.
        q (numpy.ndarray): point defining the start of line2.
        q2 (numpy.ndarray): point defining the end of line2.
    
    Attributes:
        p_o_i (tuple): point of intersection.
    
    Returns:
        returns (type): Description.

    '''

    r = p2 - p
    s = q2 - q

    rxs = np.cross(r,s)
    qp = q - p
    qpxr = np.cross(qp,r)

    try:
        
        t = np.cross(qp,s)/rxs
        u = np.cross(qp,r)/rxs
        
    except RuntimeWarning:
        
        msg = ("Cable crossing")
        module_logger.info(msg)
        

    if rxs == 0 and qpxr == 0:
        crossing = False
        ## May need to improve logical testing of collinear points
    elif rxs != 0 and (0 <= t <= 1) and (0 <= u <=1):
        crossing = True
        
        if float(p2[0]-p[0]) == 0 or float(q2[0]-q[0]) == 0:
            vertical = True
        else:
            vertical = False
        # need to check if any lines are horizontal
        if float(p2[1]-p[1]) == 0 or float(q2[1]-q[1]) == 0:
            horizontal = True
        else:
            horizontal = False
        
        # based on this, calculate the point of intersection
        if vertical and not horizontal:
            if float(p2[0]-p[0]) == 0:
                # gradient            
                m2 = gradient(q, q2)
                # y-intercept
                c2 = q2[1]-m2*q2[0]
                ## point of intersection
                x_intersect = p2[0]
                y_intersect= m2*x_intersect+c2
            else:
                # gradient
                m1 = gradient(p, p2)
                # y-intercept
                c1 = p[1]-m1*p[0]
                ## point of intersection
                x_intersect = q2[0]
                y_intersect= m1*x_intersect+c1

        elif horizontal and not vertical:
            if float(p2[1]-p[1]) == 0:
                # gradient
                m2 = gradient(q, q2)
                # y-intercept
                c2 = q2[1]-m2*q2[0]
                y_intersect= p2[1]
                x_intersect = (y_intersect-c2)/m2
            else:
                # gradient
                m1 = gradient(p, p2)
                # y-intercept
                c1 = p[1]-m1*p[0]
                y_intersect= q2[1]
                x_intersect = (y_intersect-c1)/m1

        elif vertical and horizontal:
            if float(p2[1]-p[1]) == 0:
                x_intersect = q2[0]
                y_intersect = p2[1]
            else:
                x_intersect = p2[0]
                y_intersect = q2[1]

        else:
            ## equation line 1
            # gradient
            m1 = gradient(p,p2)
            # y-intercept
            c1 = y_intercept(p,m1)
            ## equation line 2
            # gradient
            m2 = gradient(q,q2)
            # y-intercept
            c2 = q2[1]-m2*q2[0]
            ## point of intersection
            x_intersect = (c2-c1)/(m1-m2)
            y_intersect= m2*x_intersect+c2
    else:
        crossing = False

    # check that crossing does not occur at start point
    if crossing:
        p_o_i = [x_intersect, y_intersect]
        if ((np.allclose(p_o_i, np.ndarray.tolist(p),1e-3)) or
            (np.allclose(p_o_i, np.ndarray.tolist(p2),1e-3)) or
            (np.allclose(p_o_i, np.ndarray.tolist(q),1e-3)) or
            (np.allclose(p_o_i, np.ndarray.tolist(q2),1e-3))):

                crossing = False

    return crossing

def gradient (point_one, point_two):
    
    '''Calculate the gradient between two points.
    
    Args:
        point_one (list)
        point_two (list)
        
    Attributes:
        m (float) [m]: the gradient!
        
    Returns:
        m

    '''
    
    m = (point_two[1]-point_one[1])/float(point_two[0]-point_one[0])
    
    return m

def y_intercept(point, gradient):
    
    '''Calculate the y intercept of a line defined by two points.
    
    Args:
        point (list)
        gradient (float)
        
    Attributes:
        c (float) [m]: the intercept!
        
    Returns:
        c

    '''
    
    c = point[1]-gradient*point[0]
    return  c

def update_path(point, route):

    '''Description.
    
    Args:
        args (type): Description.
    
    Attributes:
        attributes (type): Description.
    
    Returns:
        returns (type): Description.

    '''
    
    route.append((point[0],point[1]))
    # Remove link to 0 point
    if (point[0],0) in route:
        route.remove((point[0],0))

    path = []
    route_copy = list(route)

    path_interim = []
    start = route_copy[0][1]
    end = route_copy[0][0]
    del route_copy[0]
    path_interim.append(start)
    path_interim.append(end)

    startcheck = []

    # check if end node equals any start nodes
    while route_copy:
        for i in route_copy:
            startcheck.append(i[1])

        if end in startcheck:
            for i in route_copy:    
                if i[1] == end:
                    path_interim.append(i[0])
                    route_copy.remove(i)
                    end = i[0]
                    startcheck = []

        else:
            path.append(path_interim)
            path_interim = []
            start = route_copy[0][1]
            end = route_copy[0][0]
            del route_copy[0]
            path_interim.append(start)
            path_interim.append(end)
            startcheck = []

    # start new list
    path.append(path_interim)
    
    return path

def run_this(saving_vector, path, route, n_oec, string_max, layout):
    
    '''Description.
    
    Args:
        args (type): Description.
    
    Attributes:
        attributes (type): Description.
    
    Returns:
        returns (type): Description.

    '''
    
    for point in saving_vector:
        if ((check_in_path(point, path) == 0) and
            (check_in_route(point, route) == 1) and
            (check_neighbour_number(route, point) == 1) and
            (check_path_capacity(path, point, string_max) == 0) and
            (crossing(point, layout, route) == 0)):
                path = update_path(point, route)

    return path

def run_this_dijkstra(saving_vector, path, route, n_oec, string_max,
                      path_array, devices, site_grid):
    
    '''Description.
    
    Args:
        args (type): Description.
    
    Attributes:
        attributes (type): Description.
    
    Returns:
        returns (type): Description.

    '''
    
    for point in saving_vector:
        if ((check_in_path(point, path) == 0) and
            (check_in_route(point, route) == 1) and
            (check_neighbour_number(route, point) == 1) and
            (check_path_capacity(path, point, string_max) == 0) and
            (crossing_dijkstra(point,
                               route, 
                               path_array,
                               devices,
                               site_grid) == 0)):

                path = update_path(point, route)

    return path


