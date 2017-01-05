# -*- coding: utf-8 -*-
"""
Created on Thu Apr 07 13:38:45 2016

@author: acollin
"""


class Cable(object):

    '''Class to collect all attributes of a cable object.

    Args:
        index (int) [-]: Cable identification number. Each cable is numbered
            sequentially from 0 to n, where n is the number of cables in the
            network.
        length (float) [m]: Cable length.
        db_key (Unknown) [-]: Reference to database object.
        marker (int) [-]: Unique network identification number.

    Attributes:
        voltage (float) [V]: The rated voltage.
        current (float) [A]: The rated current.
        r (float) [Ohm/km]: ac resistance at 90 degree.
        c (float) [uF/km]: capacitance per unit length.
        x (float) [Ohm/km]: inductive reactance per unit length.
        type_ (str) [-]: Cable type: array, export or umbilical.
        upstream_id (int) [-]: Marker of cable termination in the sea
            direction.
        downstream_id (int) [-]: Marker of cable termination in the shore
            direction.
        upstream_type (str) [-]: Type of cable termination in the sea
            direction.
        downstream_type (str) [-]: Type of cable termination in the shore
            direction.

    Note:
        The marker is associated with the network format data structure, where
        each component is given a unique marker.

    '''

    def __init__(self, index, length, db_key, marker):

        # what attributes do cables have? These are mostly db but could be used
        # for identification
        self.voltage = None
        self.current = None
        self.r = None
        self.c = None
        self.x = None

        # what attributes do we place on cables?
        self.id_ = index
        self.db_key = db_key
        self.type_ = None
        self.length = length
        self.upstream_id = None
        self.downstream_id = None
        self.upstream_type = None
        self.downstream_type = None
        self.marker = marker

    def __str__(self):

        '''Override print command to display some info.

        '''

        return ('This is ' + self.__class__.__name__ + ' ' + str(self.id_) +
                '. This cable connects from up: ' + str(self.upstream_type) +
                ' ' + str(self.upstream_id) + ' to down: ' +
                str(self.downstream_type) + ' ' + str(self.downstream_id) +
                '.\nThe cable length is: ' + str(self.length) +
                '. The cable marker is: ' + str(self.marker)
                )


class StaticCable(Cable):

    '''Static cable definition.

    Args:
        route () [-]: Cable route.

    '''

    def __init__(self, index, length, db_key, marker, route, burial,
                 split_pipe):

        super(StaticCable, self).__init__(index, length, db_key, marker)

        # what attributes do we place on static cables?
        self.route = route  # could be route class
#        self.external_protection = None  # could be route class
#        self.target_burial_depth = None  # could be route class
        self.split_pipe = split_pipe  # could be route class
        self.target_burial_depth = burial  # could be route class


class ArrayCable(StaticCable):

    def __init__(self, index, length, db_key, marker, route,  burial,
                 split_pipe, upstream_type, downstream_type, upstream_id,
                 downstream_id):

        super(ArrayCable, self).__init__(index, length, db_key, marker, route,
                                         burial, split_pipe)

        self.type_ = 'array'
        self.upstream_type = upstream_type
        self.upstream_id = upstream_id
        self.downstream_type = downstream_type
        self.downstream_id = downstream_id


class ExportCable(StaticCable):

    def __init__(self, index, length, db_key, marker, route, burial,
                 split_pipe, upstream_type, upstream_id):

        super(ExportCable, self).__init__(index, length, db_key, marker, route,
                                          burial, split_pipe)

        self.type_ = 'export'
        self.downstream_type = 'Landing point'
        self.upstream_type = upstream_type
        self.upstream_id = upstream_id


class UmbilicalCable(Cable):

    def __init__(self, index, length, db_key, marker, seabed_connection_point,
                 device, x_coordinates, z_coordinates):

        super(UmbilicalCable, self).__init__(index, length, db_key, marker)

        self.type_ = 'umbilical'
        self.seabed_termination_x = seabed_connection_point[0]
        self.seabed_termination_y = seabed_connection_point[1]
        self.seabed_termination_z = seabed_connection_point[2]
        self.device = device
        self.x_coordinates = x_coordinates
        self.z_coordinates = z_coordinates


def get_burial_depths(route, grid, target_depth):

    '''Get the target burial depths.

    Args:
        route (list) [-]: Cable route defined by grid point id.
        grid (pd) [-]: Pandas series containing only the Target burial depth
            column of the grid_pd data.

    Attributes:
        burial_depth (list) [m]: List of burial depths.

    Returns
        burial_depth

    '''

    if target_depth is not None:

        burial_depth = [target_depth]*len(route)

    else:

        burial_depth = grid[grid.index.isin(route)].tolist()

    return burial_depth


def get_split_pipes(burial_depth):

    '''Set split pipes based on burial depth.

    Args:
        burial_depth (list) [m]: List of burial depths.

    Attributes:
        split_pipe (list) [-]: List of bools indicating need or not of split
            pipe.

    Returns:
        split_pipe

    '''

    split_pipe = [False if x > 0 else True for x in burial_depth]

    return split_pipe
