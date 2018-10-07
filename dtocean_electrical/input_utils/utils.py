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
This collection of functions perform simple processes on the input data.

.. module:: utils
   :platform: Windows
   :synopsis: Input data processing.

.. moduleauthor:: Adam Collin <adam.collin@ieee.org>
.. moduleauthor:: Mathew Topper <mathew.topper@dataonlygreater.com>
"""

import operator

import numpy as np
from scipy import spatial
from shapely.geometry import Polygon, Point


def hydro_process(device_power_per_seastate,
                  seastate_probability,
                  power_factor,
                  n_wave_period,
                  n_wave_height,
                  n_wave_direction):

    '''Structural placeholder for the hydrodynamic processing to be shifted to
    the electrical design module.

    Args:
        unknown

    Attributes:
        unknown

    Returns:
        list: histogram of array power output; val1 = Bin edge [pc of array
            installed power], val2 = frequency of occurrence [pc]

    Note:
        This function needs to be completed. It now responds to the number of
        values requested in the power factor.

    '''

    # get power outputs for assessment from power_factor
    p_out = []

    for power in power_factor:

        p_out.append(power[0])

    probability = 1./len(p_out)

    return [p_out, [probability]*len(p_out)]


def seabed_range(bathy_data):

    '''Get the min and max water depth from the given bathymetry data.

    Returns:
        min_ (float) [m]
        max_ (float) [m]

    '''

    min_ = abs(bathy_data['layer 1 start']).min()
    max_ = abs(bathy_data['layer 1 start']).max()

    return min_, max_


def device_footprints_from_coords(layout, footprint):

    '''Get the device footprint using coordinate system.

    '''

    all_exclusions = []

    for key, value in layout.iteritems():

        exclusion = []

        for point in footprint:

            new_point = tuple(map(operator.add, value[:2], point[:2]))
            exclusion.append(new_point)

        all_exclusions.append(Polygon(exclusion))

    return all_exclusions


def device_footprints_from_rad(layout, radius):

    '''Get the device footprint from radius.

    '''

    all_exclusions = []

    for key, value in layout.iteritems():

        all_exclusions.append(Point(value).buffer(radius))

    return all_exclusions


def ideal_power_quantities(seastate_occurrence,
                           n_devices,
                           device_power):

    '''Calculate array power output assuming no losses. Used for efficiency
    calculations later in module.

    Args:
        seastate_occurrence (list):

    Attributes:
        bins
        bin_edges (list):

    Returns:
        ideal_annual_yield (float)
        ideal_histogram (list)

    '''

    bins = 1.0/len(seastate_occurrence)
    bin_edges = np.linspace(bins, 1, len(seastate_occurrence))

    year_hours = 365*24
    ideal_annual_yield = 0
    ideal_histogram = []
    array_power = n_devices * device_power

    for time, power in zip(seastate_occurrence, bin_edges):

        ideal_annual_yield += time*year_hours*array_power
        ideal_histogram.append(power*array_power)

    return ideal_annual_yield, ideal_histogram


def get_bin_edges(power_factor):

    '''Get bin edges of power factor var for analysis.

    Args:
        power_factor (list): structured input of user input power factor data.

    Returns:
         (list): bin edges

    '''

    return [edge[0] for edge in power_factor]


def snap_to_grid(grid_points, point):

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

    grid = np.array(grid_points[['x', 'y']])

    new_coords = grid[spatial.KDTree(grid).query(np.array(point))[1]].tolist()

    # and add z coord
    z = (grid_points[(grid_points.x == new_coords[0]) &
                     (grid_points.y == new_coords[1])]
                    ['layer 1 start'].values[0])

    grid_id = grid_points[(grid_points.x == new_coords[0]) &
                          (grid_points.y == new_coords[1])]['id'].values[0]

    new_coords.append(z)
    new_coords = [float(i) for i in new_coords]

    return (tuple(new_coords), grid_id)


def convert_df_column_type(df, ids, data_type):

    '''Convert pandas DataFrame column data types.

    Args:
        df (pd.DataFrame) [-]: DataFrame to be updated.
        ids (list) [-]: List of columns to be updated.
        data_type () []: Data type to be applied to columns.

    Returns:
        df (pd.DataFrame) [-]: Updated DataFrame.

    '''

    df[ids] = df[ids].astype(data_type)

    return df


def get_key(item):

    '''Method from http://pythoncentral.io/how-to-sort-a-list-tuple-or-
    object-with-sorted-in-python/

    '''

    return item[0]


def set_burial_from_bpi(row):

    '''Function to code the bpi for burial depths.

    Args:
        row (pd) [-]: Row of pandas dataframe.

    Attributes:
        bpi (float) [m]: Burial depth from burial protection index.

    Returns:
        bpi

    '''

    if row['layer 1 type'] in ['hard glacial till',
                               'cemented',
                               'soft rock coral',
                               'hard rock',
                               'gravel cobble']:

        bpi = 0.0

    elif 'sand' in row['layer 1 type']:

        bpi = 0.5

    elif 'clay' in row['layer 1 type']:

        bpi = 1.0
        
    else:
        
        errStr = "Sediment type '{}' is not recognised".format(
                                                        row['layer 1 type'])
        raise ValueError(errStr)

    return bpi
