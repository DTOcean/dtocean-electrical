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

.. module:: test
   :platform: Windows
   :synopsis: Input data testing.

.. moduleauthor:: Adam Collin <adam.collin@ieee.org>
.. moduleauthor:: Mathew Topper <mathew.topper@dataonlygreater.com>
"""

import logging

# Start logging
module_logger = logging.getLogger(__name__)


def check_inputs(input_data):

    '''This checks the input data against the module data structures.
    
    Args:
        input_data (object) [-]: Instance of Electrical class.
        
    Attributes:
        site_data (object) [-]: Instance of ElectricalSiteData class.
        array_data (object) [-]: Instance of ElectricalArrayData class.
        export_data (object) [-]: Instance of ElectricalExportData class.
        options (object) [-]: Instance of ConfigurationOptions class.
        database (object) [-]: Instance of ElectricalComponentDatabase class.

    Returns:
        status
        error_string
        
    '''

    errstatus = 0
    errstr = []

        # machine checks
    # check that machine technology is ok
    try:
        machine_type = input_data.array_data.machine_data.technology
        machine_type = str(machine_type)
        
        if machine_type not in ['floating', 'fixed']: 
            
            errstatus -= 1
            errstr.append("Machine type is not valid. Please specify: "
                          "\'floating\' or \'fixed\'")

    except ValueError:

        errstatus -= 1
        errstr.append('Data format must be str compatible')

    # check that machine power is in range and is float
    try:

        power = input_data.array_data.machine_data.power
        power = float(power)

        if not 3000000.0 >= power >= 100000.0:

            errstatus -= 1
            errstr.append("Machine power does not fit: "
                          "0.1 <= Poec (MW) < 3")

    except ValueError:

        errstatus -= 1
        errstr.append('data format must be float compatible')

    # check that machine connection type is ok
    try:

        connection = input_data.array_data.machine_data.connection
        connection = str(connection)

        if connection not in ['wet-mate', 'dry-mate']:

            errstatus -= 1
            errstr.append("Connection type is not valid. Please enter: "
                          "\'wet-mate\' or \'dry-mate\'")

    except ValueError:

        errstatus -= 1
        errstr.append('Data format must be str compatible')

    # check that machine voltage is float
    try:

        voltage = input_data.array_data.machine_data.voltage
        voltage = float(voltage)

    except ValueError:

        errstatus -= 1
        errstr.append('Data format must be float compatible')

#        ## array checks
#        # check that both parts of array output exist and that lists are same
#        # length
#        if len(array_output) == 2:
#            if len(array_output[0]) != len(array_output[1]):
#                errstatus -= 1
#                errstr.append('array output parts must be of equal length')
#        else:
#            errstatus -= 1
#            errstr.append('array output must have two parts')
#
#        # check that frequency of occurence for generation sums to one
#        # may be combined with the above test
#        if sum(array_output[1]) != 1.0:
#            errstatus -= 1
#            errstr.append('frequency of occ. must sum to one')

    # check that array power is in range
#    max_power_mw = 100

    array_power_mw = (input_data.array_data.n_devices *
                      input_data.array_data.machine_data.power / 1e6)

    if  array_power_mw < 0.:

        errstatus -= 1
        errstr.append("Total installed power may not be negative")
        
#    if array_power_mw > max_power_mw:
#
#        logMsg = ("Total installed power '{} MW' exceeds maximum module "
#                  "scope of {} MW").format(array_power_mw,
#                                           max_power_mw)
#        module_logger.warning(logMsg)

    # check that landing point is in the area
    # check for optional inputs: control signals, voltage limits,
    #   reactive limits

#        # check bathy area overlap
#        lease_area = set(zip(self.site_grid.all_x, self.site_grid.all_y))
#        export_area = set(zip(self.export_grid.all_x, self.export_grid.all_y))
#        if not lease_area.intersection(export_area):
#            errstatus -= 1
#            errstr.append('No intersection between lease and export area')

    ### Check user umbililcal
    if input_data.options.user_umbilical:

        # if it is not in db, this breaks
        if input_data.options.user_umbilical not in \
                input_data.database.dynamic_cable.id.tolist():
            
            errstatus -= 1
            errstr.append('Selected umbilical cable does not exist.'
                          ' Please check database key.')

        # all others are warnings
        else:
            cable = input_data.database.dynamic_cable[\
                    input_data.database.dynamic_cable.id == \
                    input_data.options.user_umbilical]

            # check voltage rating...
            if cable.v_rate.item() != \
                    input_data.array_data.machine_data.voltage:

                msg = ("Selected umbilical cable voltage rating, {} kV, is not"
                       " equal to the device voltage rating, {} kV".format(
                       cable.v_rate.item(),
                       input_data.array_data.machine_data.voltage))

                errstatus -= 1
                errstr.append(msg)

            # ...and current rating
            if cable.a_air.item() < \
                    input_data.array_data.machine_data.max_current:

                msg = ("Selected umbilical cable current rating, {} A, is less"
                       " than the current of the device at maximum power "
                       "output, {} A".format(
                       cable.a_air,
                       input_data.array_data.machine_data.max_current))

                errstatus -= 1
                errstr.append(msg)

    return errstatus, errstr
