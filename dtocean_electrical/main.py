# -*- coding: utf-8 -*-
"""
This module defines the main module of the DTOcean electrical subsystems.

.. module:: main
   :platform: Windows
   :synopsis: Main module of the DTOcean electrical subsystems module
   
.. moduleauthor:: Adam Collin <a.collin@ed.ac.uk>
"""

#import os, sys
#sys.path.insert(1, os.path.join(sys.path[0], '..')) # added for local testing

# Start logging
import logging
module_logger = logging.getLogger(__name__)

from output import plot_devices
from grid.grid import Grid, GridPoint
from grid.grid_processing import grid_processing
from optim_codes.optimiser import (RadialNetwork,
                                   StarNetwork,
                                   SelectInstallationTool)
import networkx as nx
from shapely.geometry import Polygon

from input_utils.utils import snap_to_grid
from input_utils.input_tests import check_inputs

class Electrical(object):

    '''Main class of the electrical design module.
    
    Attributes:

        site_data
        array_data
        export_data
        options
        database

    '''

    def __init__(self, site_data, array_data, export_data, options, database):

        '''Unpackage inputs, assign to local variable names and do some
        processing of these.

        Note:
            check_inputs(): Checks the input data.
            _fix_features_to_grid(): Finds the nearest neighbour for each
                device and the landing point.

        '''

        self.site_data = site_data
        self.array_data = array_data 
        self.export_data = export_data
        self.options = options
        self.database = database
        self.exclusion_lines = []
#        self.status, self.error_string = (0, 0)
        self.status, self.error_string = check_inputs(self)

    def _fix_features_to_grid(self):

        '''Fix devices and landing point to the nearest grid point. This
        overwrites their location for use in the analysis of this module.

        '''

        # landing point
        self.array_data.landing_point,_ = snap_to_grid(
            self.grid.grid_pd, self.array_data.landing_point[:2])

        # devices
        device_override = {}
        device_locs= []

        for device, loc in self.array_data.layout.iteritems():

            new_loc, grid_id = snap_to_grid(self.grid.grid_pd, loc[:2])

            device_override[device] = new_loc

            device_locs.append((int(device[6:]), grid_id))

        sorted_device_locs = sorted(device_locs, key=lambda x: x[0])

        self.array_data.layout = device_override
        self.array_data.layout_grid = sorted_device_locs

        return

    def run_module(self, plot=False, iterate_tools=False):

        '''Call the Electrical module routines and return the found solution.

        Args:
            plot (bool) [-]: Boolean to control plot function.

        Attributes:
            electrical_design (object) [-]: Instance of the Optimiser class.

        Returns:
            result (object) [-]: Instance of the Network class.

        '''

        module_logger.info("Begin main run...")

        if self.status < 0:

            all_errors = ". ".join(self.error_string)
            errStr = 'Some errors in input data: {}'.format(all_errors)
            raise ValueError(errStr)

        else:

            self.grid = \
                grid_processing(self.site_data, self.export_data, self.options)

            if plot == True:
                plot_devices(self.grid,
                             self.exclusion_lines,
                             self.array_data.layout,
                             self.array_data.landing_point,
                             self.array_data.device_footprint,
                             [],
                             [],
                             [],
                             []
                             )

            self._fix_features_to_grid()

            module_logger.info("Building network...")

            # call the network connection code routine based on user option
            if self.options.network_configuration[0] == 'Radial':

                electrical_design = RadialNetwork(self, 'Radial')

            elif self.options.network_configuration[0] == 'Star':

                electrical_design = StarNetwork(self, 'Star')

            else:

                errStr = ("Network type not recognised")

                raise ValueError(errStr)
                
            # set design limits
            electrical_design.set_design_limits()

            # check for user defined installation tool
            if self.options.installation_tool:

                # need to pass the desired tool in here
                try:

                    msg = ("Checking cable routes for installation"
                           " tool: {}.".format(self.options.installation_tool))
                    module_logger.info(msg)

                    result = electrical_design.run_it(
                        self.options.installation_tool)

                    best_tool = self.options.installation_tool

                    msg = ("Solution found for installation tool: {}.".format(
                           self.options.installation_tool))
                    module_logger.info(msg)

                except (nx.NetworkXNoPath, KeyError):

                    msg = ("Could not create cable routes for selected "
                           "installation tool: {}. Trying with unfiltered "
                           "seabed.".format(
                           self.options.installation_tool))

                    module_logger.warning(msg)

                    try:

                        result = electrical_design.run_it()

                        best_tool = None
                        
                        msg = ("Solution found for installation tool: "
                               "{}".format(best_tool))
                        module_logger.info(msg)

                    except nx.NetworkXNoPath:

                        errStr = ("Could not create cable routes. "
                                  "Check exclusions.")

                        raise RuntimeError(errStr)

            else:

                # get tools from compatibility matrix
                all_tools = \
                    self.options.equipment_soil_compatibility.index.tolist()

                if iterate_tools:

                    # create var to store solutions
                    multi_tool_results = \
                        self.iterate_multi_tools(electrical_design, all_tools)

                else:

                    set_tool_order = \
                        SelectInstallationTool(self.grid.soil_coverage,
                                               all_tools,
                                               self.options.installation_rates)

                    tool_order = set_tool_order.tool_order

                    multi_tool_results = \
                        self.iterate_multi_tools(electrical_design,
                                                 tool_order,
                                                 ordered=True)

                # unpackage results
                solution = multi_tool_results['flag']
                all_solutions_cost = multi_tool_results['cost']
                all_solutions_data = multi_tool_results['data']
                
                if solution:

                    # compare tool solutions and select best for return
                    clean_solutions = \
                        {k: v for k, v in all_solutions_cost.items() if v}

                    best_tool = \
                        min(clean_solutions, key=clean_solutions.get)

                    result = all_solutions_data[best_tool]

                else:

                    errStr = ("Could not create cable routes. "
                              "Check exclusions.")

                    raise RuntimeError(errStr)

            if plot == True:

                plot_devices(self.grid,
                             self.exclusion_lines,
                             self.array_data.layout,
                             self.array_data.landing_point,
                             self.array_data.device_footprint,
                             result.collection_points,
                             result.umbilical_cables,
                             result.array_cables,
                             result.export_cables
                             )

            module_logger.info("Electrical network design complete...")

        return result, best_tool

    def iterate_multi_tools(self, electrical_design, tools, ordered=False):
        
        all_solutions_cost = {key: None for key in tools}

        all_solutions_data = {key: None for key in tools}

        solution_found = False

        for tool in tools:

            msg = ("Checking cable routes for installation"
                       " tool: {}.".format(tool))
            module_logger.info(msg)

            try:

                tool_result = electrical_design.run_it(tool)
                solution_found = True

                msg = ("Solution found for installation tool: "
                       "{}.".format(tool))
                module_logger.info(msg)

                all_solutions_data[tool] = tool_result
                all_solutions_cost[tool] = tool_result.total_cost

                if ordered is True:

                    break  # Stop when solution is found

            except (nx.NetworkXNoPath, KeyError):

                msg = ("Could not find cable routes for "
                       "installation tool: {}.".format(tool))
                module_logger.warning(msg)

        final_result = {'cost': all_solutions_cost,
                        'data': all_solutions_data,
                        'flag': solution_found}

        return final_result
