# -*- coding: utf-8 -*-
"""
Created on Thu Apr 07 13:32:59 2016

@author: acollin
"""

import numpy as np
import pandas as pd
from copy import deepcopy
from collections import Counter
from collection_point import PassiveHub, Substation
from cable import (ArrayCable, ExportCable, UmbilicalCable, get_burial_depths,
                   get_split_pipes)
from connector import WetMateConnector, DryMateConnector
from pprint import pformat
import logging

# Start logging
module_logger = logging.getLogger(__name__)


class Network(object):

    '''Data structure for the network description. This is composed of a number
    of network component objects and contains all data required to describe the
    network structure and performance.

    Args:
        index (int) [-]: Unique index of the instance of the Network object.
        configuration (str) [-]: The network configuration, defined as either:
            radial or star.
        power_histogram (list) [pc]: The probability of occurrence of each
            power bin.
        array_power_output (list) [MW]: The array power output for each for
            each bin edge in the power histogram.

    Attributes:
        index (int)
        configuration (str)
        n_cp (int) [-]: The number of collection points.
        n_export (int) [-]: The number of export cables.
        export_voltage (float) [kV]: The export cable voltage.
        array_voltage (float) [kV]: The array cable voltage.
        shore_to_device (nparray) [-]: Shore to device connection matrix.
        device_to_device (nparray) [-]: Device to device connection matrix.
        shore_to_cp (nparray) [-]: Shore to collection point connection matrix.
        cp_to_cp (nparray) [-]: Collection point to collection point connection
            matrix.
        cp_to_device (nparray) [-]: Collection point to device connection
            matrix.
        export_cables (list) [-]: List of ExportCable objects.
        array_cables (list) [-]: List of ArrayCable objects.
        collection_points (list) [-]: List of CollectionPoint objects.
        wet_mate (list) [-]: List of WetMateConnector objects.
        dry_mate (list) [-]: List of DryMateConnector objects.
        seastate_occurrence () [-]
        array_power_output () [-]
        b_o_m (pd.DataFrame) [-]: Network bill of materials;
            db ref (int) [-]: Component database key.
            install_type (str) [-]: Component type.
            marker (int) [-]: Component unique marker.
            quantity (float): Unit quantity of the component. Cable lengths
                given in [m], all others are dimensionless.
            utm_x (float) [m]: UTM coordinate in the x direction.
            utm_y (float) [m]: UTM coordinate in the y direction.
        economics_data (pd.DataFrame) [-]: Total component cost of the network;
            db ref (int) [-]: Component database key.
            quantity (float) [-]: Total unit quantity of all uses of component
                in the network. Cable lengths given in [m], all others are
                dimensionless.
            cost (float) [E]: Total cost of all uses of component in the
                network.
            year (int) [-]: Installation year. Assumed zero for all.
        total_cost (float) [E]: Total cost, i.e. sum, of all components in the
            network.
        all_connections (dict) [-]: Structure which carries all data required
            for hierarchy and network_design. Each component is represented by
            the database id and the unique component marker.
        hierarchy (dict): Structure to carry the component-to-component
            connection relationship. See notes for further information.
        network_design (dict) [-]: Network design in the dictionary format
            required for downstream analysis. See notes for further
            information.
        cable_routes (pd.DataFrame): Cable route information as points;
            db ref (int) [-]: Component database key.
            marker (int) [-]: Component unique marker.
            grid_id (int) [-]: Point object id traversed by cable.
            burial_depth (float) [m]: Burial depth at the specified point.
            split pipe (bool) [-]: Presence of split pipe protection at the
                specified point.
        collection_points_design (pd.DataFrame) [-]: Collection point
            specification for foundation design;
            centre_of_gravity (np.array) [m]: Centre of gravity with respect to
                local coordinate system, as [x,y,z].
            dry beam area (float) [m2]: Dry beam area.
            dry frontal area (float) [m2]: Dry frontal area.
            foundation type (str) [-]: Predefined foundation type, either
                'gravity' or 'pile'.
            foundation locations (np.array) [m]: Foundation location with
                respect to the origin point, as [x,y,z].
            height (float) [m]: Unit height.
            length (float) [m]: Unit length.
            marker (int) [-]: Collection point unique marker.
            mass (float) [kg]: Total unit weight.
            orientation angle (float) [deg]: Device orientation angle.
            origin (np.array) [m]: Collection point origin point, UTM
                coordinates as [x,y].
            profile (str) [-]: Shape, either 'cylindrical' or 'rectangular'.
            surface roughness (float) [m]: Surface roughness.
            type_ (str) [-]: Collection point type, either 'subsea passive',
                'subsea substation' or 'surface substation'.
            volume (float) [m3]: Submerged volumne.
            wet beam area (float) [m2]: Wet beam area.
            wet frontal area (float) [m2]: Wet frontal area.
            width (float) [m]: Unit width.
        annual_yield (float) [Wh]: Array power output for a year period
            less electrical losses.
        annual_losses (float) [Wh]: Electrical losses for a year period.
        annual_efficiency (float) [pc]: Network efficiency for a year period.
        histogram_losses (list) [Wh]: Electrical losses at each power
            generation level.
        histogram_efficiency (list) [pc]: Electrical efficiency at each power
            generation level.

    Returns:
        none

    Note
        The structure of hierarchy is divided by two types of key. 'array'
        denotes the system level connection: 'Export cable' is the connection
        from the onshore landing point to the 'Substation'. 'layout' denotes
        the series/parallel connections of the oec: within the list, comma
        separated values represent series connections, brackets separated by a
        comma indicate a new branch. The 'device' keys indicate the connection
        between the last system and the device. All values are database ids.

        network_design follows the description above but specifies the number
        of unique componens in a given system. Each component is assigned a
        unique marker.

    '''

    def __init__(self,
                 index,
                 configuration,
                 power_histogram,
                 array_power_output,
                 floating,
                 substation,
                 export_constraints,
                 array_constraints):

        # network characteristics
        self.index = index
        self.configuration = configuration
        self.floating = floating
        self.substation = substation
        self.n_cp = 0
        self.n_export = 0
        self.export_voltage = 0.
        self.array_voltage = 0.
        self.shore_to_device = None
        self.device_to_device = None
        self.shore_to_cp = None
        self.cp_to_cp = None
        self.cp_to_device = None

        # network components
        self.export_cables = []
        self.array_cables = []
        self.umbilical_cables = []
        self.collection_points = []
        self.wet_mate = []
        self.dry_mate = []

        # assessment states
        self.power_histogram = power_histogram
        self.array_power_output = array_power_output

        self.array_constraints = array_constraints
        self.export_constraints = export_constraints

        # high level description
        self.b_o_m = None
        self.economics_data = None
        self.total_cost = None
        self.lcoe = None
        self.all_connections = None
        self.hierarchy = None
        self.network_design = None
        self.cable_routes = None
        self.collection_points_design = None
        self.umbilical_cable_design = None
        self.annual_yield = None
        self.annual_losses = None
        self.annual_efficiency = None
        self.histogram_losses = None
        self.histogram_efficiency = None

    def __str__(self):

        '''Override the print command for this object to display some useful
        information.

        '''

        return ('This ' + 'network has: ' +
                str(self.n_cp) + ' collection point(s), ' +
                str(len(self.array_cables)) + ' array cable(s) and ' +
                str(len(self.export_cables)) + ' export cable(s).')

    def get_db_keys(self, component):

        '''Get all database keys for the component class.

        Args:
            component (str): the component class under consideration.

        Returns:
            db_keys (list): list of database keys.

        Note:
            Can this be improved by using component to directly set the
            component class considered, something like str('self.' + component)
            but not a string?

        '''

        db_keys = set()

        if component == 'array_cables':
            for item in self.array_cables:
                db_keys.add(item.db_key)
        elif component == 'export_cables':
            for item in self.export_cables:
                db_keys.add(item.db_key)

        return list(db_keys)

    def cable_bom(self):

        '''Get the cable data for the bill of materials.

        Attributes:
            cable_keys (list): list of db keys.

        Returns:
            none.

        '''

        cable_keys = []
        for type_ in ['array_cables', 'export_cables']:
            cable_keys.append(self.get_db_keys(type_))

        return

    def add_collection_point(self, n_cp, cp_loc, db_key, db):

        '''Add collection point object(s) to the network object.

        Args:
            n_cp (int): number of collection points to be added.
            cp_loc (list): list of collection point locations.

        Returns:
            none.

        '''

        self.n_cp = n_cp

        data = db[db.id == db_key]

        for cp in range(n_cp):

            if data.v1.values[0] == data.v2.values[0]:

                self.collection_points.append(
                    PassiveHub(cp, cp_loc[cp], db_key, data))

            else:

                self.collection_points.append(
                    Substation(cp, cp_loc[cp], db_key, data))

        return

    def add_cables_cp_three(self, cp_device_distance,
                                  cp_cp_distance,
                                  device_connection,
                                  device_layout,
                                  cp_device_paths,
                                  cp_cp_paths,
                                  export_route,
                                  export_length,
                                  umbilical_data,
                                  components,
                                  burial_depths,
                                  burial_array,
                                  burial_export):

        '''Add cables to the network.  This also adds connectors at cable ends.
        This also produces the hierarchy and network design dictionaries.

        Args:
            distance_array (nparray): Distance between devices and central
                location.
            device_connection (str): Type of connector.
            device_layout (dict): Device locations.
            path_array (nparray): Seabed paths between devices and central
                location.
            export_route (tuple): Export cable route, defined by grid point
                ids.
            export_length (float): Export cable length.
            umbilical_data
            components
            burial_depths
            burial_array (float) [-]: User defined array cable burial depth.
                Can be None.
            burial_export (float) [-]: User defined export cable burial depth.
                Can be None.

        Attributes:
            marker (int): Unique component marker.
            device_to_device (nparray): Copy of the device to device connection
                matrix.
            visited_nodes (list): List of visisted devices.
            export_idx (int): Export cable index.
            array_idx (int): Array cable index.
            wet_mate_idx (int): Wet mate connector index.
            dry_mate_idx (int): Dry mate connector index.
            hierarchy (dict): Network connection hierarchy for downstream
                analysis.
            array (list): Container for array level information, to be stored
                in heirarchy.

        Returns:
            none.

        '''

        marker = 0
        export_idx = 0
        array_idx = 0
        wet_mate_idx = 0
        dry_mate_idx = 0
        umbilical_idx = 0

        # vars for dictionary structures
        hierarchy = {}
        array = []

        for connection in np.where(self.shore_to_cp > 0)[0]:
            cluster = {}
            cluster['layout'] = []
#                cable_length = distance_array[0][connection+1]
            db_key = components['export']

            burial = (
                get_burial_depths(export_route, burial_depths, burial_export))

            split_pipe = get_split_pipes(burial)

            self.export_cables.append(ExportCable(export_idx,
                                                  export_length,
                                                  db_key,
                                                  marker,
                                                  export_route,
                                                  burial,
                                                  split_pipe,
                                                  'collection point',
                                                  connection))

            cluster['Export cable'] = [(db_key, marker)]
            marker += 1
            # create substation - need to make sure that each is given a
            # marker
            # Then update the marker of the collection point

            # need to handle differently if there is no substation
            visited_nodes = []
            layout = []

            cp_to_cp = deepcopy(self.cp_to_cp)

            if self.substation:

                # need to add reference to export side connector for
                # installation

                export_connector = \
                    self.collection_points[connection].input_connectors

                db_key, wet_mate_idx, dry_mate_idx = (
                    self._add_connector(
                        components, export_connector, wet_mate_idx,
                        dry_mate_idx,
                        self.collection_points[connection].db_key,
                        marker,
                        (self.collection_points[connection].utm_x,
                         self.collection_points[connection].utm_y),
                        0))

                cluster['Export cable'].append((db_key, marker))
                marker += 1

                self.collection_points[connection].marker = marker
                cluster['Substation'] = [
                    (self.collection_points[connection].db_key,
                     self.collection_points[connection].marker)]

            else:

                cluster['Substation'] = ['Ideal']
                layout.append('subhub' + str(connection).zfill(3))
                hierarchy['subhub' + str(connection).zfill(3)] = {}

                self.collection_points[connection].marker = marker

                hierarchy['subhub' + str(connection).zfill(3)].update(
                    {'Elec sub-system': [],
                     'Substation': [
                         (self.collection_points[connection].db_key,
                          self.collection_points[connection].marker)]})

                visited_nodes.append(connection)

            marker += 1

            if self.configuration == 'Star':

                cp_to_cp[:, connection] = 0

                cp_layout = []

                for next_cp in np.where(cp_to_cp[connection] > 0)[0]:

                    visited_nodes.append(next_cp)
                    cp_to_cp[:, next_cp] = 0

                    # Link to previous cp
                    # need to add reference to previous cp side connector for
                    # installation - keep as ideal
                    array_connector = \
                        self.collection_points[connection].output_connectors

                    db_key, wet_mate_idx, dry_mate_idx = (
                        self._add_connector(
                            components, array_connector, wet_mate_idx,
                            dry_mate_idx,
                            self.collection_points[connection].db_key,
                            marker,
                            (self.collection_points[connection].utm_x,
                             self.collection_points[connection].utm_y),
                            0))

                    link_to_cp = []
                    link_to_cp.append((db_key, marker))
                    marker += 1

                    cable_length = cp_cp_distance[connection][next_cp]
                    route = cp_cp_paths[connection][next_cp]
                    db_key = components['array']

                    burial = (
                        get_burial_depths(
                            route, burial_depths, burial_array))
                    split_pipe = get_split_pipes(burial)

                    self.array_cables.append(ArrayCable(array_idx,
                                                        cable_length,
                                                        db_key,
                                                        marker,
                                                        route,
                                                        burial,
                                                        split_pipe,
                                                        'collection point',
                                                        'collection point',
                                                        connection,
                                                        next_cp))

                    link_to_cp.append((db_key, marker))
                    marker += 1

                    # Link at current cp
                    # need to add reference to previous cp side connector for
                    # installation - keep as ideal
                    array_connector = \
                        self.collection_points[next_cp].input_connectors

                    db_key, wet_mate_idx, dry_mate_idx = (
                        self._add_connector(
                            components, array_connector, wet_mate_idx,
                            dry_mate_idx,
                            self.collection_points[connection].db_key,
                            marker,
                            (self.collection_points[connection].utm_x,
                             self.collection_points[connection].utm_y),
                            0))

                    link_to_cp.append((db_key, marker))

                    layout.append('subhub' + str(next_cp).zfill(3))
                    hierarchy['subhub' + str(next_cp).zfill(3)] = {}

                    self.collection_points[next_cp].marker = marker + 1

                    hierarchy['subhub' + str(next_cp).zfill(3)].update(
                        {'Elec sub-system': link_to_cp,
                         'Substation': [
                            (self.collection_points[next_cp].db_key,
                             self.collection_points[next_cp].marker)]})

                    marker += 2  # +2 to account for collection point
                    array_idx += 1

                    chain = True

                    start_node = next_cp

                    while chain is True:

                        if np.any(cp_to_cp[next_cp] > 0):

                            next_cp = int(np.where(cp_to_cp[next_cp] > 0)[0])
                            visited_nodes.append(cp_to_cp)
                            cp_to_cp[:, next_cp] = 0

                            cable_length = cp_cp_distance[start_node][next_cp]
                            route = cp_cp_paths[start_node][next_cp]
                            db_key = components['array']

                            burial = (
                                    get_burial_depths(
                                        route, burial_depths, burial_array))
                            split_pipe = get_split_pipes(burial)

                            self.array_cables.append(
                                ArrayCable(array_idx, cable_length, db_key,
                                           marker, route, burial, split_pipe,
                                           'collection point',
                                           'collection point', start_node,
                                           next_cp))

                            layout.append('subhub' + str(next_cp).zfill(3))
                            hierarchy['subhub' + str(next_cp).zfill(3)] = {}

                            self.collection_points[next_cp].marker = marker + 1

                            hierarchy['subhub' + str(next_cp).zfill(3)].update(
                                {'Elec sub-system': [(db_key, marker)],
                                 'Substation': [
                                    (self.collection_points[next_cp].db_key,
                                     self.collection_points[next_cp].marker)]})

                            marker += 2  # +2 to account for collection point
                            array_idx += 1

                            start_node = next_cp

                        else:

                            chain = False
                            cp_layout.append(layout)
                            layout = []

                cluster['layout'] = cp_layout

            cp_to_device = deepcopy(self.cp_to_device)
            device_to_device = deepcopy(self.device_to_device)

            visited_nodes = []

            if self.floating:

                for idx, cp in enumerate(self.cp_to_device):

                    if self.n_cp > 1:

                        sub_hub_layout = []

                    for connection in np.where(cp > 0)[0]:

                        layout = []
                        link_to_cp = []
                        # find where this connects
                        visited_nodes.append(connection)
                        # then connect from here
                        chain = True
                        chain_step = connection
            #            device_to_device[:,connection] = 0

                        # Link to cp
                        # need to add reference to export side connector for
                        # installation - keep as ideal
                        array_connector = \
                            self.collection_points[idx].output_connectors

                        db_key, wet_mate_idx, dry_mate_idx = (
                            self._add_connector(
                                components, array_connector, wet_mate_idx,
                                dry_mate_idx,
                                self.collection_points[idx].db_key,
                                marker,
                                (self.collection_points[idx].utm_x,
                                 self.collection_points[idx].utm_y),
                                0))

                        link_to_cp = []
                        link_to_cp.append((db_key, marker))
                        marker += 1

                        if self.n_cp > 1:

                            cable_length = cp_device_distance[idx][chain_step]
                            route = cp_device_paths[idx][chain_step]

                        else:

                            cable_length = \
                                cp_device_distance[idx][chain_step+1]

                            route = cp_device_paths[idx][chain_step+1]

                        burial = get_burial_depths(route,
                                                   burial_depths,
                                                   burial_array)

                        split_pipe = get_split_pipes(burial)

                        db_key = components['array']
                        self.array_cables.append(ArrayCable(array_idx,
                                                            cable_length,
                                                            db_key,
                                                            marker,
                                                            route,
                                                            burial,
                                                            split_pipe,
                                                            'connector',
                                                            'collection point',
                                                            marker+1,
                                                            idx))

                        # add static cable to layout
                        link_to_cp.append((db_key, marker))
                        layout.append('device' + str(connection+1).zfill(3))
                        hierarchy['device' + str(connection+1).zfill(3)] = {}
#                        hierarchy['device' + str(connection+1).zfill(3)]\
#                            ['Elec sub-system'] = [(db_key, marker)]

                        marker += 1
                        array_idx += 1

                        # add connector to layout
                        location = (umbilical_data[
                            'Device' + str(connection+1).zfill(3)]
                            ['termination'])

                        db_key, wet_mate_idx, dry_mate_idx = (
                            self._add_connector(
                                components, device_connection, wet_mate_idx,
                                dry_mate_idx, db_key, marker, location,
                                connection))
                        link_to_cp.append((db_key, marker))
#                        hierarchy['device' + str(connection+1).zfill(3)]\
#                            ['Elec sub-system'].append((db_key, marker))

                        marker += 1
                        # add dynamic cable to layout
                        cable = umbilical_data['Device' +
                                               str(connection+1).zfill(3)]
                        self.umbilical_cables.append(UmbilicalCable(
                            umbilical_idx,
                            cable['length'],
                            cable['db_key'],
                            marker,
                            cable['termination'],
                            cable['device'],
                            cable['x coords'],
                            cable['z coords']))

                        link_to_cp.append((cable['db_key'], marker))
                        umbilical_idx += 1
                        marker += 1

                        # add device connector to layout
                        location = device_layout['Device' +
                                                 str(connection+1).zfill(3)]

                        db_key, wet_mate_idx, dry_mate_idx = (
                            self._add_connector(
                                components, device_connection, wet_mate_idx,
                                dry_mate_idx, db_key, marker, location,
                                connection))

                        link_to_cp.append((db_key, marker))
                        marker += 1

                        (hierarchy['device' + str(connection+1).zfill(3)]
                            ['Elec sub-system']) = link_to_cp

                        start = connection
                        cp_to_device[idx][connection] = 0

                        while chain is True:

                            elec_sub_system = []

                            if np.any(device_to_device[chain_step] > 0):

                                next_device = np.where(
                                    device_to_device[chain_step] > 0)[0]

                                # filter against visited nodes
                                for node in next_device:

                                    if node not in visited_nodes:

                                        next_device = node

                                layout.append('device' +
                                              str(next_device+1).zfill(3))

                                # add static cable between connectors
                                chain_step = int(next_device)

                                cable_length = \
                                    cp_device_distance[start+1][chain_step+1]

                                route = cp_device_paths[start+1][chain_step+1]

                                burial = (
                                    get_burial_depths(
                                        route, burial_depths, burial_array))

                                split_pipe = get_split_pipes(burial)

                                db_key = components['array']
                                self.array_cables.append(ArrayCable(
                                    array_idx,
                                    cable_length,
                                    db_key,
                                    marker,
                                    route,
                                    burial,
                                    split_pipe,
                                    'connector',
                                    'connector',
                                    marker+1,
                                    marker-3))

                                elec_sub_system.append((db_key, marker))

                                marker += 1
                                array_idx += 1
                                # add connectr to layout
                                location = (umbilical_data[
                                    'Device' + str(next_device+1).zfill(3)]
                                    ['termination'])

                                db_key, wet_mate_idx, dry_mate_idx = (
                                    self._add_connector(
                                        components, device_connection,
                                        wet_mate_idx, dry_mate_idx, db_key,
                                        marker, location, connection))

                                elec_sub_system.append((db_key, marker))
                                marker += 1

                                # add dynamic cable to layout
                                cable = umbilical_data[
                                    'Device' + str(next_device+1).zfill(3)]

                                self.umbilical_cables.append(UmbilicalCable(
                                    umbilical_idx,
                                    cable['length'],
                                    cable['db_key'],
                                    marker,
                                    cable['termination'],
                                    cable['device'],
                                    cable['x coords'],
                                    cable['z coords']))

                                elec_sub_system.append((cable['db_key'],
                                                        marker))

                                marker += 1

                                umbilical_idx += 1
                                marker += 1

                                # add device connector to layout
                                location = device_layout[
                                    'Device' + str(next_device+1).zfill(3)]

                                db_key, wet_mate_idx, dry_mate_idx = (
                                    self._add_connector(
                                        components, device_connection,
                                        wet_mate_idx, dry_mate_idx, db_key,
                                        marker, location, next_device))

                                elec_sub_system.append((db_key, marker))
                                marker += 1
                                hierarchy['device' +
                                          str(next_device+1).zfill(3)] = {}

                                (hierarchy['device' +
                                           str(next_device+1).zfill(3)]
                                 ['Elec sub-system']) = elec_sub_system

                                array_idx += 1
                                visited_nodes.append(chain_step)
                                device_to_device[start][chain_step] = 0
                                device_to_device[chain_step][start] = 0
                                start = chain_step

                            else:

                                chain = False

                        if self.n_cp > 1:

                            sub_hub_layout.append(layout)

                        else:

                            cluster['layout'].append(layout)

                    if self.n_cp > 1:

                        if self.substation:

                            if idx == 0:

                                pass

                            else:

                                # +1 offset for substation
                                try:

                                    hierarchy['subhub' +
                                              str(idx).zfill(3)]['layout'] = \
                                                  sub_hub_layout

                                except KeyError:

                                    msg = ("cant find subhub")
                                    module_logger.warning(msg)

                        else:

                            try:

                                hierarchy['subhub' +
                                          str(idx).zfill(3)]['layout'] = \
                                              sub_hub_layout

                            except KeyError:

                                msg = ("cant find subhub")
                                module_logger.warning(msg)

            else:

                for idx, cp in enumerate(self.cp_to_device):

                    if self.n_cp > 1:

                        sub_hub_layout = []

                    for connection in np.where(cp > 0)[0]:

                        layout = []
                        # find where this connects
                        visited_nodes.append(connection)
                        # then connect from here
                        chain = True
                        chain_step = connection
            #            device_to_device[:,connection] = 0

                        # Link to cp
                        # need to add reference to export side connector for
                        # installation - keep as ideal
                        array_connector = \
                            self.collection_points[idx].output_connectors

                        db_key, wet_mate_idx, dry_mate_idx = (
                            self._add_connector(
                                components, array_connector, wet_mate_idx,
                                dry_mate_idx,
                                self.collection_points[idx].db_key,
                                marker,
                                (self.collection_points[idx].utm_x,
                                 self.collection_points[idx].utm_y),
                                0))

                        link_to_cp = []
                        link_to_cp.append((db_key, marker))
                        marker += 1

                        if self.n_cp > 1:

                            cable_length = cp_device_distance[idx][chain_step]
                            route = cp_device_paths[idx][chain_step]

                        else:

                            cable_length = \
                                cp_device_distance[idx][chain_step+1]
                            route = cp_device_paths[idx][chain_step+1]

                        db_key = components['array']

                        try:

                            burial = (get_burial_depths(route,
                                                        burial_depths,
                                                        burial_array))

                            split_pipe = get_split_pipes(burial)

                        except TypeError:

                            errStr = (
                                "Components spacing too close to make cable "
                                "route. Consider changing device spacing "
                                "constraints, number of devices or selecting "
                                "a different network topology.")

                            raise RuntimeError(errStr)

                        self.array_cables.append(ArrayCable(array_idx,
                                                            cable_length,
                                                            db_key,
                                                            marker,
                                                            route,
                                                            burial,
                                                            split_pipe,
                                                            'device',
                                                            'collection point',
                                                            connection,
                                                            idx))

                        # add device to layout
                        link_to_cp.append((db_key, marker))
                        layout.append('device' + str(connection+1).zfill(3))
                        hierarchy['device' + str(connection+1).zfill(3)] = {}
                        (hierarchy['device' + str(connection+1).zfill(3)]
                            ['Elec sub-system']) = link_to_cp

                        marker += 1

                        # add device connector to layout
                        location = device_layout['Device' +
                                                 str(connection+1).zfill(3)]

                        db_key, wet_mate_idx, dry_mate_idx = (
                            self._add_connector(
                                components, device_connection, wet_mate_idx,
                                dry_mate_idx, db_key, marker, location,
                                connection))

                        (hierarchy['device' + str(connection+1).zfill(3)]
                            ['Elec sub-system']).append((db_key, marker))

                        marker += 1
                        array_idx += 1

                        start = connection
                        cp_to_device[idx][connection] = 0

                        while chain is True:

                            elec_sub_system = []

                            if np.any(device_to_device[chain_step] > 0):

                                next_device = \
                                    np.where(
                                        device_to_device[chain_step] > 0)[0]

                                # filter against visited nodes
                                for node in next_device:

                                    if node not in visited_nodes:

                                        next_device = node

                                layout.append('device' +
                                              str(next_device+1).zfill(3))

                                chain_step = int(next_device)

                                cable_length = \
                                    cp_device_distance[start+1][chain_step+1]
                                route = cp_device_paths[start+1][chain_step+1]

                                burial = (get_burial_depths(route,
                                                            burial_depths,
                                                            burial_array))

                                split_pipe = get_split_pipes(burial)

                                db_key = components['array']
                                self.array_cables.append(ArrayCable(
                                    array_idx,
                                    cable_length,
                                    db_key,
                                    marker,
                                    route,
                                    burial,
                                    split_pipe,
                                    'device',
                                    'device',
                                    chain_step,
                                    start))

                                elec_sub_system.append((db_key, marker))
                                marker += 1

                                # add device connector to layout
                                location = \
                                    device_layout['Device' +
                                                  str(next_device+1).zfill(3)]

                                db_key, wet_mate_idx, dry_mate_idx = (
                                    self._add_connector(
                                        components, device_connection,
                                        wet_mate_idx, dry_mate_idx, db_key,
                                        marker, location, next_device))

                                elec_sub_system.append((db_key, marker))
                                marker += 1
                                hierarchy['device' +
                                          str(next_device+1).zfill(3)] = {}

                                (hierarchy['device' +
                                           str(next_device+1).zfill(3)]
                                 ['Elec sub-system']) = elec_sub_system

                                array_idx += 1
                                visited_nodes.append(chain_step)
                                device_to_device[start][chain_step] = 0
                                device_to_device[chain_step][start] = 0
                                start = chain_step

                            else:

                                chain = False

                        if self.n_cp > 1:

                            sub_hub_layout.append(layout)

                        else:

                            cluster['layout'].append(layout)

                    if self.n_cp > 1:

                        if self.substation:

                            if idx == 0:

                                pass

                            else:

                                # +1 offset for substation
                                try:

                                    hierarchy['subhub' +
                                              str(idx).zfill(3)]['layout'] = \
                                                  sub_hub_layout

                                except KeyError:

                                    msg = ("cant find subhub {}".format(idx))
                                    module_logger.warning(msg)

                        else:

                            try:

                                hierarchy['subhub' +
                                          str(idx).zfill(3)]['layout'] = \
                                              sub_hub_layout

                            except KeyError:

                                msg = ("cant find subhub {}".format(idx))
                                module_logger.warning(msg)

        array.append(cluster)

        hierarchy['array'] = array
        self.all_connections = hierarchy

        return

    def make_hierarchy(self, n_device):

        '''Make the network hierarchy for downstream analysis.

        Args:
            n_device (int): The number of OECs in the array.

        Atttributes:
            hier (dict): Network connection hierarchy for downstream analysis.
            array (list): Container for array level information, to be stored
                in heirarchy.
            sub_array (dict): Container for sub array level information, to be
                stored in array.

        Returns:
            none.

        '''

        hier = {}
        index = 'db'

        # devices
        for oec in range(n_device):

            hier['device' + str(oec+1).zfill(3)] = {}
            hier['device' + str(oec+1).zfill(3)]['Elec sub-system'] = []

            local_system = []

            for item in self.all_connections[
                    'device' + str(oec+1).zfill(3)]['Elec sub-system']:

                local_system.append(item[0])

            hier['device' + str(oec+1).zfill(3)]['Elec sub-system'].append(
                local_system)

        # array
        for item in self.all_connections['array']:

            sub_array = dict.fromkeys(['Substation', 'Export cable', 'layout'])

            sub_array['layout'] = item['layout']

            for system in ['Substation', 'Export cable']:

                sub_array[system] = \
                    [self.get_network_components(system, item, index)]

        hier['array'] = sub_array

        # sub hubs
        for key, val in self.all_connections.iteritems():

            if 'subhub' in key:

                sub_array = \
                    dict.fromkeys(['Substation', 'Elec sub-system', 'layout'])

                sub_array['layout'] = val['layout']

                for system in ['Substation', 'Elec sub-system']:

                    sub_array[system] = \
                        [self.get_network_components(system, val, index)]

                hier[key] = sub_array

        self.hierarchy = hier

        return

    def make_network_design(self, n_device):

        '''Make the network design table for downstream analysis.

        '''

        design = {}

        for oec in range(n_device):

            design['device' + str(oec+1).zfill(3)] = {}
            design['device' + str(oec+1).zfill(3)] = {'marker': []}

            local_system = []

            for item in (self.all_connections[
                    'device' + str(oec+1).zfill(3)]['Elec sub-system']):

                local_system.append(item[1])

            design['device' + str(oec+1).zfill(3)]['marker'].append(
                    local_system)

            design['device' + str(oec+1).zfill(3)].update(
                {'quantity': Counter(component[0] for
                                     component in
                                     self.all_connections['device' +
                                     str(oec+1).zfill(3)]['Elec sub-system'])})

        # sub hubs
        for key, val in self.all_connections.iteritems():

            if 'subhub' in key:

                design[key] = {}
                design[key] = {'marker': []}

                local_markers = []
                local_keys = []

                for item in self.all_connections[key]['Elec sub-system']:

                    local_markers.append(item[1])
                    local_keys.append(item[0])

                for item in self.all_connections[key]['Substation']:

                    local_markers.append(item[1])
                    local_keys.append(item[0])

                design[key]['marker'].append(local_markers)

                counter = Counter()

                for item in local_keys:

                    counter[item] += 1

                design[key].update({'quantity': counter})

        # array
        index = 'marker'

        for item in self.all_connections['array']:

            sub_array = dict.fromkeys(['Substation', 'Export cable'], {})

            for system in ['Substation', 'Export cable']:

                components = \
                    self.get_network_components(system, item, index)

                sub_array[system] = {'marker': [components]}

                count = Counter(component[0] for component in item[system])

                sub_array[system].update({'quantity': count})

        design['array'] = sub_array

        self.network_design = design

        return

    def get_network_components(self, system_name, system, index):

        if index == 'db':

            i = 0

        else:

            i = 1

        local_system = [component[i] for component in system[system_name]]

        return local_system

    def make_bom(self):

        '''Make the network bill of materials for downstream analysis. For each
        component get the marker, db ref, type, utm x, utm y and quantity. Then
        collate in pandas table.

        Attributes:
            markers (list): List of all component markers.
            db_key (list): List of all component database keys.
            install_type (list): List of all component types.
            utm_x (list): List of all component locations, x coordinate.
            utm_y (list): List of all component locations, y coordinate.
            quantity (list): Quantity of all components. This is '1 unit' for
                all components except cables which gives the cable length.
            b_o_m_dict (dict): Structured dictionary for converting to pandas
                dataframe.

        Note:
            This could be improved by making the list creation a function?

        '''

        markers = []
        db_key = []
        install_type = []
        utm_x = []
        utm_y = []
        quantity = []

        # for all possible components
        for item in self.wet_mate:
            markers.append(item.marker)
            db_key.append(item.db_key)
            install_type.append(item.type_)
            utm_x.append(item.utm_x)
            utm_y.append(item.utm_y)
            quantity.append(1)

        for item in self.dry_mate:
            markers.append(item.marker)
            db_key.append(item.db_key)
            install_type.append(item.type_)
            utm_x.append(item.utm_x)
            utm_y.append(item.utm_y)
            quantity.append(1)

        for item in self.export_cables:
            markers.append(item.marker)
            db_key.append(item.db_key)
            install_type.append(item.type_)
            utm_x.append('None')
            utm_y.append('None')
            quantity.append(item.length)

        for item in self.array_cables:
            markers.append(item.marker)
            db_key.append(item.db_key)
            install_type.append(item.type_)
            utm_x.append('None')
            utm_y.append('None')
            quantity.append(item.length)

        for item in self.collection_points:
            markers.append(item.marker)
            db_key.append(item.db_key)
            install_type.append(item.type_)
            utm_x.append(item.utm_x)
            utm_y.append(item.utm_y)
            quantity.append(1)

        for item in self.umbilical_cables:
            markers.append(item.marker)
            db_key.append(item.db_key)
            install_type.append(item.type_)
            utm_x.append('None')
            utm_y.append('None')
            quantity.append(item.length)

        # make the pandas table
        b_o_m_dict = {"marker": markers,
                      "db ref": db_key,
                      "install_type": install_type,
                      "utm_x": utm_x,
                      "utm_y": utm_y,
                      "quantity": quantity}

        self.b_o_m = pd.DataFrame(b_o_m_dict)

        return

    def _get_db_keys_from_pd(self):

        '''Get all database keys of all components used in the array.

        Attributes:
            keys (set) [-]: DB keys.

        Returns:
            list [-]: List keys.

        '''

        keys = set(self.b_o_m['db ref'])

        return list(keys)

    def set_economics_data(self, db, onshore_cost):

        '''Compile network design data into economics bill of materials.

        Args:
            db (pd) [-]: The component database.

        Attributes:
            network_keys (list) [-]: Database keys of all components used in
                the array.
            quantity (list) [-]: Total quantity of each unique component
                database key used in the array.
            economics_dict (dict) [-]: Structured dictionary for converting to
                pandas dataframe.

        Returns:
            none.

        '''

        network_keys = self._get_db_keys_from_pd()

        quantity = []
        type_ = []

        for key in network_keys:

            type_.append(
                self.b_o_m[
                    self.b_o_m['db ref'] == key].install_type.values.tolist())

            quantity.append(
                self.b_o_m[self.b_o_m['db ref'] == key].sum()["quantity"])

        type_ = [item[0] for item in type_]  # get item type without using set
        type_ = self._map_component_types(type_)

        cost = self._get_costs_from_db(db, network_keys, type_, quantity)

        if onshore_cost:

            network_keys.append(None)
            type_.append('Onshore')
            quantity.append(1)
            cost.append(onshore_cost)

        economics_dict = {"db ref": network_keys,
                          "quantity": quantity,
                          "cost": cost,
                          "year": [0]*len(network_keys)}

        self.economics_data = pd.DataFrame(economics_dict)

        return

    def _map_component_types(self, type_):

        '''Map component types. Required to ensure compatibility between names
        used in the install modules and the electrical module.

        Args:
            type_ (list) [-]: Component type list.

        Attributes:
            name_map (dict) [-]: Key = installation module label,
                                 value = electrical database label.

        Returns:
            list [-]: Component types to match electrical database.

        '''

        name_map = {'export': 'static_cable',
                    'array': 'static_cable',
                    'wet-mate': 'wet_mate_connectors',
                    'dry-mate': 'dry_mate_connectors',
                    'substation': 'collection_points',
                    'passive hub': 'collection_points',
                    'umbilical': 'dynamic_cable',
                    }

        return ([name_map[item] for item in type_])

    def _get_costs_from_db(self, db, keys, type_, quantity):

        '''For each component, extract unit cost from the database.

        Args:
            db () [-]:
            keys () [-]:
            type_ () [-]:
            quantity () [-]:

        Attributes:
            all_cost (list) [E]: Total component cost.

        Returns:
            all_cost

        Note: The commented-out code will calculate the total cost of each
              component based on quantity and unit cost.

        '''

        all_cost = []

#        for component in zip(keys, type_, quantity):
#            db_dict = getattr(db, component[1])
#            cost = db_dict[db_dict.id == component[0]].cost.values[0]
#            all_cost.append(cost * component[2])

        for component in zip(keys, type_):
            db_dict = getattr(db, component[1])
            all_cost.append(db_dict[db_dict.id == component[0]].cost.values[0])
#            all_cost.append(cost * component[2])
        return all_cost

    def total_network_cost(self):

        '''Calculate total network cost and set total cost attribute.

        Returns:
            none.

        '''

#        self.total_cost = self.economics_data['cost'].sum(axis=0)
        self.total_cost = sum(
            self.economics_data.cost * self.economics_data.quantity)

        return

    def make_cable_routes(self, grid, all_x, all_y):

        '''Collect the cable routes in pd.DataFrame for downstream analysis.

        Args:
            grid
            all_x
            all_y

        '''

        marker = []
        db_ref = []
        grid_id = []
        burial_depth = []
        split_pipe = []

        for cable in (self.array_cables + self.export_cables):

            marker += [cable.marker] * len(cable.route)
            db_ref += [cable.db_key] * len(cable.route)
            burial_depth += cable.target_burial_depth
            split_pipe += cable.split_pipe
            grid_id += cable.route

        cable_x, cable_y = \
            self._convert_path_to_coordinates(grid_id, all_x, all_y)
            
        indexed_grid = grid.set_index("id")
        depth_cols = ['layer 1 start'] * len(grid_id)
        type_cols = ['layer 1 type'] * len(grid_id)
        
        cable_depth = indexed_grid.lookup(grid_id, depth_cols)
        cable_type = indexed_grid.lookup(grid_id, type_cols)

        cable_dict = {"marker": marker,
                      "db ref": db_ref,
                      "burial_depth":  burial_depth,
                      "split pipe": split_pipe,
                      "x": cable_x,
                      "y": cable_y,
                      'layer 1 start': cable_depth,
                      'layer 1 type': cable_type
                      }
                
        self.cable_routes = pd.DataFrame(cable_dict)

        return

    def make_collection_point_design(self):

        '''Collection point output data structure.

        Args:
            none

        Attributes:
            centre_of_gravity (np.array) [m]: Centre of gravity with respect to
                local coordinate system, as [x,y,z].
            dry beam area (float) [m2]: Dry beam area.
            dry frontal area (float) [m2]: Dry frontal area.
            foundation type (str) [-]: Predefined foundation type, either
                'gravity' or 'pile'.
            foundation locations (np.array) [m]: Foundation location with
                respect to the origin point, as [x,y,z].
            height (float) [m]: Unit height.
            length (float) [m]: Unit length.
            marker (int) [-]: Collection point unique marker.
            mass (float) [kg]: Total unit weight.
            orientation angle (float) [deg]: Unit orientation angle.
            origin (np.array) [m]: Collection point origin point, UTM
                coordinates as [x,y].
            profile (str) [-]: Shape, either 'cylindrical' or 'rectangular'.
            surface roughness (float) [m]: Surface roughness.
            type_ (str) [-]: Collection point type, either 'subsea passive',
                'subsea substation' or 'surface substation'.
            volume (float) [m3]: Submerged volumne.
            wet beam area (float) [m2]: Wet beam area.
            wet frontal area (float) [m2]: Wet frontal area.
            width (float) [m]: Unit width.

        Returns:
            none

        Note:
            All attributes within container list for conversion to pandas
            DataFrame.

        '''

        marker = []
        origin = []
        operating_environment = []
        foundation_type = []
        type_ = []
        mass = []
        volume = []
        centre_of_gravity = []
        wet_frontal_area = []
        wet_beam_area = []
        dry_frontal_area = []
        dry_beam_area = []
        length = []
        width = []
        height = []
        profile = []
        surface_roughness = []
        orientation_angle = []
        foundation_locations = []

        for cp in self.collection_points:

            marker.append(cp.marker)
            origin.append(np.array((cp.utm_x, cp.utm_y)))
            foundation_type.append(cp.foundation_type)
            operating_environment.append(cp.operating_environment)
            mass.append(cp.mass)
            length.append(cp.length)
            width.append(cp.width)
            height.append(cp.height)
            volume.append(cp.volume)
            centre_of_gravity.append(np.array((0., 0., 0.)))
            profile.append(cp.profile)
            wet_frontal_area.append(cp.wet_frontal_area)
            wet_beam_area.append(cp.wet_beam_area)
            dry_frontal_area.append(cp.dry_frontal_area)
            dry_beam_area.append(cp.dry_beam_area)
            surface_roughness.append(cp.surface_roughness)
            orientation_angle.append(cp.orientation_angle)
            foundation_locations.append(np.array((0., 0., 0.)))

        collection_point_dict = {"marker": marker,
                                 "origin": origin,
                                 "type": operating_environment,
                                 "mass": mass,
                                 "volume": volume,
                                 "length": length,
                                 "width": width,
                                 "height": height,
                                 "centre_of_gravity": centre_of_gravity,
                                 "profile": profile,
                                 "wet frontal area": wet_frontal_area,
                                 "wet beam area": wet_beam_area,
                                 "dry frontal area": dry_frontal_area,
                                 "dry beam area": dry_beam_area,
                                 "surface roughness": surface_roughness,
                                 "orientation angle": orientation_angle,
                                 "foundation locations": foundation_locations}

        self.collection_points_design = pd.DataFrame(collection_point_dict)

        return

    def calculate_power_quantities(self, ideal_yield, ideal_histogram):

        '''Processing of power quantities for network assessment.

        Args:
            ideal_yield (float): Array power output for a year period assuming
                no electrical losses.

            ideal_histogram (list): Array power output at each power generation
                level assuming no electrical losses.

        Returns:
            none.

        '''

        self.annual_yield = self.calculate_annual_yield()
        self.annual_losses = self.calculate_annual_losses(ideal_yield)
        self.annual_efficiency = self.annual_yield / ideal_yield
        self.histogram_losses = self.calculate_histogram_losses(
                                                            ideal_histogram)
        self.histogram_efficiency = self.calculate_histogram_efficiency(
                                                            ideal_histogram)

        return

    def calculate_annual_yield(self):

        '''Calculate annual energy yield.

        Returns:
            annual_yield (float):  Array power output for a year period with
                electrical losses.

        '''

        year_hours = 365 * 24
        annual_yield = 0

        for time, power in zip(self.power_histogram,
                               self.array_power_output):
            
            # Convert from MW to W and ignore directionality
            power_w = abs(power * 1e6)

            annual_yield += time * year_hours * power_w
        
        # Correct for bad calculations
        if np.isnan(annual_yield): annual_yield = 0.

        return annual_yield

    def calculate_annual_losses(self, ideal_yield):

        '''Calculate annual energy losses by comparing against ideal.

        Args:
            ideal_yield (float): Array power output for a year period assuming
                no electrical losses.

        Returns:
            annual_losses (float): Electrical losses for a year period.

        '''

        annual_losses = ideal_yield - self.annual_yield

        return annual_losses

    def calculate_histogram_losses(self, ideal_histogram):

        '''Calculate histogram losses by comparing against ideal.

        Args:
            ideal_histogram (list): Array power output at each power generation
                level assuming no electrical losses.

        Returns:
            histogram_losses (list): Electrical losses at each power
                generation level.

        '''

        histogram_losses = [ideal-(abs(actual)*1000000) for
                            ideal, actual in
                            zip(ideal_histogram, self.array_power_output)]

        return histogram_losses

    def calculate_histogram_efficiency(self, ideal_histogram):

        '''Calculate array efficiency at each bin in the power histogram.

        Args:
            ideal_histogram (list): Array power output at each power generation
                level assuming no electrical losses.

        Returns:
            histogram_efficiency (list): Electrical efficiency at each power
                generation level.

        '''

        histogram_efficiency = [(abs(actual)*1000000)/ideal for
                                actual, ideal in
                                zip(self.array_power_output, ideal_histogram)]

        return histogram_efficiency
    
    def calculate_lcoe(self):

        '''Simplified LCOE for comparison of electrical networks.
        '''
        
        if self.annual_yield == 0.:
            self.lcoe = np.inf
        else:
            self.lcoe = self.total_cost / self.annual_yield * 1e3
            
        module_logger.debug("LCOE: {}".format(self.lcoe))

        return

    def make_umbilical_table(self):

        '''Quick fix to make the umbilical data table.

        For each umbilical get: the marker, db ref, device, seabed connection
        point and length. Then Collate in a pandas table.

        Attributes:
            marker (list): Umbilical cable markers.
            db_key (list): Umbilical cable database keys.
            device (list): Associated device.
            seabed_connection_point (list): Umbilical cable connection points.
            length (list): Umbilical cable lengths.
            umbilical_dict (dict): Structured dictionary for converting to
                pandas dataframe.

        Note:
            This could be improved by making the list creation a function?

        '''

        marker = []
        db_key = []
        device = []
        seabed_connection_point = []
        length = []

        for cable in self.umbilical_cables:

            marker.append(cable.marker)
            db_key.append(cable.db_key)
            device.append(cable.device)
            seabed_connection_point.append(
                np.array([cable.seabed_termination_x,
                          cable.seabed_termination_y,
                          cable.seabed_termination_z]))
            length.append(cable.length)

        umbilical_dict = {"marker": marker,
                          "db ref": db_key,
                          "device": device,
                          "seabed_connection_point": seabed_connection_point,
                          "length": length}

        self.umbilical_cable_design = pd.DataFrame(umbilical_dict)

        return

    def _make_result_str(self):

        msg = "\n"
        msg += "Annual yield: {}\n\n".format(self.annual_yield)
        msg += "Bill of Materials:\n\n"
        msg += "{}\n\n".format(self.economics_data)
        msg += "Component Data:\n\n"
        msg += "{}\n\n".format(self.b_o_m)
        msg += "Hierarchy:\n\n"
        msg += "{}\n\n".format(pformat(self.hierarchy))
        msg += "Network design:\n\n"
        msg += "{}\n\n".format(pformat(self.network_design))
        msg += "Cable routes:\n\n"
        msg += "{}\n\n".format(self.cable_routes)
        msg += "Collection points:\n\n"
        msg += "{}\n\n".format(self.collection_points_design)
        msg += "Umbilical cables:\n\n"
        msg += "{}\n\n".format(self.umbilical_cable_design)

        return msg

    def print_result(self):

        msg = self._make_result_str()
        print msg

        return

    def log_result(self):

        msg = self._make_result_str()
        module_logger.info(msg)

        return

    def _convert_path_to_coordinates(self, grid_id, all_x, all_y):

        '''Convert a list of grid points into x and y coordinates.

        Args:
            grid_id (list) [-]: List of grid point ids.
            all_x (list) [m]: List of all x coordinates in the area.
            all_y (list) [m]: List of all y coordinates in the area.

        Attributes:
            x (list) [m]: List of x coordinates traversed by cables.
            y (list) [m]: List of y coordinates traversed by cables.

        Returns:
            x
            y

        Note:
            Faster to perform two list separate list comprehensions?
            Grid point id is set at input of module for internal use only.

        '''

        (x, y) = zip(*[(all_x[point], all_y[point]) for point in grid_id])

        return x, y

    def _add_connector(self, components, device_connection, wet_mate_idx,
                       dry_mate_idx, db_key, marker, location, start):

        '''Define connector in the network.

        Args:

        Attributes:

        Returns:

        '''

        db_key = components['connector']

        if device_connection == 'wet-mate':
            self.wet_mate.append(WetMateConnector(wet_mate_idx,
                                                  db_key,
                                                  marker,
                                                  location))

            wet_mate_idx += 1

        elif device_connection == 'dry-mate':
            self.dry_mate.append(DryMateConnector(dry_mate_idx,
                                                  db_key,
                                                  marker,
                                                  location))

            dry_mate_idx += 1

        return db_key, wet_mate_idx, dry_mate_idx
