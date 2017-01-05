# -*- coding: utf-8 -*-
"""
This module converts the electrical network into a format compatible with
PyPower.

.. module:: power_flow_v2
   :platform: Windows
   :synopsis: Convert network into PyPower object.

.. moduleauthor:: Adam Collin <a.collin@ed.ac.uk>

"""
import numpy as np
import copy
from pypower import runpf, ppoption
import itertools

class PyPower(object):
    
    '''Container class for all information required to construct and simulate 
    an electrical network using PyPower.
    
    Args:
        n_cp (int) [-]: The number of collection points in the network.
        n_devices (int) [-]: The number of devices in the network.
        network_connections (dict) [-]: Network connections, keys:
            shore_to_device, shore_to_cp, cp_to_cp, cp_to_device,
            device_to_device.
        export_voltage (float) [V]: Export system voltage.
        array_voltage (float) [V]: Array system voltage 
        device_voltage (float) [V]: Device system voltage.
        substation (bool) [-]: True = Substation present, False = Substation
            not present.
        topology (str) [-]: Type of network to be analysed, either 'Radial' or
            'Star'.
        floating (bool) [-]: True = Floating devices, False= Fixed devices.

    Attributes:
        shore_to_device (np.array) [-]: Matrix of shore to device connections.
            Dimensions: [n_devices x 1].
        shore_to_cp (np.array) [-]: Matrix of shore to collection point
            connections. Dimensions: [n_cp x 1].
        cp_to_cp (np.array) [-]: Matrix of collection point to collection point
            connections. Dimensions: [n_cp x n_cp].
        cp_to_device (np.array) [-]: Matrix of collection point to device
            connections. Dimensions: [n_cp x n_devices].
        device_to_device (np.array) [-]: Matrix of device to device
            connections. Dimensions: [n_devices x n_devices].
        n_branch (int) [-]: The number of branches in the PyPower network.
        n_bus (int) [-]: The number of buses in the PyPower network.
        gen_data (np.array) [-]: PyPower generator data structure.
        bus_data (np.array) [-]: PyPower bus data structure.
        branch_data (np.array) [-]: PyPower branch data structure.
        flag (list) [-]: Boolean operator to indicate power flow succees;
            1 = True, 0 = False.
        onshore_active_power
        onshore_reactive_power
        busbar_voltages
        busbar_angles
        gen_active_power 
        gen_reactive_power 

    Returns:

    Note:
        This has been developed for PyPower v5.0.1. Further information on the
        PyPower data structures is available at:
        http://www.pserc.cornell.edu/matpower/MATPOWER-manual.pdf

    '''

    def __init__(self,
                 n_cp,
                 n_devices,
                 network_connections,
                 export_voltage, 
                 array_voltage,
                 device_voltage,
                 substation,
                 topology,
                 floating):

        self.n_cp = n_cp
        self.n_devices = n_devices
        self.shore_to_device = network_connections['shore_to_device']
        self.device_to_device = network_connections['device_to_device']
        self.shore_to_cp = network_connections['shore_to_cp']
        self.cp_to_cp = network_connections['cp_to_cp']
        self.cp_to_device = network_connections['cp_to_device']
        self.export_voltage = export_voltage
        self.array_voltage = array_voltage
        self.device_voltage = device_voltage
        self.substation = substation
        self.topology = topology
        self.floating = floating
        self.n_branch = None
        self.n_bus = None
        self.gen_data = None
        self.bus_data = None
        self.branch_data = None
        self.flag = None
        self.onshore_active_power = None
        self.onshore_reactive_power = None
        self.busbar_voltages = None
        self.busbar_angles = None
        self.gen_active_power = None
        self.gen_reactive_power = None
        self.MVAb = 100e6 # fixed system base
        self.all_results = None
        self.export_branches = []
        self.array_branches = []
        self.device_branches = []

    def build_network(self, z_export, z_array, z_device, T_export_array, 
                      T_array_device, z_umbilical):
        
        '''Build the network for use in PyPower here by calling connection
        functions. Impedance values are passed to this function as a network
        configuration can be built using different cables.

        Args:
            z_export (list) [pu]: Impedance of export cable, val1 = resistance,
                val2 = reactance, val3 = susceptance.
            z_array (list) [pu]: Impedance of array cable, val1 = resistance,
                val2 = reactance, val3 = susceptance.
            z_device (list) [pu]: Impedance of device cable, val1 = resistance,
                val2 = reactance, val3 = susceptance.
            T_export_array(float) [pu]: Reactance of export to array
                transformer.
            T_array_device (float) [pu]: Reactance of array to device
                transformer.

        Returns:

        Note:
            This calls the following methods:
                branch_count,
                set_generators,
                label_buses,
                connect_branches

        '''

        self.n_branch = self.branch_count()
        self.n_bus = self.n_branch+1 # offset is for the slack bus
        self.gen_data = self.set_generators()
        self.bus_data = self.label_buses()
        self.branch_data = self.connect_branches(z_export, z_array, z_device,
                                                 T_export_array,
                                                 T_array_device, z_umbilical)

        return

    def connect_branches(self, z_export, z_array, z_device, T_export_array,
                         T_array_device, z_umbilical):
        
        '''Logic control to guide network branch connection process.
        
        Args:
            z_export (list) [pu]: Impedance of export cable, val1 = resistance,
                val2 = reactance, val3 = susceptance.
            z_array (list) [pu]: Impedance of array cable, val1 = resistance,
                val2 = reactance, val3 = susceptance.
            z_device (list) [pu]: Impedance of device cable, val1 = resistance,
                val2 = reactance, val3 = susceptance.
            T_export_array(float) [pu]: Reactance of export to array
                transformer.
            T_array_device (float) [pu]: Reactance of array to device
                transformer.
        
        Attributes:
            branch_data (np.array) [-]: PyPower branch data structure.

        Returns:
            branch_data

        Note:
            This calls the following methods:
                shore_to_device_to_device,
                shore_to_cp_fn,
                cp_to_cp_fn,
                cp_to_device_to_device            

        '''

        branch_data = np.array([[0]*13]*(self.n_branch), dtype=float)

        if self.n_cp == 0:

            branch_data, static_ends, device_order, n = \
                self.shore_to_device_to_device(branch_data, z_export, z_device)

        else:

            branch_data, cp_ends, n = self.shore_to_cp_fn(branch_data,
                                                          z_export,
                                                          T_export_array)

            if self.n_cp > 1:

                branch_data, cp_ends, n = self.cp_to_cp_fn(branch_data, 
                                                           cp_ends, n, z_array,
                                                           T_array_device)

            branch_data, static_ends, device_order, n  = \
                self.cp_to_device_to_device(branch_data, cp_ends, n, z_device)

        if self.floating:
            
            branch_data = \
                self.add_umbilicals(branch_data, static_ends, device_order, n,
                                    z_umbilical)

        return branch_data

    def branch_count(self):

        '''Set the number of branches in the Pypower network. This is
        determined by the number of components and the number of voltage levels
        in the network.
    
        Args:

        Attributes:
            n_voltage_levels (int) [-]: The number of voltage levels within the
                network.
            n_cp (int) [-]: The number of collection points in the network.
            n_devices (int) [-]: The number of devices in the network.
            shore_to_cp (np.array) [-]: Connection matrix of shore to cp
                system.
            shore_links (int) [-]: Number of cables between shore and cp.
            n_branch (int) [-]: The number of branches in the network.

        Returns:
            n_branch

        Note:
            cp = collection points.

        '''

        n_voltage_levels = len(set([self.export_voltage,
                                    self.array_voltage,
                                    self.device_voltage]))

        if self.n_cp == 0:
            # if no cp, all devices are single serving connections.
            n_branch = self.n_devices

        elif self.n_cp == 1:
            # if one cp, check number of voltage levels
            n_branch = self.n_devices + n_voltage_levels

        else:
            # if more than one cp, check number of voltage levels
            if n_voltage_levels == 1:
                # export == array == device
                n_branch = self.n_cp + self.n_devices

            elif n_voltage_levels == 3:
                # export != array != device
                n_branch = self.n_cp*2 + self.n_devices

            else:

                shore_links = np.sum(self.shore_to_cp)

                if self.export_voltage == self.array_voltage:
                    # export == array != device
                    n_branch = ((self.n_cp - shore_links)*2 + self.n_devices +
                        shore_links)

                else:
                    # export != array == device
                    n_branch = ((self.n_cp - shore_links) + self.n_devices + 
                        shore_links*2)

        if self.floating == True:
            
            n_branch = n_branch + self.n_devices

        return n_branch

    def label_buses(self):
        
        '''Define busbars by setting the base voltage. Three voltage areas are
        defined: export, array and device.

        Args:

        Attributes:
            last_export_bus (int) [-]: Bus identifier for the last bus in the
                export area.
            first_device_bus (int) [-]: Bus identifier for the first bus in the
                device area.
            bus_data (np.array) [-]: PyPower bus data structure.

        Returns:
            bus_data

        '''

        last_export_bus = 1+max(int(np.sum(self.shore_to_cp)),
                                int(np.sum(self.shore_to_device)))

        first_device_bus = (self.n_bus - self.n_devices)

        bus_data = np.array([[0]*13]*self.n_bus, dtype=float)

        # slack bus
        bus_data[0][[1,9]] = [3, self.export_voltage]

        # set export bus voltages
        bus_data[:,9][1:last_export_bus] = self.export_voltage

        # set array bus voltages - may need to look at this for array settings
        bus_data[:,9][1+int(np.sum(self.shore_to_cp)):first_device_bus] = (
            self.array_voltage)

        # set device bus voltages
        bus_data[:,9][first_device_bus:self.n_bus] = self.device_voltage

        # bus constant non zero values
        bus_data[:,0] = range(1,self.n_bus + 1) 
        bus_data[:,1][1:self.n_bus] = 1 # all buses are PQ buses
        bus_data[:,[6, 7, 10, 11, 12]] = 1, 1.0, 1, 1.1, 0.9

        return bus_data

    def set_generators(self):
        
        '''Add generators to buses. Generators are connected to the last N
        buses, where N is the number of generators. A slack generator is always
        connected to bus 0.

        Args:

        Attributes:
            gen_data (np.array) [-]: PyPower generator data structure.
            connect_bus (int) [-]: Bus identifier for the oec device.

        Returns:
            gen_data

        '''

        # initialise empty data structure
        gen_data = np.array([[0]*21]*(self.n_devices+1), dtype=float)
        # set generator at slack bus        
        gen_data[0][:2] = 1
        # set oec bus
        connect_bus = self.n_bus - self.n_devices + 1 # offset for slack bus
        gen_data[:,0][1:self.n_devices+1] = range(connect_bus, self.n_bus + 1)
        # generator constant values
        gen_data[:,5:8] = [1, 100, 1]

        return gen_data

    def shore_to_device_to_device(self, branch_data, z_export, z_device):

        '''Connect from shore to device (to device). Update branch data with
        these connections.

        Args:
            branch_data (np.array) [-]: PyPower branch data structure.
            z_export (tuple) [pu]: Export cable impedance.
            z_device (np.array) [pu]: Impedances expressed in matrix form.

        Attributes:
            device_to_device (np.array) [-]: Matrix of device to device
                connections.
            branch_count (int) [-]: Global counter for branches.
            branch_start (int) [-]: Branch start busbar number.
            branch_end (int) [-]: Branch end busbar number.
            visited_nodes (list) [-]: List of nodes visited for ignore.
            chain_step (int) [-]: Indicates device to be considered.

        Returns:
            branch_data

        '''

        device_to_device = copy.deepcopy(self.device_to_device)
        branch_count = 0 
        branch_start = 1
        branch_end = 2
        visited_nodes = []
        static_ends = []

        for connection in np.where(self.shore_to_device > 0)[0]:
            # find where this connects
            visited_nodes.append(connection)
            # then connect from here
            chain_step = connection
            device_to_device[:,connection] = 0
            # create branch
            branch_data[branch_count][:5] = [branch_start, branch_end,
                                             z_export[0], z_export[1],
                                             z_export[2]]

            static_ends.append(branch_end)
            branch_count += 1

            (branch_end, branch_start, visited_nodes, branch_count,
             device_to_device, branch_data, static_ends) = \
                 self.device_to_device_fn(
                     device_to_device, chain_step, visited_nodes, branch_start,
                     branch_end, branch_data, branch_count, True, z_device,
                     static_ends)

        # set branch constant values
        branch_data[:,[5,6,7,10,11,12]] = [0, 0, 0, 1, -360, -360]

        return branch_data,static_ends, visited_nodes, branch_count

    def device_to_device_fn(self, device_to_device, chain_step, visited_nodes,
                            branch_start, branch_end, branch_data,
                            branch_count, shore, z_matrix, static_ends):

        ''' Connect device to device. Update branch data with these
        connections.

        Args:
            device_to_device (np.array) [-]: Matrix of device to device
                connections.
            chain_step (int) [-]: Indicates device to be considered.
            visited_nodes (list) [-]: List of nodes visited for ignore.
            branch_start (int) [-]: Branch start busbar number.
            branch_end (int) [-]: Branch end busbar number.
            branch_data (np.array) [-]: PyPower branch data structure.
            branch_count (int) [-]: Global counter for branches.
            shore (bool) [-]: Boolean operator indicating presence of shore
                connection.
            z_matrix (np.array) [pu]: Impedances expressed in matrix form.

        Attributes:
            chain (bool) [-]: Boolean operator indicating presence of chain of
                devices.
            next_devices (list) [-]: List of connected devices.
            next_device (int) [-]: Reduced version of next_devices.

        Returns:
            branch_end
            branch_start
            visited_nodes
            branch_count
            device_to_device
            branch_data

        '''

        chain = True

        while chain == True:

            if np.any(device_to_device[chain_step] > 0):

                next_devices = np.where(device_to_device[chain_step] > 0)[0]
                #filter against visited nodes
                for node in next_devices:

                    if node not in visited_nodes:

                        next_device = node

                branch_start = branch_end
                branch_end += 1

                # create branch
                z_device = z_matrix[chain_step+1][next_device+1]# +1 substation
                branch_data[branch_count][:5] = [branch_start, branch_end,
                                                 z_device[0], z_device[1],
                                                 z_device[2]]

                self.array_branches.append(branch_end)

                static_ends.append(branch_end)
                branch_count += 1
                chain_step = int(next_device)
                visited_nodes.append(chain_step)
                device_to_device[:,chain_step] = 0

            else:

                chain = False

                if shore == True:

                    branch_start = 1

                branch_end += 1

        return (branch_end, branch_start, visited_nodes, branch_count,
                device_to_device, branch_data, static_ends)

    def shore_to_cp_fn(self, branch_data, z_export, T_export_array):
        
        '''Connect from shore to collection point(s). Update branch data with
        these connections.

        Args:
            branch_data (np.array) [-]: PyPower branch data structure.
            z_export (list) [pu]: Impedance of export cable.
            T_export_array (float) [pu]: Impedance of export to array
                transformer.

        Attributes:
            n_shore_links (int) [-]: Number of cables between shore and cp.
            branch_start (int) [-]: Branch start busbar number.
            branch_end (int) [-]: Branch end busbar number.
            branch_count (int) [-]: Global counter for branches.
            cp_ends (list) [-]: List of busbar numbers at 'end' of cp, i.e. at
                point of connection to cable.

        Returns:
            branch_data
            cp_ends
            branch_count

        '''

        n_shore_links = np.sum(self.shore_to_cp)
        branch_start = 1
        branch_end = 2
        branch_count = 0
    
        cp_ends = [0]*np.size(self.shore_to_cp)
    
        for connection in np.where(self.shore_to_cp > 0)[0]:
            
            branch_data[branch_count][:5] = \
                [branch_start, branch_end, z_export[0],z_export[1],z_export[2]]

            cp_ends[connection] = branch_end
            
            self.export_branches.append(branch_end)
            
            branch_end += 1
            branch_count += 1
            
            
    
        branch_start = 2
        
        if not self.export_voltage == self.array_voltage:

            for cp in range(n_shore_links):
                # add transformers to these collection points
                # 0 is for resistance, i.e. ideal transformer
                branch_data[branch_count][:4] = \
                    [branch_start, branch_end, 0, T_export_array]

                branch_start += 1
                branch_end += 1
                branch_count += 1

            cp_ends = [x+n_shore_links if x > 0 else 0 for x in cp_ends]

        return branch_data, cp_ends, branch_count
    
    def cp_to_cp_fn(self, branch_data, cp_ends, branch_count, z_matrix,
                    T_array_device):

        '''Connect from collection point to collection point. Update branch
        data with these connections.

        Args:
            branch_data (np.array) [-]: PyPower branch data structure.
            cp_ends (list) [-]: List of busbar numbers at 'end' of cp, i.e. at
                point of connection to cable.
            branch_count (int) [-]: Global counter for branches.
            z_matrix (np.array) [pu]: Impedances expressed in matrix form.
            T_array_device (float) [pu]: Impedance of array to device
                transformer.

        Attributes:
            cp_to_cp (np.array) [-]: Matrix of cp to cp connections.
            branch_start (int) [-]: Branch start busbar number.
            branch_end (int) [-]: Branch end busbar number.
            chain_step (int) [-]: Indicates device to be considered.
            z_array (list) [-]: Impedace of local connection.

        Returns:
            branch_data
            cp_ends
            branch_count

        '''

        cp_to_cp = copy.deepcopy(self.cp_to_cp)
        branch_end = max(cp_ends)+1
        #start from a cp connected to shore
        for connection in np.where(np.array(cp_ends) > 0)[0]:
    
            chain_step = connection
            branch_start = cp_ends[chain_step]
    
            while np.any(cp_to_cp[chain_step] > 0):
    
                next_cp = np.where(cp_to_cp[chain_step] > 0)[0]
    
                for cp in next_cp:
                    
                    if not cp_ends[cp] == 0: branch_end = cp_ends[cp]
                        
                    z_array = z_matrix[chain_step][cp]
    
                    branch_data[branch_count][:5] = [branch_start, branch_end, 
                                                     z_array[0], z_array[1],
                                                     z_array[2]]
                    cp_to_cp[chain_step][cp] = 0
                    cp_to_cp[cp][chain_step] = 0

                    cp_ends[cp] = branch_end
                    
                    self.array_branches.append(branch_end)

                    branch_end += 1
                    branch_count += 1

                chain_step = cp

        if not self.array_voltage == self.device_voltage:

            index = [x for x, val in enumerate(self.shore_to_cp) if val == 0]

            for cp in index:
                # add transformers to these collection points
                # 0 is for resistance, i.e. ideal transformer
                branch_start = cp_ends[cp]
                branch_data[branch_count][:4] = \
                    [branch_start, branch_end, 0, T_array_device]

                cp_ends[cp] = branch_end
                self.array_branches.append(branch_end)                
                
                branch_end += 1
                branch_count += 1

        return branch_data, cp_ends, branch_count

    def cp_to_device_to_device(self, branch_data, cp_ends, branch_count,
                               z_matrix):
        
        '''Connect from collection point to device (and then to device). Update
        branch data with these connections.

        Args:
            branch_data (np.array) [-]: PyPower branch data structure.
            cp_ends (list) [-]: List of busbar numbers at 'end' of cp, i.e. at
                point of connection to cable.
            branch_count (int) [-]: Global counter for branches.
            z_matrix (np.array) [pu]: Impedances expressed in matrix form.

        Attributes:
            device_to_device (np.array) [-]: Matrix of device to device
                connections.
            branch_start (int) [-]: Branch start busbar number.
            branch_end (int) [-]: Branch end busbar number.
            visited_nodes (list) [-]: List of nodes visited for ignore.
            cp (int) [-]: cp under consideration.
            chain_step (int) [-]: Indentifier of next device.

        Returns:
            branch_data

        '''

        device_to_device = copy.deepcopy(self.device_to_device)
        branch_end = max(cp_ends)+1
        visited_nodes = []
        static_ends = []

        for cp,_ in enumerate(self.cp_to_device):

            for connection in np.where(self.cp_to_device[cp] > 0)[0]:

                branch_start = cp_ends[cp]
                # find where this connects
                visited_nodes.append(connection)
                # then connect from here
                chain_step = connection
                device_to_device[:,connection] = 0

                # create branch
                if self.topology == 'Star':

                    z_device = z_matrix[cp][connection]

                elif self.substation:

                    # offset for radial substation
                    z_device = z_matrix[cp][connection + 1] 

                branch_data[branch_count][:5] = [branch_start, branch_end,
                                                 z_device[0], z_device[1],
                                                 z_device[2]]
                                                 
                self.array_branches.append(branch_end)

                static_ends.append(branch_end)
                branch_count += 1

                if chain_step < self.n_devices:

                    (branch_end, branch_start, visited_nodes, branch_count,
                     device_to_device, branch_data, static_ends) = \
                         self.device_to_device_fn(
                             device_to_device, chain_step, visited_nodes,
                             branch_start, branch_end, branch_data, 
                             branch_count, False, z_matrix, static_ends)

        # set branch constant values
        branch_data[: , [5, 6, 7, 10, 11, 12]] = [0, 0, 0, 1, -360, -360]

        return (branch_data, static_ends, visited_nodes, branch_count)

    def add_umbilicals(self, branch_data, static_ends, device_order,
                       branch_count, z_umbilical_all):

        branch_end = max(static_ends) + 1
        
        for device, cable_end in zip(device_order, static_ends):

            # add umbilical here
            z_umbilical = z_umbilical_all[device - 1]

            branch_start = cable_end

            z_umbilical = z_umbilical_all[device]

            branch_data[branch_count][:5] = [branch_start, branch_end,
                                             z_umbilical[0], z_umbilical[1],
                                             z_umbilical[2]]

            branch_count += 1
            branch_end += 1

        return branch_data

    def run_pf(self, power_factor, power):

        '''Run the power flow. This is repeated for each power bin contained in
        the power factor data structure.

        Args:
            power_factor (list) [-]: List of tuples, val1 = array power output
                in put, val2 = power factor.
            power (float) [W]: Rated power of oec.
        
        Attributes:
            ppopt (dict) [-]: PyPower options.
            ppc (dict) [-]: PyPower case definition. This contains all newtork
                information.
            result (dict) [-]: PyPower results.
            all_success (list) [int]: Binary flag indicating if power solved, 
                1 = successful solution, 0 = unsuccessful solution.
            all_active_powers (list) [MW]: Active power delivered onshore.
            all_reactive_powers (list) [MVAr]: Reactive power delivered
                onshore.
            all_busbar_voltage (list) [pu]: Busbar voltage magnitudes.
            all_busbar_angle (list) [deg]: Busbar voltage angles.
            
        Returns:

        Note:
            The power flow is repeated for each defined power level. The length
            of the lists above should equal this length.

        '''

        ppopt = ppoption.ppoption(VERBOSE=0, OUT_ALL=0)
        ppc = {"version": '2'}   
        ppc["baseMVA"] = 100.0
        ppc["bus"] = self.bus_data
        ppc["branch"] = self.branch_data
        all_active_powers = []
        all_reactive_powers = []
        all_success = []
        all_busbar_voltage = []
        all_busbar_angle = []
        gen_active_powers = []
        gen_reactive_powers = []

        # for each generator output
        for output_power in power_factor:
            gen_power = output_power[0] * power/1000000 # power in MW
            gen_reactive_power = gen_power * np.tan(np.arccos(output_power[1]))

            for gen in range(1, self.n_devices+1):

                self.gen_data[gen][1] = gen_power
                self.gen_data[gen][2] = gen_reactive_power

            ppc["gen"] = self.gen_data
            result, success = runpf.runpf(ppc, ppopt)
            all_active_powers.append(result['gen'][0][1])
            all_reactive_powers.append(result['gen'][0][2])
            all_success.append(success)
            all_busbar_voltage.append(result['bus'][:,7])
            all_busbar_angle.append(result['bus'][:,8])
            gen_active_powers.append(result['gen'][:,1][1:])
            gen_reactive_powers.append(result['gen'][:,2][1:])

        self.flag = all_success
        self.onshore_active_power = all_active_powers
        self.onshore_reactive_power = all_reactive_powers
        self.busbar_voltages = all_busbar_voltage
        self.busbar_angles = all_busbar_angle
        self.gen_active_powers = gen_active_powers
        self.gen_reactive_powers = gen_reactive_powers
        self.all_results = result

        return result
        
    def calculate_impedances(self, distance_matrix, impedance, local_system):

        '''Convert cable impedance of given network system into pu value for
        use in power flow. Valid local_systems are: 'device', 'array' and
        'export'.

        Args:
            distance_matrix (np.array) [m]: Distance between components.
            impedance (tuple) [ohm/km]: Impedance of connection.
            local_system (str) [-]: Network system under consideration.

        Attributes:
            impedance_matrix (np.array) [ohm]: Impedances expressed in matrix
                form.
            shape (tuple) [-]: Shape of impedance_matrix.
            impedance_base (float) [ohm]: Impedance base for the local_system.
            impedance_pu (np.array) [pu]: impedance_matrix as pu values.

        Returns:
            impedance_pu

        '''

        impedance_matrix = self.calculate_impedance_matrix(
            distance_matrix, impedance)

        shape = impedance_matrix.shape

        impedance_base = self.calculate_impedance_base(local_system)

        shunt_base = 1/impedance_base

        impedance_matrix_flat = impedance_matrix.flatten()
        
        impedance_matrix_shunt_corrected = \
            self.calculate_susceptance(impedance_matrix_flat, skip = True)
            
        
        impedance_pu = self.convert_to_pu(impedance_matrix_shunt_corrected,
                                          impedance_base, shunt_base)
        
        impedance_pu = impedance_pu.reshape(shape)
        
        return impedance_pu

    def convert_to_pu(
            self, impedance_matrix_shunt_corrected, impedance_base, shunt_base):
        
        '''Calculate per unit quantities and merge back into list.

        '''

        r_pu = impedance_matrix_shunt_corrected[::3]/impedance_base
        x_pu = impedance_matrix_shunt_corrected[1::3]/impedance_base
        b_pu = impedance_matrix_shunt_corrected[2::3]/shunt_base

        impedance_pu = \
            list(itertools.chain.from_iterable(itertools.izip(r_pu,x_pu,b_pu)))

        return np.asarray(impedance_pu)

    def calculate_impedance_matrix(self, distance_matrix, impedance):
        
        '''Cobmine distance matrix and cable impedance to create an impedance
        matrix.

        Args:
            distance_matrix (np.array) [m]: Distance between components.
            impedance (tuple) [Ohm/km]: Impedance of connection.

        Attributes:
            row (int) [-]: Number of rows in distance_matrix.
            col (int) [-]: Number of columns in distance_matrix.
            distance_list (np.array) [km]: distance_matrix expressed in 1D.
            impedance_list (list) [ohm]: List of tuples expressing total
                impedance of connections.
            local_z (tuple) [ohm]: Impedance of local connection.
            impedance_matrix (list) [ohm]: Impedances expressed in matrix form.

        Returns:
            np.array

        Note:
            Distances are converted here from m to km to be compatible with
            units of cable impedance (/km).

        '''

        row, col = distance_matrix.shape
       
        distance_list = distance_matrix.flatten()
        
        distance_list = distance_list/1000

        impedance_list = []
        
        for item in distance_list:

            local_z = []

            for val in impedance:
                
                vis = item * val
                local_z.append(vis)

            local_z = tuple(local_z)
            impedance_list.append(local_z)

        impedance_matrix = \
            [impedance_list[col*i : col*(i+1)] for i in range(row)]

        return np.asarray(impedance_matrix)

    def calculate_impedance_base(self, local_system):
        
        '''Calculate impedance of given network system. Valid local_systems
        are: 'device', 'array' and 'export'.

        Args:
            local_system (str) [-]: Network system under consideration.

        Attributes:
            MVAb (float) [VA]: System base.
            voltage (float) [V]: Voltage of the local_system.
            impedance_base (float) [ohm]: Impedance base for the local_system.

        Returns:
            impedance_base

        '''

        

        # equation here
        if local_system == 'device':

            voltage = self.device_voltage

        elif local_system == 'array':

            voltage = self.array_voltage

        elif local_system == 'export':

            voltage = self.export_voltage

        elif local_system == 'umbilical':

            voltage = self.device_voltage

        else:
            # error catch here
            pass

        impedance_base = np.square(voltage)/self.MVAb

        return impedance_base
        
    def calculate_export_impedance(self, distance, impedance):
        
        '''Convert export cable impedance into pu value for use in power flow.
        
        Args:
            distance (float) [m]: Export cable distance.
            impedance (tuple) [ohm]: Impedance of export cable.

        Attributes:
            impedance_base (float) [ohm]: Impedance base of the export system.
            export_impedance_pu (list) [pu]: Export system impedance expressed
                in pu.

        Returns:
            export_impedance_pu

        Note:
            Distances are converted here from m to km to be compatible with
            units of cable impedance (/km).

        '''

        impedance_base = self.calculate_impedance_base('export')
        shunt_base = 1 / impedance_base
        
        export_r_pu = ((distance / 1000) * impedance[0]) / impedance_base
        export_x_pu = ((distance / 1000) * impedance[1]) / impedance_base
        export_c = ((distance / 1000) * impedance[2])
        
        export_b = self.susceptance_formula(export_c)
        export_b_pu = export_b / shunt_base

        export_impedance_pu = [export_r_pu, export_x_pu, export_b_pu]
        
        return export_impedance_pu

    def calculate_umbilical_impedance(self, umbilical_impedance):
        
        '''Converts umbilical impedance from ohm to pu.

        '''

        impedance_base = self.calculate_impedance_base('umbilical')
        shunt_base = 1 / impedance_base

        z_umbilical = [(val[0]/impedance_base,
                        val[1]/impedance_base,
                        self.susceptance_formula(val[2])/shunt_base)
                        for val in umbilical_impedance]

        return z_umbilical

    def calculate_susceptance(self, impedance_list, skip):
        
        '''Calculate susceptance for all impedances.
        
        '''

        if skip == True:

            C_values = impedance_list[2::3]
            B_values = map(self.susceptance_formula, C_values)
            impedance_list[2::3] = B_values
            
            return impedance_list

        else:

            return map(self.susceptance_formula, impedance_list)

    def susceptance_formula(self, C):
        
        '''Calculate susceptance for pypower format.

        Args:
            C (float) [uF]: Capacitance.

        Attributes:
            B (float) [Mho]: Susceptance.
            pi_model (float) [Mho]: Susceptance for pi model.

        Returns:
            pi_model.

        Note:
            frequency is assumed as 50 Hz.

        '''

        f = 50.

        B = 2*np.pi*f*C*1e-6 # convert from uF to F

        pi_model = B/2

        return pi_model

    def transformer_impedance(self, voltage1, voltage2, transformer_db):

        '''Get transformer impedance and rating from db and convert to
        system_base.
        
        Args:
            voltage1 (float) [kV]: Primary winding voltage.
            voltage2 (float) [kV]: Secondary winding voltage.
            transformer_db (pd.DataFrame) [-]: Transformer database.
            
        Attributes:
            voltages (list) [V]: Voltages for pd.DataFrame search.
            transformer_base (float) [MVA]: Transformer base.
            transformer_pu_own_base (float) [pc]: Transformer series reactance
                in put to be converted.
            transformer_pu_system_base (float) [pu]: Transformer series
                reactance converted to system base.

        Returns:
            transformer_pu_system_base.
            
        Note:
            voltage1 and voltage2 are combined into voltages to avoid possible
            error with order of voltages.

        '''

        voltages = (voltage1 , voltage2 )
        # get impedance from db
        transformer_base, transformer_pu_own_base = \
            transformer_db[(transformer_db.v1.isin(voltages)) & \
                           (transformer_db.v2.isin(voltages))]\
                           [['rating', 'impedance']].values[0]

        transformer_pu_system_base = \
            transformer_pu_own_base /100. * (transformer_base/self.MVAb)

        return transformer_pu_system_base
        
class ComponentLoading(object):

    def __init__(self, system, voltage):

        self.flag = False
        self.system = system
        self.voltage = voltage

    def check_component_loading(self,
                                pf_result,
                                cable_rating,
                                nodes,
                                allow_overload = False):

        '''Get power flow in branch elements and compare against ratings from
        db. If loading exceeds branch rating, raise flag.

        '''

        bus_data = zip(pf_result['bus'][:, 0].astype(int),
                       pf_result['bus'][:, 9])

        P_injection = pf_result['branch'][: , 14] # results in MW
        Q_injection = pf_result['branch'][: , 15] # results in MVAr

        S_injection = [np.sqrt(np.square(P) + np.square(Q))
                       for P, Q
                       in zip(P_injection, Q_injection)] # result in MVA

        S_injection_branched = \
            zip(pf_result['branch'][:, 1].astype(int), S_injection)

        S_injection = [t[1] for t in S_injection_branched if t[0] in nodes]

        V = [self.voltage] * len(S_injection)

        I_injection = map(self.current_formula, S_injection, V)

        self.constraint_breach = [i for i in I_injection if i > cable_rating]

        if self.constraint_breach: 

            self.flag = True

        return

    def current_formula(self, S, V):
        
        '''Current injections from apparent power and voltage.
        
        Args:
            S (float) [MVA]: Apparent power.
            V (float) [V]: Voltage.

        Attributes:
            current (float) [A]: Current injection.

        Returns:
            current.       

        '''

        current = (S * 1e6) / (float(V) * 1000)

        return current