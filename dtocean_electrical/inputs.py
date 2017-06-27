# -*- coding: utf-8 -*-
"""
This module defines the DTOcean electrical subsystems module inputs.

.. module:: input
   :platform: Windows
   :synopsis: Input module to DTOcean electrical subsystems module.
   
.. moduleauthor:: Adam Collin <a.collin@ed.ac.uk>
"""

# External module import
import numpy as np
from shapely.geometry import Polygon
from copy import deepcopy
from input_utils.utils import (seabed_range,
                               device_footprints_from_coords,
                               device_footprints_from_rad,
                               ideal_power_quantities,
                               get_bin_edges,
                               convert_df_column_type,
                               get_key,
                               set_burial_from_bpi)

class ElectricalComponentDatabase(object):

    '''Container class for the electrical component database.

    Args:
        static_cable (pd.DataFrame) [-]: static cable data;
            id (int) [-]: a unique key identifier.
            n (int) [-]: the number of conductors.
            v_rate (float) [V]: the rated voltage.
            a_air (float) [A]: the rated current operating in air.
            a_bury (float) [A]: the rated current operating buried.
            a_jtube (float) [A]: the rated current operating in jtube.
            r_dc (float) [Ohm/km]: dc resistance at 20 degree.
            r_ac (float) [Ohm/km]: ac resistance at 90 degree.
            xL (float) [Ohm/km]: inductive reactance per unit length.
            c (float) [uF/km]: capacitance per unit length.
            colour (string) [-]: serving colour.
            dry_mass (float) [kg/km]: cable weight in air per unit length.
            wet_mass (float) [kg/km]: cable weight in water per unit length.
            diameter (float) [m]: cable diameter.
            mbr (float) [m]: minimum bend radius.
            mbl (float) [N]: minimum break load.
            fibre (Bool) [-]: fibre optic cable. True = yes, False = no.
            cost (float) [E/m]: unit cost per metre.
            max_operating_temp (float) [deg]: maximum temperature.
            environmental_impact (X) [X]: data yet to be formatted.

        dynamic_cable (pd.DataFrame) [-]: dynamic cable data;
            id (int) [-]: a unique key identifier.
            n (int) [-]: the number of conductors.
            v_rate (float) [V]: the rated voltage.
            a_air (float) [A]: the rated current operating in air.
            a_bury (float) [A]: the rated current operating buried.
            a_jtube (float) [A]: the rated current operating in jtube .
            r_dc (float) [Ohm/km]: dc resistance at 20 degree.
            r_ac (float) [Ohm/km]: ac resistance at 90 degree.
            xL (float) [Ohm/km]: inductive reactance per unit length.
            c (float) [uF/km]: capacitance per unit length.
            colour (string) [-]: serving colour.
            dry_mass (float) [kg/km]: cable weight in air per unit length.
            wet_mass (float) [kg/km]: cable weight in water per unit length.
            diameter (float) [m]: cable diameter.
            mbr (float) [m]: minimum bend radius.
            mbl (float) [N]: minimum break load.
            fibre (Bool) [-]: fibre optic cable. True = yes, False = no.
            cost (float) [E/m]: unit cost per metre.
            max_operating_temp (float) [deg]: maximum temperature.
            environmental_impact (X) [X]: data yet to be formatted.

        wet_mate_connectors (pd.DataFrame) [-]: wet mate connector data;
            id (int) [-]: a unique key identifier.
            n (int) [-]: the number of contacts.
            v_rate (float) [V]: the rated voltage.
            a_rate (float) [A]: the rated current.
            dry_mass (float) [kg/km]: unit weight in air.
            height (float) [m]: unit height.
            width (float) [m]: unit width.
            depth (float) [m]: unit depth.
            mating (float) [N]: force required to make the connection.
            demating (float) [N]: force required to make the disconnection.
            colour (string) [-]: outer colour of the unit.
            outer_coating (string) [-]: material of the connector housing.
            fibre (int) [-]: the number of fibre optical channels. Can be zero.
            cost (float) [E]: unit cost.
            max_water_depth (float) [m]: the maximum specified operating depth.
            min_temperature (float) [degree celcius]: minimum operational 
                temperature. 
            max_temperature (float) [degree celcius]: maximum operational 
                temperature.
            min_cable_csa (float) [mm^2]: Minimum cable size compatibility.
            max_cable_csa (float) [mm^2]: Maximum cable size compatibility.

        dry_mate_connectors (pd.DataFrame) [-]: dry mate connector data;
            id (int) [-]: a unique key identifier.
            n (int) [-]: the number of contacts.
            v_rate (float) [V]: the rated voltage.
            a_rate (float) [A]: the rated current.
            dry_mass (float) [kg/km]: unit weight in air.
            height (float) [m]: unit height.
            width (float) [m]: unit width.
            depth (float) [m]: unit depth.
            mating (float) [N]: force required to make the connection.
            demating (float) [N]: force required to make the disconnection.
            colour (string) [-]: outer colour of the unit.
            outer_coating (string) [-]: material of the connector housing.
            fibre (int) [-]: the number of fibre optical channels. Can be zero.
            cost (float) [E]: unit cost.
            max_water_depth (float) [m]: the maximum specified operating depth.
            min_temperature (float) [degree celcius]: minimum operational 
                temperature. 
            max_temperature (float) [degree celcius]: maximum operational 
                temperature.
            min_cable_csa (float) [mm^2]: Minimum cable size compatibility.
            max_cable_csa (float) [mm^2]: Maximum cable size compatibility.

        transformers (pd.DataFrame) [-]: power transformer data;
            id (int) [-]: a unique key identifier.
            windings (int) [-]: the number of windings.
            rating (float) [VA]: transformer power rating.
            v1 (float) [V]: rated voltage of the primary winding.
            v2 (float) [V]: rated voltage of the secondary winding.
            v3 (float) [V]: rated voltage of the tertiary winding.
            dry_mass (float) [kg/km]: unit weight in air.
            height (float) [m]: unit height.
            width (float) [m]: unit width.
            depth (float) [m]: unit depth.
            colour (string) [-]: outer colour of the unit.
            outer_coating (string) [-]: material of the transformer housing.
            cost (float) [E]: unit cost.
            cooling (string) [-]: type of cooling system.
            max_water_depth (float) [m]: the maximum specified operating depth.
            min_temperature (float) [degree celcius]: minimum operational 
                temperature. 
            max_temperature (float) [degree celcius]: maximum operational 
                temperature.
            impedance (float) [pc]: transformer series reactance, given to own
                base (specified by rating).

        collection_points (pd.DataFrame) [-]: collection point data;
            id (int) [-]: a unique key identifier.
            input (int) [-]: number of input connections.
            output (int) [-]: number of output connections.
            input_connector (str) [-]: type of input connector(s).
            output_connector (str) [-]: type of output connector(s).
            v1 (float) [V]: rated voltage of the primary winding.
            v2 (float) [V]: rated voltage of the secondary winding.
            a_rate (float) [A]: the rated current.
            dry_mass (float) [kg/km]: unit weight in air.
            wet_mass (float) [kg/km]: unit weight in water.
            height (float) [m]: unit height.
            width (float) [m]: unit width.
            depth (float) [m]: unit depth.
            colour (string) [-]: outer colour of the unit.
            outer_coating (string) [-]: material of the unit housing.
            foundation (str) [-]: foundation type.
            busbar (str) [-]: busbar configuration.
            cost (float) [E]: unit cost.
            cooling (string) [-]: type of cooling system.
            max_water_depth (float) [m]: the maximum specified operating depth.
            fibre (int) [-]: the number of fibre optical channels. Can be zero.
            gravity_centre (float) [m]: the unit centre of gravity.
            operating_environment (str) [-]: defines if the unit is prepared
                for subsea, offshore or onshore operation.
            min_temperature (float) [degree celcius]: minimum operational 
                temperature. 
            max_temperature (float) [degree celcius]: maximum operational 
                temperature.
            wet_frontal_area (float) [m^2]: Unit wet frontal area.
            dry_frontal_area (float) [m^2]: Unit dry frontal area.
            wet_beam_area (float) [m^2]: Unit wet beam area.
            dry_beam_area (float) [m^2]: Unit dry beam area.
            orientation_angle (float) [degree]: Unit orientation angle.
            foundation_loc (list(tuple)) [m]: Foundation locations in local
                coordinates for N foundations.

        switchgear (pd.DataFrame) [-]: switchgear equipment class data;
            id (int) [-]: a unique key identifier.
            v_rate (float) [V]: the rated voltage.
            a_rate (float) [A]: the rated current.
            dry_mass (float) [kg/km]: unit weight in air.
            height (float) [m]: unit height.
            width (float) [m]: unit width.
            depth (float) [m]: unit depth.
            max_water_depth (float) [m]: the maximum specified operating depth.
            operating_environment (str) [-]: defines if the unit is prepared
                for subsea, offshore or onshore operation.
            outer_coating (string) [-]: material of the unit housing.
            colour (string) [-]: outer colour of the unit.
            min_temperature (float) [degree celcius]: minimum operational 
                temperature. 
            max_temperature (float) [degree celcius]: maximum operational 
                temperature.
            cost (float) [E]: unit cost.

        power_quality (pd.DataFrame) [-]: power quality equipment class data;
            id (int) [-]: a unique key identifier.
            v_rate (float) [V]: the rated voltage
            reactive_power (float) [VAr]:reactive power capability of the unit.
            n_control (int) [-]: the number of discrete control steps.
            var_per_step {float) [-]: the reactive power per discrete step.
            dry_mass (float) [kg/km]: unit weight in air.
            operating_environment (str) [-]: defines if the unit is prepared
                for subsea, offshore or onshore operation.
            height (float) [m]: unit height.
            width (float) [m]: unit width.
            depth (float) [m]: unit depth.
            cooling (string) [-]: type of cooling system.
            outer_coating (string) [-]: material of the unit housing
            colour (string) [-]: outer colour of the unit
            min_temperature (float) [degree celcius]: minimum operational 
                temperature. 
            max_temperature (float) [degree celcius]: maximum operational 
                temperature.
            max_water_depth (float) [m]: the maximum specified operating depth.
            cost (float) [E]: unit cost.

    Attributes:
        static_cable (pd.DataFrame)
        dynamic_cable (pd.DataFrame)
        wet_mate_connectors (pd.DataFrame)
        dry_mate_connectors (pd.DataFrame)
        transformers (pd.DataFrame)
        collection_points (pd.DataFrame)
        switchgear (pd.DataFrame)
        power_quality (pd.DataFrame)

    Returns:

    '''

    def __init__(self,
                 static_cable,
                 dynamic_cable,
                 wet_mate_connectors,
                 dry_mate_connectors,
                 transformers,
                 collection_points,
                 switchgear,
                 power_quality):

            self.static_cable = static_cable
            self.dynamic_cable = dynamic_cable
            self.wet_mate_connectors = wet_mate_connectors
            self.dry_mate_connectors = dry_mate_connectors
            self.transformers = transformers
            self.collection_points = collection_points
            self.switchgear = switchgear
            self.power_quality = power_quality

class ElectricalSiteData(object):

    '''Define the electrical systems site data object. This includes all
    geotechnical and geophysical data.
    
    Args:
        bathymetry (pd.dataframe) [m]: the vertical profile of the sea bottom
            and geotechnical layers at each (given) UTM coordinate within the
            lease area; expressed as [id,i,j,x,y,layer1depth,layer1type,..,..,
            layerNdepth, layerNtype], where: i and j are local indices; x,y are
            grid coordinates and layerndepth, laterntype define the start depth
            and the type of layer n (from n=1:N).
        exclusion_zones (list) [-]: list containing the UTM coordinates of the
                                    exclusion zone polygons within the lease
                                    area.
        max_temp (float) [degree]: the maximum seabed temperature recorded in
                                   the lease area.
        max_soil_res (float) [K.m/W]: the maximum soil resistivity recorded in
                                      the lease area.
        tidal_current_direction (float) [rad]: the tidal stream current
                                               prevalent direction in the lease
                                               area.
        tidal_current_flow (float) [m/s]: maximum tidal stream current velocity
                                          in the lease area.
        wave_direction (float) [rad]: the prevalent wave direction in the lease
                                      area.
        shipping (numpy.ndarray) [-]: histogram of shipping activity in the
                                      lease area; expressed as [val1, val2]
                                      where: val1 is the bin edge of the vessel
                                      deadweight (T) val2 is the frequency (pc)

    Attributes:
        bathymetry (pd.dataframe)           
        exclusion_zones (list)
        max_temp (float)
        max_soil_res (float)
        tidal_current_direction (float)
        tidal_current_flow (float)
        wave_direction (float)
        shipping (numpy.ndarray)
        max_water_depth (float) [m]: max water depth obtained from bathy data.
        min_water_depth (float) [m]: max water depth obtained from bathy data.

    Returns:

    '''

    def __init__(self,
                 bathymetry,
                 exclusion_zones,
                 max_temp,
                 max_soil_res,
                 tidal_current_direction,
                 tidal_current_flow,
                 wave_direction,
                 shipping):

        self.bathymetry = self.check_bathy_data(bathymetry)
        self.exclusion_zones = exclusion_zones
        self.max_temp = max_temp
        self.max_soil_res = max_soil_res
        self.tidal_current_direction = tidal_current_direction
        self.tidal_current_flow = tidal_current_flow
        self.wave_direction = wave_direction
        self.shipping = shipping
        self.min_water_depth, self.max_water_depth = \
            seabed_range(self.bathymetry)

        self.set_target_burial_depth()

    def check_bathy_data(self, bathymetry):

        '''Assign bathymetry and convert column type of i and j.

        Args:
            bathymetry

        Returns:
            (DataFrame) [-]: Updated DataFrame

        '''

        return convert_df_column_type(bathymetry, ['i', 'j'], np.int64)
    
          
    def set_target_burial_depth(self):
        
        self.bathymetry[u'Target burial depth'] = self.bathymetry.apply(
            set_burial_from_bpi, axis = 1)
        
        return

class ElectricalExportData(object):
    
    '''Define the electrical systems export data object. This includes all
    geotechnical and geophysical data.

    Args:
        bathymetry (pd.dataframe) [m]: the vertical profile of the sea bottom
            and geotechnical layers at each (given) UTM coordinate within the
            export area; expressed as [id,i,j,x,y,layer1depth,layer1type,..,..,
            layerNdepth, layerNtype], where: i and j are local indices; x,y are
            grid coordinates and layerndepth, laterntype define the start depth
            and the type of layer n (from n=1:N).
        exclusion_zones (list) [-]: list containing the UTM coordinates of the
            exclusion zone polygons within the export area.
        max_temp (float) [degree]: the maximum seabed temperature recorded in
            the export area.
        max_soil_res (float) [K.m/W]: the maximum soil resistivity recorded in
            the export area.
        tidal_current_direction (float) [rad]: the tidal stream current
            prevalent direction in the export area.
        tidal_current_flow (float) [m/s]: maximum tidal stream current velocity
            in the export area.
        wave_direction (float) [rad]: the prevalent wave direction in the
            export area.
        shipping (numpy.ndarray) [-]: histogram of shipping activity in the
            export area; expressed as [val1, val2] where: val1 is the bin edge
            of the vessel deadweight (T) val2 is the frequency (pc).

    Attributes:
        bathymetry (pd.dataframe)
        exclusion_zones (list) 
        max_temp (float)
        max_soil_res (float)
        tidal_current_direction (float)
        tidal_current_flow (float)
        wave_direction (float)
        shipping (numpy.ndarray)
        max_water_depth (float) [m]: max water depth obtained from bathy data.
        min_water_depth (float) [m]: max water depth obtained from bathy data.

    Returns:

    '''

    def __init__(self,
                 bathymetry,
                 exclusion_zones,
                 max_temp,
                 max_soil_res,
                 tidal_current_direction,
                 tidal_current_flow,
                 wave_direction,
                 shipping):

        self.bathymetry = self.check_bathy_data(bathymetry)
        self.exclusion_zones = exclusion_zones
        self.max_temp = max_temp
        self.max_soil_res = max_soil_res
        self.tidal_current_direction = tidal_current_direction
        self.tidal_current_flow = tidal_current_flow
        self.wave_direction = wave_direction
        self.shipping = shipping
        self.min_water_depth, self.max_water_depth = \
            seabed_range(self.bathymetry)
            
        self.set_target_burial_depth()

    def check_bathy_data(self, bathymetry):
        
        '''Assign bathymetry and convert column type of i and j.
        
        Args:
            bathymetry.
            
         Returns:
            (DataFrame) [-]: Updated DataFrame
        
        '''

        return convert_df_column_type(bathymetry, ['i', 'j'], np.int64)


    def set_target_burial_depth(self):
        
        self.bathymetry[u'Target burial depth'] = self.bathymetry.apply(
            set_burial_from_bpi, axis = 1)
        
        return

class ElectricalMachineData(object):

    '''Container class to carry the OEC device object.

    Args:
        technology (str) [-]: floating or fixed
        power (float) [W]: OEC rated power output
        voltage (float) [V]: OEC rated voltage at point of network connection
        connection (str) [-]: Type of connection, either 'wet-mate', 'dry-mate'
            or 'hard-wired'.
        variable_power_factor (list) [-]: List of tuples for OEC power factor;
                                 val1 = power in pu, val2 = pf.
        constant_power_factor (float) [-]: A power factor value to be applied
            at every point of analysis.
        footprint_radius (float) [m]: The device footprint defined by radius.
        footprint_coords (list) [m]: The device footprint by utm [x,y,z]
            coordinates.
        connection_point (tuple) [m]: Location of electrical connection, as
            (x, y, z) coordinates in local coordinate system.
        equilibrium_draft (float) [m]: Device equilibrium draft without mooring
            system.

    Attributes:
        technology (int)
        power (float)
        voltage (float)
        connection (int)
        power_factor (list): OEC power factor;
                                 val1 = power in pu,
                                 val2 = pf,
                                 val3 = pf_angle
        footprint_type (string): Either 'radius' or 'coords'.
        footprint (-) [m]: Value of footprint_radius or footprint_coords.
        connection_point (tuple)
        draft (float) [m]: See equilibrium_draft.

    Returns:
        None

    Notes:
        footprint_radius and footprint_coords are conflicting data and only one
        should be supplied. A logic test is included.

    '''

    def __init__(self,
                 technology,
                 power,
                 voltage,
                 connection,
                 variable_power_factor,
                 constant_power_factor,
                 footprint_radius,
                 footprint_coords,
                 connection_point,
                 equilibrium_draft,
                 ):

        self.technology = technology
        self.power = power
        self.voltage = voltage
        self.connection = connection
        
        self.power_factor = self._check_power_factor_type(
            constant_power_factor, variable_power_factor)

        self.footprint_type, self.footprint = self._check_footprint_type(
            footprint_radius, footprint_coords)
        
        self.connection_point = connection_point
        self.draft = equilibrium_draft
        self.floating = self._set_floating_flag()
        
        self.max_current = self.get_current_ratings(constant_power_factor,
                                                    variable_power_factor)

    def _check_power_factor_type(self, constant, variable):
        
        '''Process all power factor data. This checks to see what values have
        been supplied.
        
        Args:
            constant (float) [-]: Single power factor value to be applied at
                all points of analysis.
            variable (list) [-]: List of power factors and the corresponding
                power value.

        Note:
            If no data is supplied a constant unity power factor is assumed.
            If both are supplied then the variable values are assumed.

        '''

        if variable is not None:

            power_factor = variable

        elif constant is not None:

            power_factor = constant
            
        else:

            power_factor = 1.0

        return power_factor

    def _calculate_power_factor_angle(self, power_factor):

        '''This calculates the power factor angle in radians for each specified
        power factor.
        
        Args:
            power_factor (list) [-]: OEC power factor;
                                     val1 = power in pu,
                                     val2 = pf,
                                     val3 = pf_angle
        
        Attributes:
            power_factor (list)
                                 
        Returns:
            power_factor

        '''

        j = []
        for i in power_factor:
            
            power_factor_angle = np.arccos(i[1])
            
            i.append(power_factor_angle)
            
            j.append(i)

        power_factor = j

        return power_factor

    def _check_footprint_type(self, radius, coordinates):

        '''Logic test against footprint_radius and footprint_coords.

        Args:
            radius (float) [m]: The device footprint defined by radius.
            coordinates (list) [m]: The device footprint by utm [x,y,z]
                coordinates.

        Attributes:
            area_radius (float) [m2]:
            coordinate_radius (float) [m2]:

        Returns:
            type_ (str) [-]: The radius type., either 'radius' or 'footprint'.
            value (-) [m]: The radius value.

        Note:
            Decision to be taken on what to do if both or none are defined. At
            the moment the largest area is used if both are defined. If the
            areas are equal, coordinates are used.
            None set to coordinates to avoid crashing the code.

        '''

        if radius and coordinates:

            area_radius = np.pi*(np.square(radius))

            coordinates_radius = Polygon(coordinates).area

            if area_radius > coordinates_radius:

                type_ = 'radius'
                value = radius

            else:

                type_ = 'coordinates'
                value= coordinates

        elif radius:

            type_ = 'radius'
            value = radius

        elif coordinates:

            type_ = 'coordinates'
            value = coordinates

        else:

            type_ = 'coordinates'
            value= coordinates

        return type_, value
    
    def _set_floating_flag(self):
        
        '''Check for floating device.

        Args:
            technology (str) [-]: Device technology.

        Attributes:
            flag (bool): True for floating device, False for fixed.

        Return:
            flag.

        '''

        if 'floating' in self.technology.lower():

            flag = True

        else:

            flag = False

        return flag

    def get_current_ratings(self, constant, variable):

        ''' Get current rating for each power factor to identify worst case
        loading condition.

        '''

        current = []
        
        if variable is not None:

            for power, pf in self.power_factor:
    
                current.append((power*self.power) / (self.voltage * pf))
                max_current = max(current)

        else:
            
            max_current = (1.*self.power) / (self.voltage * self.power_factor)

        return max_current

class ElectricalArrayData(object):

    '''Container class to carry the array object. This inherets the machine.

    Args:
        ElectricalMachineData (class) [-]: Class containing the machine
                                           specification
        landing_point (tuple) [m]: UTM coordinates of the landing areas;
                                   expressed as [x,y, id]
        layout (dict) [m]: OEC layout in dictionary from WP2;
                           key = device id,
                           value = UTM coordinates, as [x,y,z]
        n_devices (int) [-]: The number of OECs in the array.
        total_power (float) [W]: Installed power of the array.
        array_output (numpy.ndarray) [pc]: The total array power output in
            histogram form.
        control_signal_type (str) [-]: The type of control signal used in the
                                       array, accepts 'fibre optic'.
        control_signal_cable (bool) [-]: Defines if the control signal is to 
                                         packaged in the power cable (True) or
                                         not (False)
        control_signal_channels (int) [-]: Defines the number of control signal
                                           pairs per device
        voltage_limit_min (float) [pu]: The minimum voltage allowed in the
                                        offshore network
        voltage_limit_max (float) [pu]: The maximum voltage allowed in the
                                        offshore network
        offshore_reactive_limit (list) [-]: The target power factor at the
                                            offshore collection point. This is
                                            a list of pairs: val1 = power [pu],
                                            val2 = reactive power [pu]
        onshore_infrastructure_cost (float) [E]: Cost of the onshore
            infrastructure, for use in LCOE calculation.
        onshore_losses (float) [pc]: Electrical losses of the onshore
            infrastructure, entered as percentage of annual energy yield.
        orientation_angle (float) [degree]: Device orientation angle.

    Attributes:
        machine_data (Class)
        landing_point (tuple)
        layout (dict)
        n_devices (int)
        total_power (float)
        array_output (numpy.ndarray)
        onshore_infrastructure_cost (float)
        onshore_losses (float)
        control_signal_type (str)
        control_signal_cable (Bool)
        control_signal_channels (int)
        voltage_limit_min (float)
        voltage_limit_max (float)
        offshore_reactive_limit (list)
        orientation_angle (float)

    Returns:
        none

    Notes:
        The device footprints are also set here at an array level. The bin
        width of array_output is assumed equal.

    '''

    def __init__(self,
                 ElectricalMachineData,
                 landing_point,
                 layout,
                 n_devices,
                 array_output,
                 onshore_infrastructure_cost=0.,
                 onshore_losses=0.,
                 control_signal_type='fibre_optic',
                 control_signal_cable=True,
                 control_signal_channels=2,
                 voltage_limit_min=0.9,
                 voltage_limit_max=1.1,
                 offshore_reactive_limit=None,
                 orientation_angle=0.):

        self.machine_data = ElectricalMachineData
        self.landing_point = landing_point
        self.layout = layout
        self.n_devices = n_devices
        self.total_power = n_devices * ElectricalMachineData.power
        self.array_output = array_output        
        self.onshore_infrastructure_cost = onshore_infrastructure_cost
        self.onshore_losses = onshore_losses

        self.device_footprint = self._set_footprints()

        self.control_signal_type = control_signal_type
        self.control_signal_cable = control_signal_cable
        self.control_signal_channels = control_signal_channels
        self.voltage_limit_min = voltage_limit_min
        self.voltage_limit_max = voltage_limit_max
        self.offshore_reactive_limit = offshore_reactive_limit

        # set power factor based on array output bin edges
        self._set_histogram_edges(ElectricalMachineData.power_factor,
                                  array_output)

        # calculate ideal power quantities
        self.ideal_annual_yield, self.ideal_histogram = ideal_power_quantities(
            array_output, n_devices, ElectricalMachineData.power)
            
        self.orientation_angle = orientation_angle

    def _set_footprints(self):
        
        '''Set footprint of each device in the array.
        
        Attributes:
            device_footprint (list): footprint area of each device.
        
        Returns:
            list of footprint areas

        '''

        if self.machine_data.footprint_type is 'radius':
            device_footprint = (
                device_footprints_from_rad(
                    self.layout,self.machine_data.footprint))
        elif self.machine_data.footprint_type is 'coordinates':
            device_footprint = (
                device_footprints_from_coords(
                    self.layout,self.machine_data.footprint))

        return device_footprint
        
    def _set_histogram_edges(self, power_factor, array_output):

        '''Set the histogram edges depending on if a variable power factor is
        passed or not. This also updates the power_factor edges in the constant
        power factor case to provide the tuple pairs required for power flow
        analysis.

        Args:
            power_factor (unknown) [-]: Either single float or list of tuple
                pairs, where: val1 = device power, val2 = power factor.
            array_output (list) [-]: The total array power output in
                histogram form.

        Attributes:
            bins (int) [-]: The number of power bins in the histogram.
            edges (list) [pu]: List of array power outputs for analysis
                expressed as histogram bin edges.
            pf_range (tuple) [-]: Sets a range of powers to which a power
                factor is applied, where: val1 = bin min, val2 = bin max,
                val3 = power factor of bin.
            power_factor_range (list) [-]: List of pf_ranges.
            modified_power_factor (list) [-]: New list of power factors mapped
                to length and powers of array_output.

        Returns:

        '''

        bin_width = 1.0 / len(array_output)
        bin_offset = bin_width / 2
        centers = np.linspace(bin_offset, 1 - bin_offset, len(array_output))

        if type(power_factor) == list:
            
            power_factor = sorted(power_factor, key=get_key, reverse = True)
            power_factor_range = []
            # make range
            for i, item in enumerate(power_factor):
                if i != (len(power_factor)-1):
                    pf_range = (power_factor[i+1][0], item[0],  item[1])
                    power_factor_range.append(pf_range)
                else:
                    if item[0] == 0:
                        break
                    else:
                        pf_range = (0, item[0],  item[1])
                        power_factor_range.append(pf_range)

            modified_power_factor = []
            # apply range
            for c in centers:
                for range_ in power_factor_range:
                    if range_[0]< c <= range_[1]:
                        modified_power_factor.append((c, range_[2]))
                        break
                    
            self.machine_data.power_factor = modified_power_factor

        else:

            self.machine_data.power_factor = (
                zip(centers, [power_factor]*len(array_output))
                )

        return

class ConfigurationOptions(object):

    '''Container class for the configuration options defined in the core. These
    can be specificed by the user at GUI interface or by the core during the
    global optimisation process.

    Args:
        network_configuration (list, str) [-]: list of networks to evaluate:
            radial or star.
        export_voltage (float) [V]: export cable voltage.
        export_cables (int) [-]: number of export cables.
        ac_power_flow (Bool) [-]: run full ac power flow (True) or dc (False).
        target_burial_depth_array (float) [m]: array cable burial depth.
        target_burial_depth_export (float) [m]: export cable burial depth.
        connector_type (string) [-]: 'wet mate' or 'dry mate'. This will be
            applied to all connectors.
        collection_point_type (string) [-]: 'subsea' or 'surface'. This will be
            applied to all collection points.
        devices_per_string (int) [-]: number of devices per string.
        equipment_gradient_constraint (float) [deg]: the maximum seabed
            gradient considered by the cable routing analysis.
        equipment_soil_compatibility (pd.DataFrame) [m/s]: the equipment soil
            installation compatibility matrix.

###
        export_cable_length (float) [m]: user defined export cable length.
        max_water_depth (float) [m]: user defined maximum water depth.
        export_cable_protection (string) [-]: user defined cable protection
            technique, from: concrete_mattress, rock_bag, none.
        export_cable_derate (int) [%]: derating factor of export cable.
###

        umbilical_safety_factor (float) [-]: Umbilical safety factor from
            DNV-RP-F401.
        gravity (float) [m/s^2]: Acceleration due to gravity.
        user_umbilical  (int) [-]: Database reference of the user predefined
            umbilical.
        installation_tool (str) [-]: Predefined installation tool.

    Attributes:
        network_configuration (list, str)
        export_voltage (float)
        export_cables (int)
        ac_power_flow (Bool)
        target_burial_depth_array (float)
        target_burial_depth_export (float)
        connector_type (string)
        collection_point_type (string)
        devices_per_string (int)
        equipment_gradient_constraint (float)
        equipment_soil_compatibility (pd.DataFrame)

    Returns:

    Notes:
        Equipment_soil_compatibility is also used by the logistics and
        installations modules.

    '''
    
    def __init__(self,
                 network_configuration,
                 export_voltage,
                 export_cables,
                 ac_power_flow,
                 target_burial_depth_array,
                 target_burial_depth_export,
                 connector_type,
                 collection_point_type,
                 devices_per_string,
                 equipment_gradient_constraint,
                 equipment_soil_compatibility,
                 installation_tool = None,
                 umbilical_safety_factor = 1.4925,                 
                 gravity = 9.80655,
                 user_umbilical = None
                 ):

        self.network_configuration = network_configuration
        self.export_voltage = export_voltage
        self.export_cables = export_cables
        self.ac_power_flow = ac_power_flow
        self.target_burial_depth_array = target_burial_depth_array
        self.target_burial_depth_export = target_burial_depth_export
        self.connector_type = connector_type
        self.collection_point_type = collection_point_type
        self.devices_per_string = devices_per_string
        self.equipment_gradient_constraint = equipment_gradient_constraint
        self.installation_rates = equipment_soil_compatibility
        self.equipment_soil_compatibility = self.binary_compatibility_matrix(
            deepcopy(equipment_soil_compatibility))
        self.equipment_soil_compatibility_dict = \
            self.dict_compatibility_matrix(
                deepcopy(equipment_soil_compatibility))
        self.installation_tool = installation_tool
        self.umbilical_safety_factor = umbilical_safety_factor
        self.gravity = gravity
        self.user_umbilical = user_umbilical

    def binary_compatibility_matrix(self, matrix):

        '''Convert the equipment-soil compatibility matrix into binary.

        Args:
            matrix (pd.DataFrame) [-]: The equipment-soil compatibility matrix.

        Return:
            matrix (pd.DataFrame) [-]: The binary matrix.

        '''

        matrix[matrix > 0] = 1

        return matrix

    def dict_compatibility_matrix(self, matrix):

        '''Make a dictionary of soil types associated with each installation
        equipment.

        Args:
            matrix (pd.DataFrame) [-]:

        Attributes:
            soil_types (list) [-]:
            equipment_dictionary (dict) [-]:

        Returns:
            equipment_dictionary

        '''

        soil_types = matrix.columns.tolist()
        equipment_dictionary = {}

        for technique in matrix.index.values:

            equipment_dictionary[str(technique)] = (
                self.get_equipment_install_lists(soil_types, technique,matrix))

        return equipment_dictionary

    def get_equipment_install_lists(self, soil_types, technique, matrix):

        '''

        '''

        valid_soils = []
        for soil in soil_types :
            if matrix.loc[technique][soil] > 0:
                valid_soils.append(str(soil))

        return valid_soils
