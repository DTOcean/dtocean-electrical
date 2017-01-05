# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 15:56:17 2016

@author: acollin
"""
import pandas as pd
import numpy as np
import pytest
import os
from shapely.geometry import Polygon

#from dtocean_electrical.hydro.hydro import emulate_wave_output
from dtocean_electrical.main import Electrical

from dtocean_electrical.inputs import (ElectricalComponentDatabase,
                    ElectricalMachineData,
                    ElectricalArrayData,
                    ConfigurationOptions,
                    ElectricalSiteData,
                    ElectricalExportData)

@pytest.mark.skipif(True,
                    reason="forced test to skip")
def test_main():
    ## Set test data directory
    #dir_ = 'C:/Users/acollin/Desktop/python/testdata/'
    this_dir = os.path.dirname(os.path.realpath(__file__))
    parent_dir = os.path.dirname(this_dir)
    dir_ = parent_dir + "\\sample_data\\"

    ## Get database
    file_ = 'mock_db_v2'
    xls_file = pd.ExcelFile(dir_ + file_ + '.xlsx', encoding = 'utf-8')
    sheet_names = xls_file.sheet_names
    static_cables = xls_file.parse(sheet_names[0])
    dynamic_cables = xls_file.parse(sheet_names[1])
    wet_mate = xls_file.parse(sheet_names[2])
    dry_mate = xls_file.parse(sheet_names[3])
    transformer = xls_file.parse(sheet_names[4])
    collection_point = xls_file.parse(sheet_names[5])
    switchgear = xls_file.parse(sheet_names[6])
    power_quality = xls_file.parse(sheet_names[7])

    database = ElectricalComponentDatabase(static_cables,
                                           dynamic_cables,
                                           wet_mate,
                                           dry_mate,
                                           transformer,
                                           collection_point,
                                           switchgear,
                                           power_quality)

    ## Define lease area
    # load bathy data from excel
    file_ = 'lease_area'
    xls_file = pd.ExcelFile(dir_ + file_ + '.xlsx', encoding = 'utf-8')
    sheet_names = xls_file.sheet_names
    lease_bathymetry = xls_file.parse(sheet_names[0])

    lease_exclusion_zones = [Polygon([(40.,40.),
                                      (40.,110.),
                                      (110.,110.),
                                      (110.,40.)])]

    lease_max_temp = 10.
    lease_max_soil_res = 10.
    lease_tidal_current_direction = 10.
    lease_tidal_current_flow = 10.
    lease_wave_direction = 10.
    lease_shipping = np.asarray([[1.,2.],[.5,.5]])

    site = ElectricalSiteData(lease_bathymetry,
                              lease_exclusion_zones,
                              lease_max_temp,
                              lease_max_soil_res,
                              lease_tidal_current_direction,
                              lease_tidal_current_flow,
                              lease_wave_direction,
                              lease_shipping)

    ## Define export area
    # load bathy data from excel         
    file_ = 'export_area'
    xls_file = pd.ExcelFile(dir_ + file_ + '.xlsx', encoding = 'utf-8')
    sheet_names = xls_file.sheet_names
    export_bathymetry = xls_file.parse(sheet_names[0])

    export_exclusion_zones = [Polygon([(40.,40.),
                                       (40.,110.),
                                       (110.,110.),
                                       (110.,40.)])]

    export_max_temp = 10.
    export_max_soil_res = 10.
    export_tidal_current_direction = 10.
    export_tidal_current_flow = 10.
    export_wave_direction = 10.
    export_shipping = np.array([[],[]])

    export = ElectricalExportData(export_bathymetry,
                                  export_exclusion_zones,
                                  export_max_temp,
                                  export_max_soil_res,
                                  export_tidal_current_direction,
                                  export_tidal_current_flow,
                                  export_wave_direction,
                                  export_shipping)

    ## Define the device
    technology = 'fixed'
    power = 1500000.
    voltage = 11000.
    connection = 'wet-mate'
    footprint_radius = None
    footprint_coords = [(0.0,25.0,0.), (-25.0,-25.0,0.), (25.,-25.,0.)]
    variable_power_factor = [[1,1],[0.5,1],[0.8,1]]
    constant_power_factor = None
    connection_point = (0,0,0)
    equilibrium_draft = 0.
    
    machine = ElectricalMachineData(technology,
                                    power,
                                    voltage,
                                    connection,
                                    variable_power_factor,
                                    constant_power_factor,
                                    footprint_radius,
                                    footprint_coords,
                                    connection_point,
                                    equilibrium_draft,
                                    )

    ## Define the array
    landing_point = (0.,0.,0.)
    layout = {'Device001': (100.,900.,0.),
              'Device002': (300.,900.,0.),
              'Device003': (400.,600.,0.)}
    
    number_of_devices = len(layout)
    array_output = [0.6, 0.3, 0.1]
    onshore_infrastructure_cost = 1000000.
    onshore_losses = 0.
    control_signal_type = 'fibre_optic'
    control_signal_cable = True
    control_signal_channels = 2
    voltage_limit_min = 0.9
    voltage_limit_max = 1.1
    offshore_reactive_limit = None
    orientation_angle = 0.

    array = ElectricalArrayData(machine,
                                landing_point,
                                layout,
                                number_of_devices,
                                array_output,
                                onshore_infrastructure_cost,
                                onshore_losses,
                                control_signal_type,
                                control_signal_cable,
                                control_signal_channels,
                                voltage_limit_min,
                                voltage_limit_max,
                                offshore_reactive_limit,
                                orientation_angle,
                                )
    
    ## Configuration options
    network_configuration = ['radial', 'star']
    export_voltage = None
    export_cables = None
    ac_power_flow = True
    target_burial_depth_array = None
    target_burial_depth_export = None
    connector_type = None
    collection_point_type = None
    devices_per_string = None
    equipment_gradient_constraint = 14.0
    
    # equipment soil installation compatibility matrix is loaded from file
    file_name = 'equipment_compatibility_matrix.xlsx'
    xls_file = pd.ExcelFile(os.path.join(dir_, file_name), encoding = 'utf-8')
    sheet_names = xls_file.sheet_names
    equipment_soil_compatibility = xls_file.parse(sheet_names[0])

    umbilical_safety_factor = 1.4925
    gravity = 9.880665
    
    options = ConfigurationOptions(network_configuration,
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
                                   umbilical_safety_factor,
                                   gravity,
                                   )

    ## Create object and run
    Elec = Electrical(site, array, export, options, database)
    Elec.run_module(plot=False)
    
    assert True

