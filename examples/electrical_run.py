# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 15:56:17 2016

@author: acollin
"""
import os

import pandas as pd
import numpy as np
from shapely.geometry import Polygon
import matplotlib.pyplot as plt

from dtocean_electrical import start_logging
#from dtocean_electrical.hydro.hydro import emulate_wave_output
from dtocean_electrical.inputs import (ElectricalComponentDatabase,
                                       ElectricalMachineData,
                                       ElectricalArrayData,
                                       ConfigurationOptions,
                                       ElectricalSiteData,
                                       ElectricalExportData)

from dtocean_electrical.main import Electrical

# Start the logging system
start_logging()

mod_path = os.path.realpath(__file__)
mod_dir = os.path.dirname(mod_path)

## Set test data directory
data_dir = os.path.join(mod_dir, "..", "sample_data")

## Get DB
file_name = 'mock_db_v2.xlsx'
xls_file = pd.ExcelFile(os.path.join(data_dir, file_name), encoding = 'utf-8')
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
                                       static_cables,
                                       dynamic_cables,
                                       wet_mate,
                                       dry_mate,
                                       transformer,
                                       collection_point,
                                       switchgear,
                                       power_quality)

## Define lease area
# load bathy data from excel
file_name = 'lease_area_0709.xlsx'
xls_file = pd.ExcelFile(os.path.join(data_dir, file_name), encoding = 'utf-8')
sheet_names = xls_file.sheet_names
lease_bathymetry = xls_file.parse(sheet_names[0])

lease_exclusion_zones = [Polygon([(40.,40.),(40.,110.),(110.,110.),(110.,40.)])]

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
file_name = 'export_area_0709.xlsx'
xls_file = pd.ExcelFile(os.path.join(data_dir, file_name), encoding = 'utf-8')
sheet_names = xls_file.sheet_names
export_bathymetry = xls_file.parse(sheet_names[0])

export_exclusion_zones = [Polygon([(40.,40.),(40.,110.),(110.,110.),(110.,40.)])]

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
                              export_shipping,
                              )

## Define the device
technology = 'fixed'
power = 1000000.
voltage = 11000.
connection = 'wet-mate'
footprint_radius = None
footprint_coords = [(0.0,25.0,0.), (-25.0,-25.0,0.), (25.,-25.,0.)]
#variable_power_factor = [[1,1],[0.5,1],[0.8,1]]
variable_power_factor = [[ 0.5,  0.8], [ 1. ,  1. ]]
constant_power_factor = 0.98
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
landing_point = (0.0, 1250.0)

#layout = {'Device003': (1280.618258923915, 687.5601863780051),
#          'Device002': (1200.0, 1100.0),
#          'Device001': (1119.381741076085, 1512.4398136219947),
#          'Device005': (1038.7634821521701, 1924.8796272439897), 
#          'Device004': (1498.0557936438454, 1694.6810423392592)}

layout = {'Device003': (1079.8987816601496, 1894.1817446392695),
          'Device002': (1494.017681525595, 1686.6241993763826),
          'Device001': (1139.949390830075, 1497.0908723196349),
          'Device005': (1260.050609169925, 702.9091276803653),
          'Device004': (1200.0, 1100.0)}
          
#layout = {'Device003': (1200.0, 1100.0),
#          'Device002': (1057.1060877331174, 1951.735312685243),
#          'Device001': (1128.5530438665587, 1525.8676563426216),
#          'Device004': (1271.4469561334413, 674.1323436573786)}

number_of_devices = len(layout)
array_output = [ 1.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.]
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
network_configuration = ['Star']
export_voltage = None
export_cables = None
ac_power_flow = True
target_burial_depth_array = 1.0
target_burial_depth_export = 2.0
connector_type = None
collection_point_type = None
devices_per_string = None
equipment_gradient_constraint = 14.0

# equipment soil installation compatibility matrix is loaded from file
file_name = 'equipment_compatibility_matrix.xlsx'
xls_file = pd.ExcelFile(os.path.join(data_dir, file_name), encoding = 'utf-8')
sheet_names = xls_file.sheet_names
equipment_soil_compatibility = xls_file.parse(sheet_names[0])

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
                               equipment_soil_compatibility)

## Create object and run
Electrical = Electrical(site, array, export, options, database)
solution, _ = Electrical.run_module(plot=True)
solution.print_result()

### Some plotting of these
#x_plot = [x for x in range(int(np.ceil(len(all_solutions)/2.)))]*2
#y_plot = [0]*int(np.ceil(len(all_solutions)/2.))+[1]*int(np.ceil(len(all_solutions)/2.))
#
#f, axarr = plt.subplots(2, int(np.ceil(len(all_solutions)/2.)))
#
## get devices
#device_x, device_y = zip(*[(oec[0], oec[1]) for oec in layout.itervalues()])
#
#for network, x, y in zip(all_solutions, x_plot, y_plot):
#
#    axarr[y, x].set_title('Axis [' + str(y) + ',' + str(x) + ']')
#    axarr[y, x].plot(landing_point[0], landing_point[1], 's', mew=2, markersize=15, color = 'none')
#    axarr[y, x].margins(0.1)
#    axarr[y, x].plot(device_x, device_y, 'o', mew=2, markersize=15, color = 'none')
#    axarr[y, x].grid(b=True)
#    for cp in network.collection_points:
#        axarr[y, x].plot(cp.location[0],cp.location[1],'k+', mew=2, markersize=15, color = 'black')

#x = np.linspace(0, 2 * np.pi, 400)
#y = np.sin(x ** 2)
#f, axarr = plt.subplots(2, int(np.ceil(len(all_solutions)/2.)))
#axarr[0, 0].plot(x, y)
#axarr[0, 0].set_title('Axis [0,0]')
#axarr[0, 1].scatter(x, y)
#axarr[0, 1].set_title('Axis [0,1]')
#axarr[0, 2].scatter(x, y)
#axarr[0, 2].set_title('Axis [0,2]')
#axarr[1, 0].plot(x, y ** 2)
#axarr[1, 0].set_title('Axis [1,0]')
#axarr[1, 1].scatter(x, y ** 2)
#axarr[1, 1].set_title('Axis [1,1]')
#axarr[1, 2].scatter(x, y ** 2)
#axarr[1, 2].set_title('Axis [1,2]')

#for network in all_solutions:
#    for cp in network.collection_points:
#        plt.plot(cp.location[0],cp.location[1],'k+', mew=2, markersize=15, color = 'black')
#        plt.show()
##    print network.array_cable_length
##    print network.device_to_device

#for network in all_solutions:
#    print network.hierarchy
#    print
#    print network.b_o_m
#    print
#    print network.economics_data
#    print
##    print network
#    for cable in network.export_cables:
#        print cable
#    for cable in network.array_cables:
#        print cable
#    print
##    print 'The number of wet-mate connectors is: ' + str(len(network.wet_mate))
##    print 'The number of dry-mate connectors is: ' + str(len(network.dry_mate))
#    print
##    print network.set_economics_data()
##    print
