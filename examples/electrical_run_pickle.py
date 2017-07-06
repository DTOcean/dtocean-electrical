# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 15:56:17 2016
"""
import os
import pickle

import pandas as pd
import numpy as np
from shapely.geometry import Polygon
import matplotlib.pyplot as plt

from dtocean_electrical import start_logging
from dtocean_electrical.main import Electrical
#from dtocean_electrical.hydro.hydro import emulate_wave_output


# Start the logging system
start_logging()

mod_path = os.path.realpath(__file__)
mod_dir = os.path.dirname(mod_path)

## Set test data directory
data_dir = os.path.join(mod_dir, "..", "sample_data")

# Pick up the pickled inputs
input_dict_file = os.path.join(mod_dir, "electrical_inputs.pkl")

with open(input_dict_file, "rb") as fstream:
    input_dict = pickle.load(fstream)

elec = Electrical(input_dict["site_data"],
                  input_dict["array_data"],
                  input_dict["export_data"],
                  input_dict["options"],
                  input_dict["database"])
 
## Create object and run
solution, installation_tool = elec.run_module()
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
