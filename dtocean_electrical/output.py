# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 15:56:41 2016

@author: acollin
"""
import matplotlib.pyplot as plt
from descartes import PolygonPatch

#import operator
#from shapely.geometry import Polygon, Point

def plot_devices(grid,
                 exclusion_lines,
                 layout,
                 landing_point,
                 footprint,
                 cp,
                 umbilical_cables,
                 array_cables,
                 export_cables):

    '''Show device locations, grid, exclusion zones and the landing point
    
    '''

    x_devices = []
    y_devices = []

    for key, item in layout.iteritems():
        x_devices.append(item[0])
        y_devices.append(item[1])

    fig = plt.figure()
    ax1 = fig.add_subplot(1,1,1, axisbg='lightskyblue')
    
    grid_plot = ax1.plot(grid.all_x,
                         grid.all_y,
                         'x',
                         mew=1,
                         markersize=1,
                         color='blue')
    
    ax1.plot(x_devices, y_devices, 'k+', mew=2, markersize=15, color = 'black')
    ax1.plot(landing_point[0],
             landing_point[1],
             'o',
             mew=2,
             markersize=15,
             color = 'black',
             fillstyle = 'none')
             
#    all_zs = [x.zorder for x in grid_plot]
#    print max(all_zs)

    for item in footprint:
        patch = PolygonPatch(item, zorder=3, fc='#cc00cc', ec='#555555')
        ax1.add_patch(patch)
        
    for line in exclusion_lines:
        plot_line(ax1, line, zorder=4)

    if cp:
        
        ax1.plot(cp[0].utm_x,
                 cp[0].utm_y,
                 's',
                 mew=2,
                 markersize=15,
                 color = 'black',
                 fillstyle = 'none')

    # get export
    export_x = []
    export_y = []
    
    for export in export_cables:
        
        (x, y) = zip(*[(grid.all_x[point], grid.all_y[point])
                                                for point in export.route])
            
        export_x.append(x)
        export_y.append(y)
        
    if export_x and export_y:

        # and plot
        ax1.plot(export_x[0],
                 export_y[0],
                 c='grey',
                 linewidth=2,
                 linestyle = '-')

    # get array
    for cable in array_cables:
        
        (x, y) = zip(*[(grid.all_x[point], grid.all_y[point])
                                                for point in cable.route])
        ax1.plot(x, y)
#        coords = solution.cable_routes[solution.cable_routes.marker == cable.marker][['x', 'y']]
#        ax1.plot(coords.x.tolist(), coords.y.tolist())

    ax1.margins(0.2)
    plt.axis('equal')
    plt.show()
    
def plot_line(ax, ob, zorder=1):
    x, y = ob.xy
    ax.plot(x,
            y,
            color='#cc00cc',
            linewidth=3,
            solid_capstyle='round',
            zorder=zorder)
