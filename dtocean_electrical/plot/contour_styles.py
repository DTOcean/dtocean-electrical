#!/usr/bin/env python
"""
Illustrate simple contour plotting, contours on an image with
a colorbar for the contours, and labelled contours.

See also contour_image.py.
"""

import matplotlib.pyplot as plt
import numpy as np

def plot_contour(data):
    
    '''
    
    '''
    
    x = np.arange(0, 4, 1)
    y = np.arange(0, 4, 1)
    X, Y = np.meshgrid(x, y)
    Z = np.array(data)
    
    plt.figure()
    levels = [0, 0.25, .49]
    CS = plt.contour(X, Y, Z, levels)
    plt.clabel(CS, inline=1, fontsize=10)
    plt.margins(0.2)
    plt.title('Simplest default with labels')
    plt.show()
