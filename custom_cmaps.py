# -*- coding: utf-8 -*-
"""
Creates custom cmaps for matplotlib
https://matplotlib.org/3.1.0/tutorials/colors/colormap-manipulation.html

Assumption: You'll import this into your master image creation script

author: thomas.turnage@noaa.gov
Last updated: 28 May 2019
------------------------------------------------

""" 

def make_cmap(colors, position=None, bit=False):
    """
    Creates colormaps (cmaps) for different products.
    
    Information on cmap with matplotlib
    https://matplotlib.org/3.1.0/tutorials/colors/colormap-manipulation.html
    
    Parameters
    ----------
       colors : list of tuples containing RGB values. Tuples must be either:
                - arithmetic (zero to one) - ex. (0.5, 1, 0.75)
                - 8-bit                    - ex. (127,256,192)
     position : ordered list of floats
                None: default, returns cmap with equally spaced colors
                If a list is provided, it must have:
                  - 0 at the beginning and 1 at the end
                  - values in ascending order
                  - a number of elements equal to the number of tuples in colors
          bit : boolean         
                False : default, assumes arithmetic tuple format
                True  : set to this if using 8-bit tuple format
    Returns
    -------
         cmap
                    
    """  
    import numpy as np
    bit_rgb = np.linspace(0,1,256)
    if position == None:
        position = np.linspace(0,1,len(colors))
    else:
        if len(position) != len(colors):
            sys.exit("position length must be the same as colors")
        elif position[0] != 0 or position[-1] != 1:
            sys.exit("position must start with 0 and end with 1")
    if bit:
        for i in range(len(colors)):
            colors[i] = (bit_rgb[colors[i][0]],
                         bit_rgb[colors[i][1]],
                         bit_rgb[colors[i][2]])
    cdict = {'red':[], 'green':[], 'blue':[]}
    for pos, color in zip(position, colors):
        cdict['red'].append((pos, color[0], color[0]))
        cdict['green'].append((pos, color[1], color[1]))
        cdict['blue'].append((pos, color[2], color[2]))

    cmap = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,256)
    return cmap

import matplotlib as mpl
import matplotlib.pyplot as plt
import sys
from metpy.plots import colortables
from matplotlib.colors import LinearSegmentedColormap

#-------- Begin creating custom color maps --------

#--- Spectrum Width
sw_colors = [(0,0,0),(220,220,255),(180,180,240),(50,50,150),(255,255,0),(255,150,0),(255,0,0),(255,255,255)]
sw_position = [0, 1/40, 5/40, 0.25, 15/40, 0.5, 0.75, 1]
sw_cmap=make_cmap(sw_colors, position=sw_position,bit=True)
plt.register_cmap(cmap=sw_cmap)

#--- Velocity Gradient
vg_colors = [(0, 0, 0),(30,30,30),(60,60,60),(110,70,0),(200,0,0),(1,1,1)]
vg_position = [0, 1/15, 2/15, 7/15, 8/15, 1]
vg_cmap=make_cmap(vg_colors, position=vg_position,bit=True)
plt.register_cmap(cmap=vg_cmap)

#--- Reflectivity
ref_colors = [(0,0,0),(130,130,130),(95,189,207),(57,201,105),(57,201,105),(0,40,0),(9,94,9),(255,207,0),(255,207,0),(255,207,0),(255,133,0),(255,0,0),(89,0,0),(255,245,255),(225,11,227),(164,0,247),(99,0,214),(5,221,224),(58,103,181),(255,255,255)]
ref_position = [0, 45/110, 46/110, 50/110, 51/110, 65/110, 66/110, 70/110, 71/110, 80/110, 81/110, 90/110, 91/110, 100/110, 101/110, 105/110, 106/110, 107/110, 109/110, 1]
ref_cmap=make_cmap(ref_colors, position=ref_position,bit=True)
plt.register_cmap(cmap=ref_cmap)

#--- Azimuthal Shear / Div Shear
#azdv_colors = [(0,0,0),(1,1,1),(1,0,0),(0.7,0,0),(0,0,0),(0,0,0.7),(0,0,1),(1,1,1),(0,0,0)]
azdv_colors = [(0,0,0),(1,1,1),(0,0,1),(0,0,0.7),(0,0,0),(0.7,0,0),(1,0,0),(1,1,1),(0,0,0)]
azdv_position = [0, 0.1, 0.3, 0.43, 0.5, 0.57, 0.7, 0.9, 1]
azdv_cmap=make_cmap(azdv_colors, position=azdv_position)
plt.register_cmap(cmap=azdv_cmap)

#--- Velocity - need to home grow this so I don't require metpy
v_norm, v_cmap = colortables.get_with_range('NWS8bitVel', -40, 40)
v_cmap.set_under('k')
plt.register_cmap(cmap=v_cmap)
#-------- End creating custom color maps --------

