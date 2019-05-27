# -*- coding: utf-8 -*-
"""
Creates maps from wdss-ii nedcdf files that can be saved as images
WDSS-II information : http://www.wdssii.org/

Assumption: You previously ran wdss_file_process.py to stage the netcdf files

author: thomas.turnage@noaa.gov
Last updated: 26 May 2019
------------------------------------------------

""" 

def latlon_from_radar(product,az,elevation):
    """
    Convert radar bin radial coordinates to lat/lon coordinates.
    Adapted from Brian Blaylock code
    
    Parameters
    ----------
     product : string (examples: 'Reflectivity', 'DivShear_Storm', etc.)
               Needed to determine if Reflecticty is the product because
               the number of Reflectivity range gates changes with elevation
          az : numpy array
               All the radials for that particular product and elevation
               Changes from 720 radials for super-res product cuts to 360 radials
   elevation : float
               The radar elevation slice in degrees. Needed to:
               - determine and assign the number of range gates for Reflectivity
               - calculate range gate length (gate_len) as projected on the ground
        
    Returns
    -------
         lat : array like
         lon : array like
        back : I have no idea what this is for. I don't use it.
                    
    """
    rng = None
    factor = math.cos(math.radians(elevation))
    gate_len = 250.0 * factor
    if product == 'Reflectivity':
        if elevation < 1.0:
            rng = np.arange(2125.0,(1832*gate_len + 2125.0),gate_len)
        elif elevation < 1.5:
            rng = np.arange(2125.0,(1712*gate_len + 2125.0),gate_len)
        elif elevation < 2.0:
            rng = np.arange(2125.0,(1540*gate_len + 2125.0),gate_len)
        else:
            rng = np.arange(2125.0,(1336*gate_len + 2125.0),gate_len)   

        g = Geod(ellps='clrk66')
        center_lat = np.ones([len(az),len(rng)])*dnew2.Latitude
        center_lon = np.ones([len(az),len(rng)])*dnew2.Longitude
        az2D = np.ones_like(center_lat)*az[:,None]
        rng2D = np.ones_like(center_lat)*np.transpose(rng[:,None])
        lat,lon,back=g.fwd(center_lon,center_lat,az2D,rng2D)
        return lat,lon,back
    else:
        rng = np.arange(2125.0,(1192*gate_len + 2125.0),gate_len)
        g = Geod(ellps='clrk66')
        center_lat = np.ones([len(az),len(rng)])*dnew2.Latitude
        center_lon = np.ones([len(az),len(rng)])*dnew2.Longitude
        az2D = np.ones_like(center_lat)*az[:,None]
        rng2D = np.ones_like(center_lat)*np.transpose(rng[:,None])
        lat,lon,back=g.fwd(center_lon,center_lat,az2D,rng2D)
        return lat,lon,back
    

def plot_set():
    """
    Returns nothing, just sets the matplotlib plot environment 
    for consistency among image panes.
    Includes font/axes settings and defines plot display domain

    Future work: take input xlim,ylim values for feature following zoom    
    """
    font = {'family' : 'normal',
            'weight' : 'bold',
            'size'   : 20}
    plt.titlesize : 24
    plt.labelsize : 20
    plt.grid()
    #plt.linewidth : 3
    #lines.markersize : 10
    #xtick.labelsize : 16
    #ytick.labelsize : 16
    mpl.rc('font', **font)
    plt.xlim(-85.5,-84.9)
    plt.ylim(42.2,42.8)
    #plt.axis('off')
    return

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

def radar_cut_time(t):
    """
    Creates user-friendly time strings from Unix Epoch Time:
    https://en.wikipedia.org/wiki/Unix_time
    
    Parameters
    ----------
         t : integer
        
    Returns
    -------
         fig_title_timestring    : 'DD Mon YYYY  -  HH:MM:SS UTC'  
         fig_filename_timestring : 'YYYYMMDD-HHMMSS'
                    
    """    
    newt = datetime.fromtimestamp(t, timezone.utc)
    fig_title_timestring = datetime.strftime(newt, "%d %b %Y  -  %H:%M:%S %Z")
    fig_filename_timestring = datetime.strftime(newt, "%Y%m%d-%H%M%S")
    #image_fname = img_fname_tstr + '_' + newcut + '.png'
    return fig_title_timestring,fig_filename_timestring

from pyproj import Geod
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import cartopy
import os
import math
import sys
from metpy.plots import colortables
from datetime import datetime,timezone
from matplotlib.colors import LinearSegmentedColormap

#-------- Begin creating custom color maps --------
sw_colors = [(0,0,0),(220,220,255),(180,180,240),(50,50,150),(1,1,0),(1,150,0),(1,0,0),(1,1,1)]
sw_position = [0, 1/40, 5/40, 0.25, 15/40, 0.5, 0.75, 1]
sw_cmap=make_cmap(sw_colors, position=sw_position,bit=True)
plt.register_cmap(cmap=sw_cmap)

vg_colors = [(0, 0, 0),(30,30,30),(60,60,60),(110,70,0),(200,0,0),(1,1,1)]
vg_position = [0, 1/15, 2/15, 7/15, 8/15, 1]
vg_cmap=make_cmap(vg_colors, position=vg_position,bit=True)
plt.register_cmap(cmap=vg_cmap)

ref_colors = [(0,0,0),(130,130,130),(95,189,207),(57,201,105),(57,201,105),(0,40,0),(9,94,9),(255,207,0),(255,207,0),(255,207,0),(255,133,0),(255,0,0),(89,0,0),(255,245,255),(225,11,227),(164,0,247),(99,0,214),(5,221,224),(58,103,181),(255,255,255)]
ref_position = [0, 45/110, 46/110, 50/110, 51/110, 65/110, 66/110, 70/110, 71/110, 80/110, 81/110, 90/110, 91/110, 100/110, 101/110, 105/110, 106/110, 107/110, 109/110, 1]
ref_cmap=make_cmap(ref_colors, position=ref_position,bit=True)
plt.register_cmap(cmap=ref_cmap)

#azdv_colors = [(0,0,0),(1,1,1),(1,0,0),(0.7,0,0),(0,0,0),(0,0,0.7),(0,0,1),(1,1,1),(0,0,0)]
azdv_colors = [(0,0,0),(1,1,1),(0,0,1),(0,0,0.7),(0,0,0),(0.7,0,0),(1,0,0),(1,1,1),(0,0,0)]
azdv_position = [0, 0.1, 0.3, 0.43, 0.5, 0.57, 0.7, 0.9, 1]
azdv_cmap=make_cmap(azdv_colors, position=azdv_position)
plt.register_cmap(cmap=azdv_cmap)

v_norm, v_cmap = colortables.get_with_range('NWS8bitVel', -40, 40)
v_cmap.set_under('k')
plt.register_cmap(cmap=v_cmap)
#-------- End creating custom color maps --------

# initialize these with False until all product arrays are made for
# the first iteration of image creation
vgdone = False
divdone = False
veldone = False
swdone = False
refdone = False
azdone = False

srcDir = 'C:/data/wdss-ii/stage099'
file_list = [os.path.join(dp, f) for dp, dn, fn in os.walk(os.path.expanduser(srcDir)) for f in fn]

# A very big assumption here is that the 'file_list' listing is in chronological order
# This seems to work when running the fileList script first, because that script
# creates new filenames beginning with their timestamps
for next_file in file_list:
    data = None
    myFile = next_file
    data = xr.open_dataset(myFile)
    dtype = data.TypeName
    print(dtype)
    degrees_tilt = data.Elevation
    cut_str = str(round(degrees_tilt,2))
    newcut = cut_str.replace(".", "")
    t_str, img_fname_tstr = radar_cut_time(data.Time)

    image_fname = img_fname_tstr + '_' + newcut + '.png'
    dnew2 = data.sortby('Azimuth')
    azimuths = dnew2.Azimuth.values
    lat,lon,back=latlon_from_radar(dtype,azimuths,degrees_tilt)

    if dtype == 'Velocity_Gradient_Storm':
        da = dnew2.Velocity_Gradient_Storm
        np_vg = da.to_masked_array(copy=True)
        vgdone = True
    elif dtype == 'DivShear_Storm':
        da = dnew2.DivShear_Storm
        np_div = da.to_masked_array(copy=True)
        divdone = True
    elif dtype == 'Velocity':
        da = dnew2.Velocity
        np_vel = da.to_masked_array(copy=True)
        veldone = True
    elif dtype == 'SpectrumWidth':
        da = dnew2.SpectrumWidth
        np_sw = da.to_masked_array(copy=True)
        swdone = True
    elif dtype == 'Reflectivity':
        da = dnew2.Reflectivity
        np_ref = da.to_masked_array(copy=True)
        refdone = True
    elif dtype == 'AzShear_Storm':
        da = dnew2.AzShear_Storm
        np_az = da.to_masked_array(copy=True)
        azdone = True
    refl_lat,refl_lon,refl_back=latlon_from_radar('Reflectivity',azimuths,degrees_tilt)

    # Check if arrays were created for all products to be plotted
    # This whole plotting thing needs to be much more robust
    # It's currently hard-wired to 6 panes and needs other stuff added
    # such as user defined # of panes, color bars, maps, ticks, etc.
    if (vgdone and divdone and veldone and swdone and refdone and azdone):
        plt.figure(figsize=(21,14))
        plt.suptitle(t_str + '\n' + cut_str + ' Degrees')
        #title_end = str(degree_cut) + ' Degrees ' + t_str
        
        plt.subplot(231)
        plot_set()
        plt.title('AzShear')
        plt.pcolormesh(lat,lon,np_az,cmap=azdv_cmap,vmin=-0.02, vmax=0.02)
        
        plt.subplot(232)
        plot_set()
        plt.title('DivShear')
        plt.pcolormesh(lat,lon,np_div,cmap=azdv_cmap,vmin=-0.02, vmax=0.02)
        
        plt.subplot(233)
        plot_set()
        plt.title('Velocity Gradient')
        plt.pcolormesh(lat,lon,np_vg,cmap=vg_cmap,vmin=0.001, vmax=0.025)
        
        
        plt.subplot(236)
        plot_set()
        plt.title('Spectrum Width')
        plt.pcolormesh(lat,lon,np_sw,cmap=sw_cmap,vmin=0,vmax=40)
        
        
        plt.subplot(235)
        plot_set()
        plt.title('Velocity')
        plt.pcolormesh(lat,lon,np_vel,cmap=v_cmap,vmin=-50.0,vmax=50.0)
        #plt.colorbar()
        
        plt.subplot(234)
        plot_set()
        plt.title('Reflectivity')
        plt.pcolormesh(refl_lat,refl_lon,np_ref,cmap=ref_cmap,vmin=-30,vmax=80)
        
        #plt.colorbar()
        image_dst_path = 'C:/data/wdss-ii/python-images/' + image_fname
        plt.savefig(image_dst_path,format='png')
        plt.show()
        #reset these for the next image creation pass
        vgdone = False
        divdone = False
        veldone = False
        swdone = False
        refdone = False
        azdone = False
    else:
        pass
