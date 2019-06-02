# -*- coding: utf-8 -*-
"""
Creates maps from wdss-ii nedcdf files that can be saved as images
WDSS-II information : http://www.wdssii.org/

Assumption: You previously ran wdss_stage_files.py to prep the netcdf file

author: thomas.turnage@noaa.gov
Last updated: 26 May 2019
------------------------------------------------

""" 

def latlon_from_radar(az,elevation,num_gates):
    """
    Convert radar bin radial coordinates to lat/lon coordinates.
    Adapted from Brian Blaylock code
    
    Parameters
    ----------
          az : numpy array
               All the radials for that particular product and elevation
               Changes from 720 radials for super-res product cuts to 360 radials
   elevation : float
               The radar elevation slice in degrees. Needed to calculate range 
               gate length (gate_len) as projected on the ground using simple
               trigonometry. This is a very crude approximation that doesn't
               factor for terrain, earth's curvature, or for standard beam refraction.
   num_gates : integer
               The number of gates in a radial, which varies with 
               elevation and radar product. That is why each product makes 
               an individual call to this function. 
                    
    Returns
    -------
         lat : array like
         lon : array like
        back : I have no idea what this is for. I don't use it.
                    
    """
    rng = None
    factor = math.cos(math.radians(elevation))
    gate_len = 250.0 * factor
    rng = np.arange(2125.0,(num_gates*gate_len + 2125.0),gate_len)
    g = Geod(ellps='clrk66')
    center_lat = np.ones([len(az),len(rng)])*dnew2.Latitude
    center_lon = np.ones([len(az),len(rng)])*dnew2.Longitude
    az2D = np.ones_like(center_lat)*az[:,None]
    rng2D = np.ones_like(center_lat)*np.transpose(rng[:,None])
    lat,lon,back=g.fwd(center_lon,center_lat,az2D,rng2D)
    return lat,lon,back


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
    return fig_title_timestring,fig_filename_timestring

def plot_set(plot_type):
    """
    Returns nothing, just sets the matplotlib plot environment 
    for consistency among image panes. Includes font/axes settings
    and defines plot display domain. Much of this could get imported
    with case_data.py in the future

    Future work: take input xlim,ylim values for feature following zoom    
    """
    if plot_type == 'single':
        font = {'weight' : 'normal',
                'size'   : 10}
        plt.titlesize : 16
        plt.labelsize : 8
    else:
        font = {'weight' : 'normal',
                'size'   : 12}
        plt.titlesize : 20
        plt.labelsize : 8

    plt.tick_params(labelsize=8)
    mpl.rc('font', **font)


    ymin = this_case['latmin']
    ymax = this_case['latmax']
    xmin = this_case['lonmin']
    xmax = this_case['lonmax']
  
    plt.xlim(xmin,xmax)
    plt.ylim(ymin,ymax)
    return

def npsq(ar1,ar2):
    """
    input: two masked np arrays
    takes the square root of the sum of the squares
    velocity gradient = sqrt(azshear**2 + divshear**)
    """

    ar_sq = np.square(ar1) + np.square(ar2)
    ar_final = np.sqrt(ar_sq)
    return ar_final

from case_data import this_case 
import matplotlib as mpl
import matplotlib.pyplot as plt
import cartopy
from cartopy.mpl.geoaxes import GeoAxes
from mpl_toolkits.axes_grid1 import AxesGrid
from custom_cmaps import sw_cmap,vg_cmap,ref_cmap,azdv_cmap,v_cmap
from pyproj import Geod
import numpy as np
import xarray as xr
#from cartopy import config
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import os
import math
from datetime import datetime,timezone

county_mi = '/data/GIS/counties/counties_mi/counties_MI.shp'
reader = shpreader.Reader(county_mi)
counties = list(reader.geometries())
COUNTIES_MI = cfeature.ShapelyFeature(counties, ccrs.PlateCarree())

states = cfeature.NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lines',
        scale='50m',
        facecolor='none')

lon_formatter = LongitudeFormatter(number_format='.1f',
                       degree_symbol='',
                       dateline_direction_label=True)
lat_formatter = LatitudeFormatter(number_format='.1f',
                      degree_symbol='')



# Grabs data from imported 'this_case'
ymin = this_case['latmin']
ymax = this_case['latmax']
xmin = this_case['lonmin']
xmax = this_case['lonmax']
case_date = this_case['date']
rda = this_case['rda']


# initialize these with False/None until all product arrays are made for
# the first iteration of image creation
vgdone = False
divdone = False
veldone = False
swdone = False
refdone = False
azdone = False
azsq_done = False
divsq_done = False
dv_fixed = False
az_fixed = False
az_shape = None
dv_shape = None

# plts is a dictionary with plotting instructions for eah product
plts = {}
plts['Velocity_Gradient_Storm'] = {'cmap':vg_cmap,'vmn':0.000,'vmx':0.015,'title':'Velocity Gradient','cbticks':[0,0.01,0.02],'cblabel':'s $\mathregular{^-}{^1}$'}
plts['DivShear_Storm'] = {'cmap':azdv_cmap,'vmn':-0.01,'vmx':0.01,'title':'DivShear','cbticks':[-0.02,-0.01,0,0.01,0.02],'cblabel':'s $\mathregular{^-}{^1}$'}
plts['Velocity'] = {'cmap':v_cmap,'vmn':-50,'vmx':50,'title':'Velocity','cbticks':[-40,-30,-20,-10,0,10,20,30,40],'cblabel':'kts'}
plts['SpectrumWidth'] = {'cmap':sw_cmap,'vmn':0,'vmx':40,'title':'Spectrum Width','cbticks':[0,10,15,20,25,40],'cblabel':'kts'}
plts['ReflectivityQC'] = {'cmap':ref_cmap,'vmn':-30,'vmx':80,'title':'Reflectivity','cbticks':[0,15,30,50,60],'cblabel':'dBZ'}
plts['AzShear_Storm'] = {'cmap':azdv_cmap,'vmn':-0.01,'vmx':0.01,'title':'AzShear','cbticks':[-0.02,-0.01,0,0.01,0.02],'cblabel':'s $\mathregular{^-}{^1}$'}
#print(plts)

# list of products to be added to figure. I admit 'test' is a bad name
test = ['AzShear_Storm','DivShear_Storm','Velocity_Gradient_Storm','ReflectivityQC','Velocity','SpectrumWidth']

arDict = {}

topDir = '/data/radar'
htmlDir = '/var/www/html/radar/images'
casePath = os.path.join(topDir,case_date,rda)
ncPath = os.path.join(casePath,'netcdf')


srcDir = os.path.join(casePath,'stage')
#os.system('mkdir -p ' + dstDir)
imagePath = os.path.join(topDir,'images',case_date,rda)
os.system('mkdir -p ' + imagePath)


file_list = [os.path.join(dp, f) for dp, dn, fn in os.walk(os.path.expanduser(srcDir)) for f in fn]
files = sorted(os.listdir(srcDir))
#files = files[0:5]

 
for filename in files:
    next_file = os.path.join(srcDir,filename)
    data = None
    num_gates = 0
    myFile = next_file
    #myFile = filename
    data = xr.open_dataset(myFile)
    dtype = data.TypeName
    dstDir = os.path.join(htmlDir,case_date,rda,dtype)
    try:
        os.mkdirs(dstDir)
    except:
        pass
    degrees_tilt = data.Elevation
    cut_str = str(round(degrees_tilt,2))
    newcut = cut_str.replace(".", "")
    if degrees_tilt < 10:
        newcut = '0' + newcut
    t_str, img_fname_tstr = radar_cut_time(data.Time)

    image_fname = img_fname_tstr + '_' + dtype + '_' + newcut + '.png'
    print('working on ' + image_fname[:-4])
    dnew2 = data.sortby('Azimuth')
    azimuths = dnew2.Azimuth.values
    num_gates = len(dnew2.Gate)
    lats,lons,back=latlon_from_radar(azimuths,degrees_tilt,num_gates)

    if dtype == 'Velocity_Gradient_Storm':
        da= dnew2.Velocity_Gradient_Storm
        vg_arr = da.to_masked_array(copy=True)
        arDict[dtype] = {'ar':vg_arr,'lat':lats,'lon':lons}
        vgdone = True

    elif dtype == 'DivShear_Storm':
        da = dnew2.DivShear_Storm
        dv_arr_tmp = da.to_masked_array(copy=True)
        dv_fill = dv_arr_tmp.filled()
        dv_fill[dv_fill<-1] = 0
        dv_shape = np.shape(dv_fill)
        arDict[dtype] = {'ar':dv_fill,'lat':lats,'lon':lons}
        dv_fixed = True
        divdone = True

    elif dtype == 'Velocity':
        da = dnew2.Velocity
        vel_arr = da.to_masked_array(copy=True)
        arDict[dtype] = {'ar':vel_arr,'lat':lats,'lon':lons}
        veldone = True

    elif dtype == 'SpectrumWidth':
        da = dnew2.SpectrumWidth
        sw_arr = da.to_masked_array(copy=True)
        arDict[dtype] = {'ar':sw_arr,'lat':lats,'lon':lons}
        swdone = True

    elif dtype == 'ReflectivityQC':
        da = dnew2.ReflectivityQC
        ref_arr = da.to_masked_array(copy=True)
        arDict[dtype] = {'ar':ref_arr,'lat':lats,'lon':lons}
        refdone = True

    elif dtype == 'AzShear_Storm':
        da = dnew2.AzShear_Storm
        az_arr_tmp = da.to_masked_array(copy=True)
        az_fill = az_arr_tmp.filled()
        az_fill[az_fill<-1] = 0
        az_shape = np.shape(az_fill)
        arDict[dtype] = {'ar':az_fill,'lat':lats,'lon':lons}
        az_fixed = True
        azdone = True
        vg_lats = lats
        vg_lons = lons
        azdone = True
    else:
        pass

    if dv_fixed and az_fixed:
        if (dv_shape == az_shape):
            np_arr = npsq(dv_fill,az_fill)
            vg_arr = np_arr
            vgdone = True
            dtype = 'Velocity_Gradient_Storm'
            arDict[dtype] = {'ar':vg_arr,'lat':vg_lats,'lon':vg_lons}
            az_fixed = False
            dv_fixed = False
            az_shape = None
            dv_shape = None
            
    if (divdone and veldone and swdone and refdone and azdone and vgdone):
        fig, axes = plt.subplots(2,3,figsize=(17,10),subplot_kw={'projection': ccrs.PlateCarree()})
        #axes = axes.flatten()
        plt.suptitle(t_str + '\n' + rda + '  ' + cut_str + '  Degrees')

        extent = [this_case['lonmin'],this_case['lonmax'],this_case['latmin'],this_case['latmax']]
        for y,a in zip(test,axes.ravel()):
            a.set_extent(extent, crs=ccrs.PlateCarree())
            a.add_feature(COUNTIES_MI, facecolor='none', edgecolor='gray')
            #a.add_feature(COUNTIES_CO, facecolor='none', edgecolor='gray')
            try:
                a.plot(this_case['eventloc'][0], this_case['eventloc'][1], 'wv', markersize=3)
            except:
                pass
            try:
                a.plot(this_case['eventloc2'][0], this_case['eventloc2'][1], 'wv', markersize=3)
            except:
                pass
            a.set_xticks(this_case['lon_ticks'], crs=ccrs.PlateCarree())
            a.set_yticks(this_case['lat_ticks'], crs=ccrs.PlateCarree())
            a.xaxis.set_major_formatter(lon_formatter)
            a.yaxis.set_major_formatter(lat_formatter)
            a.set_aspect(1.45)
            lon = arDict[y]['lon']
            lat = arDict[y]['lat']
            arr = arDict[y]['ar']

            cs = a.pcolormesh(lat,lon,arr,cmap=plts[y]['cmap'],vmin=plts[y]['vmn'], vmax=plts[y]['vmx'])
            cax,kw = mpl.colorbar.make_axes(a,location='right',pad=0.05,shrink=0.9)
            out=fig.colorbar(cs,cax=cax,**kw)
            label=out.set_label(plts[y]['cblabel'],size=10,verticalalignment='center')
            a.set_title(plts[y]['title'])

        mosaic_fname = img_fname_tstr + '_' + newcut + '.png'
        mosaic_dir = os.path.join(htmlDir,case_date,rda,'mosaic')
        try:
            os.makedirs(mosaic_dir)
        except FileExistsError:
            pass
        image_dst_path2 = os.path.join(mosaic_dir,mosaic_fname)
        plt.savefig(image_dst_path2,format='png')
        print(image_fname[:-4] + ' mosaic complete!')

        #reset these for the next image creation pass
        vgdone = False
        divdone = False
        veldone = False
        swdone = False
        refdone = False
        refqcdone = False
        azdone = False
        plt.close()
    else:
        pass
