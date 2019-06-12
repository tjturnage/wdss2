# -*- coding: utf-8 -*-
"""
Creates maps from wdss-ii nedcdf files that can be saved as images
WDSS-II information : http://www.wdssii.org/

Assumption: You previously ran wdss_stage_files.py to prep the netcdf file

author: thomas.turnage@noaa.gov
   Last updated : 10 Jun 2019

   New features:
            srv : function to create Storm Relative Velocity (SRV)
     make_ticks : function to create lat/lon ticks for an extent
        windows : Boolean - toggles betweeen running on windows vs. linux

        simpler way to call state shapefiles to add to maps
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
               factor for terrain, earth's curvature, or standard beam refraction.
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
    if num_gates <= 334:
        gate_len = 1000.0 * factor
    else:
        gate_len = 250.0 * factor
    rng = np.arange(2125.0,(num_gates*gate_len + 2125.0),gate_len)
    g = Geod(ellps='clrk66')
    center_lat = np.ones([len(az),len(rng)])*dnew2.Latitude
    center_lon = np.ones([len(az),len(rng)])*dnew2.Longitude
    az2D = np.ones_like(center_lat)*az[:,None]
    rng2D = np.ones_like(center_lat)*np.transpose(rng[:,None])
    lat,lon,back=g.fwd(center_lon,center_lat,az2D,rng2D)
    return lat,lon,back

def srv(velocity_data_array,storm_dir,storm_speed):
    """
    Subtracts storm motion from velocity bin values. This is based on the cosine of the angle
    between the storm direction and a given array radial direction.
    
    Example...
        storm_dir,storm_speed : 250,30 
                 storm motion : from 250 degrees at 30 knots
             max added amount : 30 kts at array's 70 degree radial
        max subtracted amount : 30 kts at array's 250 degree radial             
    
    Parameters
    ----------
    velocity_data_array : xarray
                          velocity array that will be recalculated
              storm_dir : float
                          storm motion direction in compass degrees
            storm_speed : integer or float
                          storm speed in knots
                    
    Returns
    -------
    Nothing - just adds new array to arDict
                    
    """
    da_new_speed = velocity_data_array
    # Storm motion is given as a "from" direction, so have to flip
    # this 180 degrees (equal to "pi" radians) to be consistent with
    # radial "to" direction convention
    storm_dir = math.radians(storm_dir) - math.pi

    for a in range(0,len(da_vel.Azimuth)):
        angle = math.radians(da_vel.Azimuth.values[a]) - storm_dir
        factor = math.cos(angle) * storm_speed
        da_new_speed[a] = da_new_speed[a] - factor 

    new_speed_arr = da_new_speed.to_masked_array(copy=True)
    arDict['SRV'] = {'ar':new_speed_arr,'lat':srv_lats,'lon':srv_lons}
    return


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


def calc_dlatlon_dt(starting_coords,starting_time,ending_coords,ending_time):
    """
    One-time calculation of changes in latitude and longitude with respect to time
    to implement feature-following zoom. Radar plotting software should be used
    beforehand to track a feature's lat/lon position at two different scan times
    along with the two different scan times.
    
    This needs to be executed only once at the beginning since dlat_dt and dlon_dt should
    remain constant.
    
    
    Parameters
    ----------
    starting_coords : tuple containing two floats - (lat,lon)
                      Established with radar plotting software to determine
                      starting coordinates for the feature of interest.

      starting_time : string
                      format - 'yyyy-mm-dd HH:MM:SS' - example - '2018-06-01 22:15:55'

      ending_coords : tuple containing two floats - (lat,lon)
                      Established with radar plotting software to determine
                      starting coordinates for the feature of interest.

        ending_time : string
                      format - 'yyyy-mm-dd HH:MM:SS' - example - '2018-06-01 23:10:05'
                                          
    Returns
    -------
            dlat_dt : float
                      Feature's movement in degrees latitude per second
            dlon_dt : float
                      Feature's movement in degrees longitude per second
                    
    """ 
    starting_datetime = datetime.strptime(starting_time, "%Y-%m-%d %H:%M:%S")
    ending_datetime = datetime.strptime(ending_time, "%Y-%m-%d %H:%M:%S")
    dt =  ending_datetime - starting_datetime
    dt_seconds = dt.seconds
    print('dt_seconds --------- ' + str(dt_seconds))
    dlat_dt = (ending_coords[0] - starting_coords[0])/dt_seconds
    dlon_dt = (ending_coords[1] - starting_coords[1])/dt_seconds

    return dlat_dt, dlon_dt
    

def make_ticks(this_min,this_max):
    """
    Determines range of tick marks to plot based on a provided range of degree coordinates.
    This function is typically called by 'calc_new_extent' after that function has calculated
    new lat/lon extents for feature-following zoom.
    
    Parameters
    ----------
           this_min : float
                      minimum value of either a lat or lon extent
           this_max : float
                      maximum value of either a lat or lon extent
                                          
    Returns
    -------
           tick_arr : float list
                      list of tick mark labels to use for plotting in the new extent

    """

    t_min = round(this_min) - 0.5
    t_max = round(this_max) + 0.5
    
    t_init = np.arange(t_min,t_max,0.5)    
    tick_arr = []
    for t in range(0,len(t_init)):
        if t_init[t] >= this_min and t_init[t] <= this_max:
            tick_arr.append(t_init[t])
        else:
            pass
    return tick_arr


def calc_new_extent(orig_t,t,lon_rate,lat_rate):
    """

    Shifts the map domain with each image to emulate AWIPS feature following zoom.
    Requires orig_extent [xmin,xmax,ymin,ymax] as a baseline to calculate new_extent

    
    Parameters
    ----------
             orig_t : integer
                      Epoch time (seconds) of first radar product in loop

                  t : integer
                      Epoch time (seconds) of current radar product in loop

           lon_rate : float
                      Change of longitude per second wrt time
                      determined ahead of time by 'calc_dlatlon_dt' function

           lat_rate : float
                      Change of latitude per second wrt time
                      determined ahead of time by 'calc_dlatlon_dt' function
                                          
    Returns
    -------
         new_extent : list of floats
                      [min longitude, max longitude, min latitude, max latitude]

             xticks : list of floats
                      list of longitude ticks to plot

             xticks : list of floats
                      list of latitude ticks to plot
                    
    """    
    time_shift = t - orig_t
    lon_shift = time_shift * lon_rate
    lat_shift = time_shift * lat_rate

    new_xmin = xmin + lon_shift
    new_xmax = xmax + lon_shift
    new_ymin = ymin + lat_shift
    new_ymax = ymax + lat_shift
    new_extent = [new_xmin,new_xmax,new_ymin,new_ymax]
    #call 'make_ticks' function to create ticks for new extent
    x_ticks = make_ticks(new_xmin,new_xmax)
    y_ticks = make_ticks(new_ymin,new_ymax)

    return new_extent,x_ticks,y_ticks

#import case_data
from case_data import this_case 
import matplotlib as mpl
import matplotlib.pyplot as plt
from custom_cmaps import sw_cmap,vg_cmap,ref_cmap,azdv_cmap,v_cmap
from pyproj import Geod
import numpy as np
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import os
import math
from datetime import datetime,timezone

# ----------------------------------------------

case_date = this_case['date']
rda = this_case['rda']
cut_list = this_case['cutlist']

try:
    storm_motion = this_case['storm_motion']
    storm_dir = storm_motion[0]
    storm_speed = storm_motion[1]
except:
    storm_dir = 240
    storm_speed = 30

# Running on Windows?
windows = True

if windows:
    topDir = 'C:/data'
    base_gis_dir = 'C:/data/GIS'
    base_dst_dir = os.path.join(topDir,case_date,rda)
else:
    topDir = '/data/radar'
    base_gis_dir = '/data/GIS'
    base_dst_dir = '/var/www/html/radar'

# case_dir example - C:/data/20190529/KGRR
case_dir = os.path.join(topDir,case_date,rda)
src_dir = os.path.join(case_dir,'stage')
image_dir = os.path.join(base_dst_dir,'images')
mosaic_dir = os.path.join(image_dir,'mosaic')  

try:
    os.makedirs(image_dir)
except:
    pass
#counties_list 

for ST in ['MI']:
    st = ST.lower()
    counties_dir = 'counties_' + st
    county_reader = 'county_' + st
    counties_shape = 'counties_' + ST + '.shp'
    COUNTIES_ST = 'COUNTIES_' + ST
    counties_shape_path = os.path.join(base_gis_dir,counties_dir,counties_shape)
    county_reader = counties_shape_path
    reader = shpreader.Reader(counties_shape_path)
    counties = list(reader.geometries())
    COUNTIES_ST = cfeature.ShapelyFeature(counties, ccrs.PlateCarree())

#county_mn = '/data/GIS/counties_mn/counties_MN.shp'
#reader_mn = shpreader.Reader(county_mn)
#counties_mn = list(reader_mn.geometries())
#COUNTIES_MN = cfeature.ShapelyFeature(counties_mn, ccrs.PlateCarree())

states = cfeature.NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lines',
        scale='50m',
        facecolor='none')

lon_formatter = LongitudeFormatter(number_format='.1f',
                       degree_symbol='',
                       dateline_direction_label=False,
                       zero_direction_label=False)

lat_formatter = LatitudeFormatter(number_format='.1f',
                      degree_symbol='')

# Grabs data from imported 'this_case'
ymin = this_case['latmin']
ymax = this_case['latmax']
xmin = this_case['lonmin']
xmax = this_case['lonmax']

orig_extent = [xmin,xmax,ymin,ymax]

case_date = this_case['date']
rda = this_case['rda']


# initialize these with False/None until all product arrays are made for
# the first iteration of image creation
vgdone = False
divdone = False
veldone = False
srvdone = False
swdone = False
refdone = False
azdone = False
azsq_done = False
divsq_done = False
dv_fixed = False
az_fixed = False
az_shape = None
dv_shape = None
got_orig_time = False

# plts is a dictionary with plotting instructions for each product
plts = {}
plts['Velocity_Gradient_Storm'] = {'cmap':vg_cmap,'vmn':0.000,'vmx':0.015,'title':'Velocity Gradient','cbticks':[0,0.005,0.010,0.015],'cblabel':'s $\mathregular{^-}{^1}$'}
plts['DivShear_Storm'] = {'cmap':azdv_cmap,'vmn':-0.01,'vmx':0.01,'title':'DivShear','cbticks':[-0.010,-0.005,0,0.005,0.010],'cblabel':'s $\mathregular{^-}{^1}$'}
plts['Velocity'] = {'cmap':v_cmap,'vmn':-100,'vmx':100,'title':'Velocity','cbticks':[-100,-80,-60,-40,-20,0,20,40,60,80,100],'cblabel':'kts'}
plts['SRV'] = {'cmap':v_cmap,'vmn':-100,'vmx':100,'title':'SRV','cbticks':[-100,-80,-60,-40,-20,0,20,40,60,80,100],'cblabel':'kts'}
plts['SpectrumWidth'] = {'cmap':sw_cmap,'vmn':0,'vmx':40,'title':'Spectrum Width','cbticks':[0,10,15,20,25,40],'cblabel':'kts'}
plts['ReflectivityQC'] = {'cmap':ref_cmap,'vmn':-30,'vmx':80,'title':'Reflectivity','cbticks':[0,15,30,50,60],'cblabel':'dBZ'}
plts['AzShear_Storm'] = {'cmap':azdv_cmap,'vmn':-0.01,'vmx':0.01,'title':'AzShear','cbticks':[-0.010,-0.005,0,0.005,0.010],'cblabel':'s $\mathregular{^-}{^1}$'}
#print(plts)

# list of products to be added to figure. I admit 'test' is a horribly ambiguous name
#test = ['AzShear_Storm','DivShear_Storm','Velocity_Gradient_Storm',
#        'ReflectivityQC','Velocity','SpectrumWidth']
test = ['AzShear_Storm','DivShear_Storm','Velocity_Gradient_Storm',
        'ReflectivityQC','SRV','SpectrumWidth']

arDict = {}

try:
    feature_following = this_case['feature_follow']
except:
    feature_following = False

if feature_following:
    start_latlon = this_case['start_latlon']
    end_latlon = this_case['end_latlon']
    start_time = this_case['start_time']
    end_time = this_case['end_time']
    dlat_dt,dlon_dt = calc_dlatlon_dt(start_latlon,start_time,end_latlon,end_time)
else:
    dlat_dt = 0
    dlon_dt = 0


file_list = [os.path.join(dp, f) for dp, dn, fn in os.walk(os.path.expanduser(src_dir)) for f in fn]
files = sorted(os.listdir(src_dir))
#files = files[0:5]

for filename in files:
    next_file = os.path.join(src_dir,filename)
    data = None
    num_gates = 0
    data = xr.open_dataset(next_file)
    dtype = data.TypeName
    degrees_tilt = data.Elevation
    cut_str = str(round(degrees_tilt,2))
    newcut = cut_str.replace(".", "")
    if degrees_tilt < 10:
        newcut = '0' + newcut
    t_str, img_fname_tstr = radar_cut_time(data.Time)
    this_time = data.Time

    if got_orig_time is False:
        orig_time = data.Time
        got_orig_time = True

    image_fname = img_fname_tstr + '_' + dtype + '_' + newcut + '.png'
    print('working on ' + image_fname[:-4])
    dnew2 = data.sortby('Azimuth')
    azimuths = dnew2.Azimuth.values
    num_gates = len(dnew2.Gate)
    lats,lons,back=latlon_from_radar(azimuths,degrees_tilt,num_gates)

    if dtype == 'Velocity':
        da = dnew2.Velocity * 1.944
        da_vel = da
        vel_arr = da.to_masked_array(copy=True)
        arDict[dtype] = {'ar':vel_arr,'lat':lats,'lon':lons}
        veldone = True
        # providing lats/lons for SRV (Storm Relative Velocity) 
        srv_lats = lats
        srv_lons = lons

    elif dtype == 'SpectrumWidth':
        da = dnew2.SpectrumWidth
        sw_arr = da.to_masked_array(copy=True)
        arDict[dtype] = {'ar':sw_arr,'lat':lats,'lon':lons}
        swdone = True

    elif dtype == 'ReflectivityQC':
        full_da = dnew2
        da = dnew2.ReflectivityQC
        ref_da = da
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
        azdone = True
        # providing lats/lons for Velocity Gradient 
        vg_lats = lats
        vg_lons = lons


    elif dtype == 'DivShear_Storm':
        da = dnew2.DivShear_Storm
        dv_arr_tmp = da.to_masked_array(copy=True)
        dv_fill = dv_arr_tmp.filled()
        dv_fill[dv_fill<-1] = 0
        dv_shape = np.shape(dv_fill)
        arDict[dtype] = {'ar':dv_fill,'lat':lats,'lon':lons}
        divdone = True

    else:
        pass

    if azdone and divdone:
        if (dv_shape == az_shape):
            # Velocity Gradient equals square root of (divshear**2 + azshear**2)
            ar_sq = np.square(dv_fill) + np.square(az_fill)
            vg_arr = np.sqrt(ar_sq)
            arDict['Velocity_Gradient_Storm'] = {'ar':vg_arr,'lat':vg_lats,'lon':vg_lons}
            vgdone = True

    if veldone:
        # Create SRV from V given an input storm motion
        da_new_speed = srv(da_vel,storm_dir,storm_speed)
        srvdone = True
            
    #if (divdone and veldone and swdone and refdone and azdone and vgdone):
    if (divdone and azdone and vgdone and refdone and swdone and veldone and srvdone):
        fig, axes = plt.subplots(2,3,figsize=(14,7),subplot_kw={'projection': ccrs.PlateCarree()})
        plt.suptitle(t_str + '\n' + rda + '  ' + cut_str + '  Degrees')
        font = {'weight' : 'normal',
                'size'   : 12}
        plt.titlesize : 20
        plt.labelsize : 8
        plt.tick_params(labelsize=8)
        mpl.rc('font', **font)

        extent,x_ticks,y_ticks = calc_new_extent(orig_time,this_time,dlon_dt,dlat_dt)

        for y,a in zip(test,axes.ravel()):
                this_title = plts[y]['title']
                a.set_extent(extent, crs=ccrs.PlateCarree())
                a.add_feature(COUNTIES_ST, facecolor='none', edgecolor='gray')
                a.tick_params(axis='both', labelsize=8)
                try:
                    a.plot(this_case['eventloc'][0], this_case['eventloc'][1], 'wv', markersize=3)
                except:
                    pass
                try:
                    a.plot(this_case['eventloc2'][0], this_case['eventloc2'][1], 'wv', markersize=3)
                except:
                    pass

                a.set_xticks(x_ticks, crs=ccrs.PlateCarree())
                a.set_yticks(y_ticks, crs=ccrs.PlateCarree())
                a.xaxis.set_major_formatter(lon_formatter)
                a.yaxis.set_major_formatter(lat_formatter)
 
                a.set_aspect(1.25)
                lon = arDict[y]['lon']
                lat = arDict[y]['lat']
                arr = arDict[y]['ar']
                title_test = ['AzShear','DivShear','Velocity Gradient','Spectrum Width']
                cs = a.pcolormesh(lat,lon,arr,cmap=plts[y]['cmap'],vmin=plts[y]['vmn'], vmax=plts[y]['vmx'])
                if this_title in title_test:
                    cax,kw = mpl.colorbar.make_axes(a,location='right',pad=0.05,shrink=0.9,format='%.4f')
                    cax.tick_params(labelsize=6)
                else:
                    cax,kw = mpl.colorbar.make_axes(a,location='right',pad=0.05,shrink=0.9)                    
                    cax.tick_params(labelsize=8)
                out=fig.colorbar(cs,cax=cax,**kw)
                label=out.set_label(plts[y]['cblabel'],size=10,verticalalignment='center')
                a.set_title(this_title)

        # name of figure file to be saved
        mosaic_fname = img_fname_tstr + '_' + newcut + '.png'

        # discovered it's much easier to organize by cuts
#        if not windows:
        mosaic_cut_dir = os.path.join(mosaic_dir,newcut)      
        try:
            os.makedirs(mosaic_cut_dir)
        except FileExistsError:
            pass


        image_dst_path = os.path.join(mosaic_cut_dir,mosaic_fname)
        plt.savefig(image_dst_path,format='png')
        print(mosaic_fname[:-4] + ' mosaic complete!')

        #reset these for the next image creation pass
        vgdone = False
        azdone = False
        divdone = False
        veldone = False
        srvdone = False
        swdone = False
        refdone = False
        refqcdone = False

        plt.close()
    else:
        pass
