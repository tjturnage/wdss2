# -*- coding: utf-8 -*-
"""
Creates maps from wdss-ii nedcdf files that can be saved as images
WDSS-II information : http://www.wdssii.org/
author: thomas.turnage@noaa.gov
   Last updated : 23 Jun 2019
   New features:
            srv : function to create Storm Relative Velocity (SRV)
     make_ticks : function to create lat/lon ticks for an extent
        simpler way to call state shapefiles to add to maps
------------------------------------------------
""" 

import sys
import os

try:
    os.listdir('/var/www')
    windows = False
    sys.path.append('/data/scripts/resources')
    from case_data import this_case
    case_date = this_case['date']
    rda = this_case['rda']
    base_gis_dir = '/data/GIS'

    case_dir = os.path.join('/data/radar',case_date,rda)
    base_dst_dir = os.path.join('/var/www/html/radar/images',case_date,rda)
    mosaic_dir = os.path.join(base_dst_dir,'mosaic')
except:
    windows = True
    sys.path.append('C:/data/scripts/resources')
    from case_data import this_case
    base_gis_dir = 'C:/data/GIS'
    topDir = 'C:/data'
    case_date = this_case['date']
    rda = this_case['rda']
    case_dir = os.path.join(topDir,case_date,rda)
    base_dst_dir = os.path.join(topDir,'images',case_date,rda)
    mosaic_dir = os.path.join(base_dst_dir,'mosaic') 


cut_list = this_case['cutlist']
#products = this_case['products']
#products = ['ReflectivityQC','Velocity','RhoHV','AzShear_Storm']
products = ['ReflectivityQC','SRV','RhoHV']
mosaic_size = {}
mosaic_size[2] = {'h':6,'w':15,'rows':1,'columns':2}
mosaic_size[3] = {'h':6,'w':17,'rows':1,'columns':3}
mosaic_size[4] = {'h':12,'w':14,'rows':2,'columns':2}
mosaic_size[6] = {'h':12,'w':20,'rows':2,'columns':3}
height = mosaic_size[len(products)]['h']
width = mosaic_size[len(products)]['w']
rows = mosaic_size[len(products)]['rows']
cols = mosaic_size[len(products)]['columns']
if 'SRV' in products and 'Velocity' not in products:
    prod_temp = products.append('Velocity')
product_status = {}
for p in products:
    product_status[p] = 'no'

ymin = this_case['latmin']
ymax = this_case['latmax']
xmin = this_case['lonmin']
xmax = this_case['lonmax']
orig_extent = [xmin,xmax,ymin,ymax]

# test for existence of associated shapefiles
try:
    shapelist = this_case['shapelist']
except:
    pass

# test for case time range in which to create figures
# uses integer format of YYYYMMDDHHMMSS
# example: 20190720053000
    
try:
    start_fig = this_case['start_figures']
except:
    start_fig = 0
try:
    end_fig = this_case['end_figures']
except:
    end_fig = 99999999999999999

# test if storm motion exists in case data for SRV calculation
try:
    storm_motion = this_case['storm_motion']
    storm_dir = storm_motion[0]
    storm_speed = storm_motion[1]
except:
    storm_dir = 240
    storm_speed = 30


from my_functions import latlon_from_radar, figure_timestamp, build_html
from my_functions import calc_srv, calc_new_extent, calc_dlatlon_dt, create_process_file_list
from custom_cmaps import plts
from gis_layers import shape_mini
from re import search
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import cartopy.crs as ccrs
#import cartopy.feature as cfeature
#   from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker
#from datetime import datetime

# ----------------------------------------------

# create a list of absolute filepaths for the netcdf data to process  
src_dir = os.path.join(case_dir,'netcdf')
files = create_process_file_list(src_dir,products,cut_list,windows)

# ensure mosaic directory is created
try:
    os.makedirs(mosaic_dir)
except:
    pass


# initialize these with False/None until all product arrays are made for
# the first iteration of image creation
az_shape = None
dv_shape = None
got_orig_time = False

# if case data defines feature_following as true
# read in start/end times and lat and lon positions
# use them to call "calc_dlatlon_dt" function
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


arDict = {}

for filename in files:
    next_file = os.path.join(src_dir,filename)
    fsplit = str.split(filename,'/')
    ffile = fsplit[-1]
    ftimestr = ffile[0:8] + ffile[9:15]
    data = None
    num_gates = 0

    if (int(ftimestr) >= start_fig) and (int(ftimestr) < end_fig): 
        data = xr.open_dataset(next_file)
        dtype = data.TypeName
        degrees_tilt = data.Elevation
        cut_str = str(round(degrees_tilt,2))
        newcut = cut_str.replace(".", "")
        if degrees_tilt < 10:
            newcut = '0' + newcut
        t_str, img_fname_tstr = figure_timestamp(data.Time)
        this_time = data.Time


        if got_orig_time is False:
            orig_time = data.Time
            got_orig_time = True
    
        image_fname = img_fname_tstr + '_' + dtype + '_' + newcut + '.png'
        print('working on ' + image_fname[:-4])
        dnew2,lats,lons,back=latlon_from_radar(data)
    
        if dtype == 'Velocity':
            da = dnew2.Velocity * 1.944
            da_vel = da
            vel_arr = da.to_masked_array(copy=True)
            arDict[dtype] = {'ar':vel_arr,'lat':lats,'lon':lons}
            product_status[dtype] = 'yes'
            # providing lats/lons for SRV (Storm Relative Velocity) 
            srv_lats = lats
            srv_lons = lons
            srv_arr = calc_srv(da_vel,storm_dir,storm_speed)
            arDict['SRV'] = {'ar':srv_arr,'lat':srv_lats,'lon':srv_lons}
            product_status['SRV'] = 'yes'
    
        elif dtype == 'SpectrumWidth':
            print('skipping SW')
            da = dnew2.SpectrumWidth
            sw_arr = da.to_masked_array(copy=True)
            arDict[dtype] = {'ar':sw_arr,'lat':lats,'lon':lons}
            product_status[dtype] = 'yes'
            #pass
    
        elif dtype == 'ReflectivityQC':
            da = dnew2.ReflectivityQC
            ref_da = da
            ref_arr = da.to_masked_array(copy=True)
            arDict[dtype] = {'ar':ref_arr,'lat':lats,'lon':lons}
            product_status[dtype] = 'yes'

        elif dtype == 'RhoHV':
            da = dnew2.RhoHV
            cc_da = da
            cc_arr = da.to_masked_array(copy=True)
            arDict[dtype] = {'ar':cc_arr,'lat':lats,'lon':lons}
            product_status[dtype] = 'yes'
    
        elif dtype == 'AzShear_Storm':
            da = dnew2.AzShear_Storm
            az_arr_tmp = da.to_masked_array(copy=True)
            az_fill = az_arr_tmp.filled()
            az_fill[az_fill<-1] = 0
            az_shape = np.shape(az_fill)
            arDict[dtype] = {'ar':az_fill,'lat':lats,'lon':lons}
            product_status[dtype] = 'yes'
    
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
            product_status[dtype] = 'yes'
    
        else:
            pass
    

        if 'Velocity_Gradient_Storm' in products:
            if product_status['DivShear_Storm'] == 'yes' and product_status['AzShear_Storm'] == 'yes' and 'Velocity_Gradient_Storm' in products:
                if (dv_shape == az_shape):
                # Velocity Gradient equals square root of (divshear**2 + azshear**2)
                #ar_sq = np.square(dv_fill) + np.square(az_fill)
                    vg_sq = np.square(dv_fill) + np.square(az_fill)
                    vg_arr = np.sqrt(vg_sq)
                    arDict['Velocity_Gradient_Storm'] = {'ar':vg_arr,'lat':vg_lats,'lon':vg_lons}
                    product_status['Velocity_Gradient_Storm'] = 'yes'
    
        if 'no' not in product_status.values():
            fig, axes = plt.subplots(rows,cols,figsize=(width,height),subplot_kw={'projection': ccrs.PlateCarree()})
            font = {'weight' : 'normal',
                    'size'   : 24}
            plt.suptitle(t_str + '\n' + rda + '  ' + cut_str + '  Degrees')
            font = {'weight' : 'normal',
                    'size'   : 12}
            plt.titlesize : 24
            plt.labelsize : 8
            plt.tick_params(labelsize=8)
            mpl.rc('font', **font)
    
            extent,x_ticks,y_ticks = calc_new_extent(orig_time,orig_extent,this_time,dlon_dt,dlat_dt)
    
            for y,a in zip(products,axes.ravel()):
                    this_title = plts[y]['title']
                    a.set_extent(extent, crs=ccrs.PlateCarree())
                    a.tick_params(axis='both', labelsize=8)
                    for sh in shape_mini:
                        if search('survey', str(sh)):
                            a.add_feature(shape_mini[sh], facecolor='none', edgecolor='yellow', linewidth=0.7)
                        else:
                            a.add_feature(shape_mini[sh], facecolor='none', edgecolor='gray', linewidth=0.5)                            
                    try:
                        a.plot(this_case['eventloc'][0], this_case['eventloc'][1], 'wv', markersize=3)
                    except:
                        pass
                    try:
                        a.plot(this_case['eventloc2'][0], this_case['eventloc2'][1], 'wv', markersize=3)
                    except:
                        pass
                    
                    a.set_aspect(1.25)
                    a.xformatter = LONGITUDE_FORMATTER
                    a.yformatter = LATITUDE_FORMATTER
                    gl = a.gridlines(color='gray',alpha=0.0,draw_labels=True) 
                    gl.xlabels_top, gl.ylabels_right = False, False
                    gl.xlabel_style, gl.ylabel_style = {'fontsize': 9}, {'fontsize': 9}
                    gl.xlocator = mticker.FixedLocator(x_ticks)
                    gl.ylocator = mticker.FixedLocator(y_ticks)
    
                    lon = arDict[y]['lon']
                    lat = arDict[y]['lat']
                    arr = arDict[y]['ar']
                    title_test = ['AzShear','DivShear']
                    #title_test = ['AzShear','DivShear','Velocity Gradient','Spectrum Width','Conv Shear Gradient']
                    cs = a.pcolormesh(lat,lon,arr,cmap=plts[y]['cmap'],vmin=plts[y]['vmn'], vmax=plts[y]['vmx'])
                    if this_title in title_test:
                        cax,kw = mpl.colorbar.make_axes(a,location='right',pad=0.05,shrink=0.9,format='%.4f')
                        cax.tick_params(labelsize=6)
                    else:
                        cax,kw = mpl.colorbar.make_axes(a,location='right',pad=0.05,shrink=0.9)                    
                        cax.tick_params(labelsize=6)
                    out=fig.colorbar(cs,cax=cax,**kw)
                    label=out.set_label(plts[y]['cblabel'],size=8,verticalalignment='center')
                    a.set_title(this_title)
    
            # name of figure file to be saved
            mosaic_fname = img_fname_tstr + '_' + newcut + '.png'
            mosaic_cut_dir = os.path.join(mosaic_dir,newcut)      
            try:
                os.makedirs(mosaic_cut_dir)
            except FileExistsError:
                pass
    
            image_dst_path = os.path.join(mosaic_cut_dir,mosaic_fname)
            plt.savefig(image_dst_path,format='png')
            print(mosaic_fname[:-4] + ' mosaic complete!')
    
            #reset these for the next image creation pass
            for p in products:
                product_status[p] = 'no'

            plt.show()
            plt.close()
        else:
            pass
  
    for c in cut_list:
        new_c = c.replace(".", "")
        new_c = new_c[:-1]
        image_dir = os.path.join(mosaic_dir,new_c)
        try:
            build_html(image_dir)
        except:
            pass