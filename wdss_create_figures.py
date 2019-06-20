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


#import case_data
import sys
sys.path.append('C:/data/scripts/resources')
from my_functions import get_shapefile, latlon_from_radar, figure_timestamp
from my_functions import calc_srv, calc_new_extent, calc_dlatlon_dt, create_process_file_list
from case_data import this_case 
import matplotlib as mpl
import matplotlib.pyplot as plt
from custom_cmaps import sw_cmap,vg_cmap,ref_cmap,azdv_cmap,v_cmap
import numpy as np
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import os
#from datetime import datetime

# ----------------------------------------------

case_date = this_case['date']
rda = this_case['rda']
cut_list = this_case['cutlist']
products = this_case['products']

try:
    storm_motion = this_case['storm_motion']
    storm_dir = storm_motion[0]
    storm_speed = storm_motion[1]
except:
    pass
    #storm_dir = 240
    #storm_speed = 30

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
#src_dir = os.path.join(case_dir,'stage')
src_dir = os.path.join(case_dir,'netcdf')
files = create_process_file_list(src_dir,products,cut_list)

image_dir = os.path.join(base_dst_dir,'images')
mosaic_dir = os.path.join(image_dir,'mosaic')  

try:
    os.makedirs(image_dir)
except:
    pass


base_gis_dir = 'C:/data/GIS'
shape_path = os.path.join(base_gis_dir,'counties_mi','counties_MI.shp')
COUNTIES_ST = get_shapefile(shape_path)

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
plts['Conv_Shear_Gradient'] = {'cmap':vg_cmap,'vmn':0.000,'vmx':0.015,'title':'Conv Shear Gradient','cbticks':[0,0.005,0.010,0.015],'cblabel':'s $\mathregular{^-}{^1}$'}

plts['DivShear_Storm'] = {'cmap':azdv_cmap,'vmn':-0.01,'vmx':0.01,'title':'DivShear','cbticks':[-0.010,-0.005,0,0.005,0.010],'cblabel':'s $\mathregular{^-}{^1}$'}
plts['DivShear_Neg'] = {'cmap':azdv_cmap,'vmn':0.000,'vmx':0.015,'title':'Convergent Shear','cbticks':[0,0.005,0.010,0.015],'cblabel':'s $\mathregular{^-}{^1}$'}

plts['AzShear_Storm'] = {'cmap':azdv_cmap,'vmn':-0.01,'vmx':0.01,'title':'AzShear','cbticks':[-0.010,-0.005,0,0.005,0.010],'cblabel':'s $\mathregular{^-}{^1}$'}
plts['AzShear_Pos'] = {'cmap':azdv_cmap,'vmn':-0.01,'vmx':0.01,'title':'Positive AzShear','cbticks':[-0.010,-0.005,0,0.005,0.010],'cblabel':'s $\mathregular{^-}{^1}$'}

plts['Velocity'] = {'cmap':v_cmap,'vmn':-100,'vmx':100,'title':'Velocity','cbticks':[-100,-80,-60,-40,-20,0,20,40,60,80,100],'cblabel':'kts'}
plts['SRV'] = {'cmap':v_cmap,'vmn':-100,'vmx':100,'title':'SRV','cbticks':[-100,-80,-60,-40,-20,0,20,40,60,80,100],'cblabel':'kts'}
plts['SpectrumWidth'] = {'cmap':sw_cmap,'vmn':0,'vmx':40,'title':'Spectrum Width','cbticks':[0,10,15,20,25,40],'cblabel':'kts'}
plts['ReflectivityQC'] = {'cmap':ref_cmap,'vmn':-30,'vmx':80,'title':'Reflectivity','cbticks':[0,15,30,50,60],'cblabel':'dBZ'}



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


#file_list = [os.path.join(dp, f) for dp, dn, fn in os.walk(os.path.expanduser(src_dir)) for f in fn]
#files = sorted(os.listdir(src_dir))
#files = files[0:5]

arDict = {}

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
        #print('Velocity!')
        da_vel = da
        vel_arr = da.to_masked_array(copy=True)
        arDict[dtype] = {'ar':vel_arr,'lat':lats,'lon':lons}
        veldone = True
        # providing lats/lons for SRV (Storm Relative Velocity) 
        srv_lats = lats
        srv_lons = lons

    elif dtype == 'SpectrumWidth':
        da = dnew2.SpectrumWidth
        #print('Spec!')
        sw_arr = da.to_masked_array(copy=True)
        arDict[dtype] = {'ar':sw_arr,'lat':lats,'lon':lons}
        swdone = True

    elif dtype == 'ReflectivityQC':
        da = dnew2.ReflectivityQC
        #print('Ref!')
        ref_da = da
        ref_arr = da.to_masked_array(copy=True)
        arDict[dtype] = {'ar':ref_arr,'lat':lats,'lon':lons}
        refdone = True

    elif dtype == 'AzShear_Storm':
        da = dnew2.AzShear_Storm
        #print('Az!')
        az_arr_tmp = da.to_masked_array(copy=True)
        az_fill = az_arr_tmp.filled()
        az_fill[az_fill<-1] = 0
        az_shape = np.shape(az_fill)
        arDict[dtype] = {'ar':az_fill,'lat':lats,'lon':lons}
        azdone = True

        azpos_tmp = da.to_masked_array(copy=True)  
        azpos_fill = azpos_tmp.filled()
        azpos_fill[azpos_fill<0] = 0
        azpos_lats = lats
        azpos_lons = lons
        arDict['AzShear_Pos'] = {'ar':azpos_fill,'lat':azpos_lats,'lon':azpos_lons}
        azposdone = True

        # providing lats/lons for Velocity Gradient and Conv Shear Gradient
        vg_lats = lats
        vg_lons = lons
        csg_lats = lats
        csg_lons = lons

    elif dtype == 'DivShear_Storm':
        da = dnew2.DivShear_Storm
        #print('Div!')
        dv_arr_tmp = da.to_masked_array(copy=True)
        dv_fill = dv_arr_tmp.filled()
        dv_fill[dv_fill<-1] = 0
        dv_shape = np.shape(dv_fill)
        dvneg_lats = lats
        dvneg_lons = lons
        arDict[dtype] = {'ar':dv_fill,'lat':dvneg_lats,'lon':dvneg_lons}
        divdone = True

        dvneg_tmp = da.to_masked_array(copy=True)
        dvneg_fill = dvneg_tmp.filled()
        dvneg_fill[dvneg_fill<-1] = 0
        dvneg_fill[dvneg_fill>0] = 0
        dvnegdone = True
        arDict['DivShear_Neg'] = {'ar':dvneg_fill,'lat':lats,'lon':lons}

    else:
        pass

    if azdone and divdone:
        if (dv_shape == az_shape):
            # Velocity Gradient equals square root of (divshear**2 + azshear**2)
            #ar_sq = np.square(dv_fill) + np.square(az_fill)
            vg_sq = np.square(dv_fill) + np.square(az_fill)
            vg_arr = np.sqrt(vg_sq)
            arDict['Velocity_Gradient_Storm'] = {'ar':vg_arr,'lat':vg_lats,'lon':vg_lons}
            vgdone = True

            # Conv Shear Gradient equals square root of (negative_divshear**2 + positive_azshear**2)
            #ar_sq = np.square(dv_fill) + np.square(az_fill)
            csg_sq = np.square(dvneg_fill) + np.square(azpos_fill)
            csg_arr = np.sqrt(csg_sq)
            arDict['Conv_Shear_Gradient'] = {'ar':csg_arr,'lat':csg_lats,'lon':csg_lons}
            csgdone = True

    if veldone:
        # Create SRV from V given an input storm motion
        srv_arr = calc_srv(da_vel,storm_dir,storm_speed)
        arDict['SRV'] = {'ar':srv_arr,'lat':srv_lats,'lon':srv_lons}
        srvdone = True
            
    test = ['AzShear_Storm','DivShear_Storm','Velocity_Gradient_Storm','ReflectivityQC','SRV','Conv_Shear_Gradient']
    #if (divdone and veldone and swdone and refdone and azdone and vgdone):
    if (divdone and azdone and vgdone and csgdone and refdone and veldone and srvdone):
        fig, axes = plt.subplots(2,3,figsize=(14,7),subplot_kw={'projection': ccrs.PlateCarree()})
        plt.suptitle(t_str + '\n' + rda + '  ' + cut_str + '  Degrees')
        font = {'weight' : 'normal',
                'size'   : 12}
        plt.titlesize : 20
        plt.labelsize : 8
        plt.tick_params(labelsize=8)
        mpl.rc('font', **font)

        extent,x_ticks,y_ticks = calc_new_extent(orig_time,orig_extent,this_time,dlon_dt,dlat_dt)

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
                title_test = ['AzShear','DivShear','Velocity Gradient','Spectrum Width','Conv Shear Gradient']
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
        csgdone = False

        azdone = False
        azposdone = False 

        divdone = False
        dvnegdone = False

        veldone = False
        srvdone = False
        swdone = False
        refdone = False
        refqcdone = False

        plt.close()
    else:
        pass
