# -*- coding: utf-8 -*-
"""
Stages created WDSS-II nedcdf files (from 'wdss_create_netcdfs.py') to a directory
and renames them to ensure they are sorted chronologically. This then makes it possible
to use 'wdss_create_figures.py' to render multi-panel figures while ensuring each panel 
is associated with the same time and elevation cut.

The netcdf file directory structure is:
    {top directory}\{product}\{elevation}\

Each directory contains files with this naming convention:
    YYYYMMDD-HHMMSS.netcdf

Full example:
    /data/radar/20190519/GRR/DivShear_Storm/00.90/20190519-224051.netcdf

Different products at identical elevations typically share identical filenames
because timestamps are the same.Reflectivity filenames represent the slightly
earlier times they're created.


author: thomas.turnage@noaa.gov
Last updated: 2 Jun 2019
"""

def remove_duplicates(src_list):
    """
    Converts a list (src_list) with possible duplicate values.
    Returns a sorted list with duplicate members removed.
    This is function is not a required part of the workflow.
    It's just nice to get a list of elevation cuts and scan times for reference.                        
    """
    dst_list = []
    for n in src_list:
        if n not in dst_list:
            dst_list.append(n)
    return sorted(dst_list)
    

from case_data import this_case 
import os
import re
import shutil

# after a listed path as shown in the example above is converted to a string,
# the regular expressions below match the elevation(s) and time stamp
# in the path string so matching files can be copied to the staging directory

regex_elevation = re.compile('[0-9]{1,2}[.][0-9]{1,2}') # match 00.90 from above example
regex_timestamp = re.compile('[0-9]+[\-][0-9]+')        # match 20190519-224051 from above example


case_date = this_case['date']
rda = this_case['rda']
cut_list = this_case['cutlist']
topDir = '/data/radar'
verb = ' --verbose'

casePath = os.path.join(topDir,case_date,rda)
ncPath = os.path.join(casePath,'netcdf')

srcDir = ncPath
dstDir = os.path.join(casePath,'stage')
os.system('rm -rf ' + dstDir)
os.system('mkdir -p ' + dstDir)
imagePath = os.path.join(topDir,'images',case_date,rda)
os.system('mkdir -p ' + imagePath)
 # the staging directory to later be accessed by wdss_image_process.py
#try:
#    shutil.rmtree(dstDir)
#except:
#    pass
#os.makedirs(dstDir)


product_list = ['AzShear_Storm','DivShear_Storm','ReflectivityQC','SpectrumWidth','Velocity']

# other products of interest... ['RhoHV','Zdr','PhiDP']

# Dictionary used to compose new filenames as they get copied to staging directory
productDict = {'AzShear_Storm':'AzShear','DivShear_Storm':'DivShear',
               'ReflectivityQC':'Refl','SpectrumWidth':'SpecWidth',
               'Velocity_Gradient_Storm':'VelGrad', 'Velocity':'Velocity'}

# typical list of radar cuts (i.e., elevation slices) ...
# ['00.50', '00.90', '01.30', '01.80', '02.40', '03.10', '04.00', '05.10', '06.40', '08.00', '10.00', '12.50', '15.60', '19.50']
#
# To know for sure which cuts are available, inspect elevation list after running script, i.e.:
# print(unique_radar_elevations)
#
# Typically we want just a subset of cuts for processing


radar_elevations = []
scan_times = []
new_file_list = []

#Builds a sorted list of full file paths
file_list = [os.path.join(dp, f) for dp, dn, fn in os.walk(os.path.expanduser(srcDir)) for f in fn]
sorted_path_list = sorted(file_list)

# Builds a list of all available radar elevation cuts
all_radar_elevations = []
for p in range(0,len(sorted_path_list)):
    full_path_string = str(sorted_path_list[p])
    el_search = re.search(regex_elevation, full_path_string)
    try:
        elevation_string = el_search.group(0)
        all_radar_elevations.append(elevation_string)
    except:
        pass
unique_radar_elevations = remove_duplicates(all_radar_elevations)

# As documented earlier, this section needs work - i.e., more robust matching
for path in range(0,len(sorted_path_list)):
    path_string = str(sorted_path_list[path])
    for product in product_list:
        if product in path_string and 'PRF' not in path_string:
            for cuts in cut_list:
                if cuts in path_string:
                    new_file_list.append(path_string)
                    m1 = re.search(regex_elevation, path_string)
                    m2 = re.search(regex_timestamp, path_string)
                    try:
                        elevation_string = m1.group(0)
                    except:
                        pass
                    try:                        
                        file_timestamp_string = m2.group(0)
                    except:
                        pass

                    radar_elevations.append(elevation_string)
                    scan_times.append(file_timestamp_string)
                    info = [product,elevation_string,file_timestamp_string]
                    # copy selected files to a staging directory with new, more descriptive names
                    # Important because original timestamp filenames are identical in different source directories
                    # and don't want to overwrite anything. Want to ensure timestamp is beginning of name to 
                    # help with sorting files accordingly
                    new_fname = file_timestamp_string + "_" + productDict[product] + "_" + elevation_string + ".netcdf"
                    os.path.join(dstDir,new_fname)
                    print(dstDir + ' ' + new_fname)
                    #if '-223' in new_fname:
                    shutil.copy2(path_string,os.path.join(dstDir,new_fname))

scans = remove_duplicates(scan_times)
elevations = remove_duplicates(radar_elevations)