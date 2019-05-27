# -*- coding: utf-8 -*-
"""
Stages wdss-ii nedcdf files to a directory and renames them
This makes it easier to subsequently use wdss_image_process.py to rended images
WDSS-II information : http://www.wdssii.org/

author: thomas.turnage@noaa.gov
Last updated: 26 May 2019
"""

def remove_duplicates(src_list):
    """
    Converts a list (src_list) with possible duplicate values
    Returns a sorted list with duplicate members removed 
                        
    """
    dst_list = []
    for n in src_list:
        if n not in dst_list:
            dst_list.append(n)
    return sorted(dst_list)
    

"""
The netcdf file directory structure is as follows:
{top directory}\{product}\{elevation}\

Each directory has files named with time stamps in this format:
YYYYMMDD-HHMMSS.netcdf

Example: 
{top directory}\DivShear_Storm\00.90\20190519-224051.netcdf

Different products at identical elevations typically share identical filenames because
timestamps are the same.
Reflectivity filenames represent the slightly earlier times in which they're created.
"""

import os
import re
import shutil

# after a listed path as shown in the example above is converted to a string,
# regular expressions can find the elevation and time stamp in that string
regex_elevation = re.compile('[0-9]{1,2}[.][0-9]{1,2}') # match 00.90 from above example
regex_timestamp = re.compile('[0-9]+[\-][0-9]+')        # match 20190519-224051 from above example



srcDir = "C:/data/wdss-ii/20190519_KGRR_netcdf"
dstDir = "C:/data/wdss-ii/stage099" # the staging directory to later be accessed by wdss_image_process.py
try:
    shutil.rmtree(dstDir)
except:
    pass
os.makedirs(dstDir)


# The following section needs work...
# For a given path string, sees if anything from product_list is in that string.
# Unfortunately, the way it's coded now, 'Velocity' also matches other products that begin with 'Velocity'.
# Similary,'Reflectivity' also returns 'ReflectivityPRF2'
# My temporary workaround is to change directory names...
# 'Velocity'     >> 'JustVelocity'
# 'Reflectivity' >> 'JustReflectivity'
# 
# Future work: develop more robust matching so directory names don't have to change
#
# product_list = ['AzShear_Storm','DivShear_Storm','JustReflectivity','SpectrumWidth','Velocity_Gradient_Storm','JustVelocity']
# other products of interest... ['RhoHV','Zdr','PhiDP']
product_list = ['JustVelocity']

# Dictionary used to compose new filenames as they get copied to staging directory
productDict = {'AzShear_Storm':'AzShear','DivShear_Storm':'DivShear','JustReflectivity':'Refl','SpectrumWidth':'SpecWidth','Velocity_Gradient_Storm':'VelGrad', 'JustVelocity':'Velocity'}

# typical list of radar cuts (i.e., elevation slices) ...
# ['00.50', '00.90', '01.30', '01.80', '02.40', '03.10', '04.00', '05.10', '06.40', '08.00', '10.00', '12.50', '15.60', '19.50']
#
# To know for sure which cuts are available, inspect elevation list after running script, i.e.:
# print(unique_radar_elevations)
#
# Typically we want just a subset of cuts for processing
cut_list = ['01.30']

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
        if product in path_string:
            for cuts in cut_list:
                if cuts in path_string:
                    new_file_list.append(path_string)
                    m1 = re.search(regex_elevation, path_string)
                    m2 = re.search(regex_timestamp, path_string)
                    elevation_string = m1.group(0)
                    file_timestamp_string = m2.group(0)
                    radar_elevations.append(elevation_string)
        
                    scan_times.append(file_timestamp_string)
                    info = [product,elevation_string,file_timestamp_string]
                    # copy selected files to a staging directory with new, more descriptive names
                    # Important because original timestamp filenames are identical in different source directories
                    # and don't want to overwrite anything. Want to ensure timestamp is beginning of name to 
                    # help with sorting files accordingly
                    new_fname = file_timestamp_string + "_" + productDict[product] + "_" + elevation_string + ".netcdf"
                    os.path.join(dstDir,new_fname)
                    print(new_fname)                
                    shutil.copy2(path_string,os.path.join(dstDir,new_fname))

scans = remove_duplicates(scan_times)
elevations = remove_duplicates(radar_elevations)

