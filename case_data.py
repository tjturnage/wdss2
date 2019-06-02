# -*- coding: utf-8 -*-
"""
A collection of dictionaries. Each dictionary contains metadata for a case.
Creating a dictionary for a case is the first in the whole work flow. A dictionary
gets selected with:

    this_case = cases['[letter]']    


That dictionary then gets imported into:

    wdss_create_netcdfs.py          :    make netcdfs for that case
    wdss_stage_files.py             :    copy/rename netcdfs to a staging directory
    wdss_create_figures.py          :    create images using staged netcdfs              
    

A dictionary contains the following case information:
    
        cases       : key       :  arbitrary name for case
        date        : string    : for file/directory naming conventions   
        rda         : string    : also for file/directory naming conventions
        latmax      : float     : north map plot extent
        latmin      : float     : south map plot extent
        lonmin      : float     : west map plot extent
        lonmax      : float     : east map plot extent
        lat_ticks   : list      : tick marks for map (keep within lat/lon bounds!)
        lon_ticks   : list      : tick marks for map (keep within lat/lon bounds!)
        eventloc    : tuple     : (optional) event lat/lon pair to plot as marker 
        eventloc2   : tuple     : (optional) event lat/lon pair to plot as marker 
        cutlist     : list      : used by wdss_stage_files.py to know which cuts to copy
     


author: thomas.turnage@noaa.gov
Last updated: 2 June 2019
------------------------------------------------
"""



def create_ticks(low,high):
    """
    This function currently not used. Purpose is to auto-create a range of map
    ticks based on the lat/lon ranges (lat_ticks and lon_ticks),
    but for now am just creating those manually
    """
    a = int(low)
    b = int(high) + 0.5
    if a < 0:
        ticks = np.arange(a,b,0.5)
        t =  ticks[:-1]
        return t
    else:
        ticks = np.arange(a,b,0.5)
        t = ticks[1:]        
        return t     



cases = {}
cases['a'] = {'date':'20170919',
     'rda':'KMVX',
     'latmax':48.050,
     'latmin':46.790,
     'lat_ticks': [47.0,47.5,48],
     'lonmin':-97.320,
     'lonmax':-95.760,
     'lon_ticks': [-97,-96.5,-96],
     'eventloc': (-96.53,47.5452),
     'eventloc2': (-96.5337,47.5401),     
     'cutlist': ['00.50', '00.90', '01.30', '01.80', '02.40', '03.10']
     }

cases['b'] = {'date':'20180827',
     'rda':'KMVX',
     'latmax':48.25,
     'latmin':46.750,
     'lonmin':-97.4,
     'lonmax':-95.5,
     'lon_ticks': [ -97, -96.5, -96],
     'lat_ticks': [47, 47.5, 48],
     'eventloc': (-96.12,47.44)
     }


cases['c'] = {'date':'20190519',
     'rda':'KGRR',
     'lonmin':-85.75,
     'lonmax':-85.0,
     'lon_ticks': [-85.5,-85], 
     'latmin':42.00,
     'latmax':43.00,
     'lat_ticks': [42.0,42.5,43],     
     'eventloc': (-85.17,42.53)
     #'extent':[-85.75,-85.0,39.30]
     }

cases['d'] = {'date':'20190527',
     'rda':'KGLD',
     'lonmin':-103.00,
     'lonmax':-101.8,
     'lon_ticks': [-103,-102.5,-102], 
     'latmin':38.05,
     'latmax':39.0,
     'lat_ticks': [38.5,39],     
     'eventloc': (-102.29,38.59),
     'cutlist': ['00.50', '00.90', '01.30', '01.80', '02.40']
     #'extent':[-85.75,-85.0,39.30]
     }

cases['e'] = {'date':'20080608',
     'rda':'KGRR',
     'lonmin':-86.8,
     'lonmax':-84.10,
     'lon_ticks': [-86.5,-86.0,-85.5,-85,-84.5], 
     'latmin':42.0,
     'latmax':44.165,
     'lat_ticks': [42.0,42.5,43,43.5,44],     
     'cutlist': ['00.50', '00.90', '01.30', '01.80', '02.40']
     #'extent':[-85.75,-85.0,39.30]
     }

import numpy as np
this_case = cases['e']

