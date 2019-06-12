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

Required:    
         cases : key      : arbitrary name for case
          date : string   : for file/directory naming conventions   
           rda : string   : also for file/directory naming conventions
        latmax : float    : north map plot extent
        latmin : float    : south map plot extent
        lonmin : float    : west map plot extent
        lonmax : float    : east map plot extent
       cutlist : list     : used by wdss_stage_files.py to know which cuts to copy
  
Optional:    
      eventloc : float tuple  : event lat/lon pair to plot as marker 
     eventloc2 : float tuple  : (event lat/lon pair to plot as marker 
  start_latlon : float tuple  : initial lat/lon coordinates of feature
    end_latlon : float tuple  : final lat/lon coordinates of feature
    start_time : string       : product time stamp associated with start_latlon
                              ::  example... '2017-07-20 00:07:34'
      end_time : string       : product time stamp associated with end_latlon
  storm_motion : float tuple  : storm direction in degrees, speed in knots

feature_follow : boolean      : whether figures should be rendered feature follow
     
author: thomas.turnage@noaa.gov
Last updated:
    10 June 2019 - removed lat/lon ticks since make_ticks function in create_figures.py does this now
    11 June 2019 - added 'storm_motion' as input for srv function in create_figures.py
------------------------------------------------
"""

cases = {}
cases['a'] = {'date':'20170919',
     'rda':'KMVX',
     'latmax':48.1,
     'latmin':46.05,
     'lonmin':-98.75,
     'lonmax':-96.15,
     'eventloc': (-96.53,47.5452),
     'eventloc2': (-96.5337,47.5401),     
     'start_latlon': (47.31,-96.70),
     'end_latlon': (47.39,-96.65),
     'start_time': '2017-07-20 00:07:34',
     'end_time': '2017-07-20 00:18:01',
     'feature_follow': True,
     'cutlist': ['00.50', '00.90', '01.30', '01.80', '02.40']
     }


cases['b'] = {'date':'20180827',
     'rda':'KMVX',
     'latmax':48.25,
     'latmin':46.750,
     'lonmin':-97.4,
     'lonmax':-95.5,
     'eventloc': (-96.12,47.44)
     }


cases['c'] = {'date':'20190519',
     'rda':'KGRR',
     'lonmin':-85.72,
     'lonmax':-85.17,
     'latmin':42.3,
     'latmax':42.75,
     'eventloc': (-85.17,42.53),
     'extent':[-85.75,-85.0,39.30],
     'start_latlon': (42.43,-85.44),
     'end_latlon': (42.51,-85.246),
     'start_time': '2019-05-19 22:03:39',
     'end_time': '2019-05-19 22:18:11',
     'storm_motion': (240,30),
     'feature_follow': True,
     'cutlist': ['00.50']
     }


cases['d'] = {'date':'20190527',
     'rda':'KGLD',
     'lonmin':-103.00,
     'lonmax':-101.8,
     'latmin':38.05,
     'latmax':39.0,     
     'eventloc': (-102.29,38.59),
     'cutlist': ['00.50', '00.90', '01.30', '01.80', '02.40']
     }

cases['e'] = {'date':'20080608',
     'rda':'KGRR',
     'lonmin':-87.45,
     'lonmax':-85.8,
     'latmin':43.05,
     'latmax':44.35,
     'start_latlon': (43.51,-86.05),
     'end_latlon': (43.05,-84.65),
     'start_time': '2008-06-08 17:47:51',
     'end_time': '2008-06-08 19:42:28',     
     'feature_follow': True,
     'cutlist': ['02.40']
     }

cases['f'] = {'date':'20190528',
     'rda':'KTWX',
     'lonmin':-96.26,
     'lonmax':-95.17,
     'latmax':39.15,
     'latmin':38.29,
     'eventloc': (-95.41,38.81),
     'eventloc2': (-94.93,39.06),
     'cutlist': ['00.50']
     }


cases['g'] = {'date':'20190601',
     'rda':'KGRR',
     'lonmin':-85.61,
     'lonmax':-84.90,
     'latmax':42.62,
     'latmin':42.06,
     'eventloc': (-85,19,42.31),
     #'eventloc2': (-94.93,39.06),
     'cutlist': ['00.50']
     }

			
cases['h'] = {'date':'20190601',
     'rda':'KGRRfirst',
     'lonmin':-88.0,
     'lonmax':-84.5,
      'latmax':44.5,
     'latmin':42.0,
     #'eventloc': (-95.41,38.81),
     #'eventloc2': (-94.93,39.06),
     'cutlist': ['00.50','08.00']
     }


#cases['f'] = {'date':'20190528',
#     'rda':'KTWX',
#     #'lonmin':-96.20,
#     'lonmin':-96.26,
#     #'lonmax':-94.45,
#     'lonmax':-95.11,
#     'latmax':39.8,
#     'latmin':38.40,
#     'eventloc': (-95.41,38.81),
#     'eventloc2': (-94.93,39.06),
#     'cutlist': ['00.50', '00.90', '01.30', '01.80', '02.40']
#     }


import numpy as np
this_case = cases['c']

