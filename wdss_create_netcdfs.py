# -*- coding: utf-8 -*-
"""
Passes commands to WDSS-II so netcdfs from radar data can be created.
These netcdfs are later used to create figures.

For information on general WDSS-II workflow and command syntax:

    http://www.wdssii.org/

author: thomas.turnage@noaa.gov
Last updated: 28 May 2019
------------------------------------------------
"""

import os
import sys

cmd1 = 'export PATH=/data/WDSS2/bin:$PATH'
cmd2 = 'export LD_LIBRARY_PATH=/data/WDSS2/lib:$LD_LIBRARY_PATH'
cmd3 = 'export W2_CONFIG_LOCATION=/data/w2config:/data/WDSS2/w2config'

os.system(cmd1)
os.system(cmd2)
os.system(cmd3)

windows = False

try:
    sys.path.append('/data/scripts/resources')
except:
    print('can not connect resources folder')

from case_data import this_case

case_date = this_case['date']
rda = this_case['rda']
topDir = '/data/radar'
verb = ' --verbose'

casePath = os.path.join(topDir,case_date,rda)
ncPath = os.path.join(casePath,'netcdf')
ncIndexPath  = os.path.join(ncPath,'code_index.xml')
comboPath = os.path.join(ncPath,'combinedIndex.xml')
rawPath = os.path.join(casePath,'raw')
nseIndexPath = os.path.join(topDir,case_date,'NSE','code_index.xml')
#terrainPath  = terrDir + "/" + radar + ".nc.gz"

# The first step is to make netcdf files from the radar files
# If you've already done this, then comment out the "os.system" line
# so the generated command (ldm_cmd) isn't run
#
# An example of a generated command:
#    ldm2netcdf -i /localdata/data/radar/20180827/KMVX/raw -o /localdata/data/radar/20180827/KMVX/netcdf -a -1 -p KMVX -s KMVX -H -L --verbose=severe
# Typing just ldm2netcdf at a command prompt provides syntax information

ldm2nc_io = 'ldm2netcdf -i ' + rawPath + ' -o ' + ncPath # 
ldm2nc_arg = ' -a -1 -p ' + rda + ' -s ' + rda + ' -H -L ' + verb
ldm2nc_cmd = ldm2nc_io + ldm2nc_arg
os.system(ldm2nc_cmd)

# makeIndex.pl /localdata/data/radar/20180827/KMVX/netcdf code_index.xml
make_index_cmd = 'makeIndex.pl ' + ncPath + ' code_index.xml'
os.system(make_index_cmd)

#replaceIndex -i "/localdata/data/radar/20180827/KMVX/netcdf/code_index.xml /localdata/data/radar/20180827/NSE/code_index.xml" -o /localdata/data/radar/20180827/KMVX/netcdf/combinedIndex.xml --verbose=severe
rep_index_cmd = 'replaceIndex -i ' + ncIndexPath + ' ' + nseIndexPath + ' -o ' + comboPath + verb
os.system(rep_index_cmd)


#dealiasVel -i /localdata/data/radar/20180827/KMVX/netcdf/combinedIndex.xml -o /localdata/data/radar/20180827/KMVX/netcdf -R KMVX --verbose=severe
dealias_cmd = 'dealiasVel -i ' + comboPath + ' -o ' + ncPath + ' -R ' + rda + verb
os.system(dealias_cmd)

#makeIndex.pl /localdata/data/radar/20180827/KMVX/netcdf code_index.xml
os.system(make_index_cmd)

#w2qcnn -i /localdata/data/radar/20180827/KMVX/netcdf/combinedIndex.xml -o /localdata/data/radar/20180827/KMVX/netcdf -s KMVX -R 0.25x0.5x460 -u -T NearRadarEnvironmentTable --outputProducts="ReflectivityQC" --verbose=severe
qc_io = 'w2qcnn -i '  + comboPath + ' -o ' + ncPath
qc_arg = ' -s ' + rda + ' -R 0.25x0.5x460 -u -T NearRadarEnvironmentTable --outputProducts="ReflectivityQC"' + verb
qc_cmd = qc_io + qc_arg
os.system(qc_cmd)


#makeIndex.pl /localdata/data/radar/20180827/KMVX/netcdf code_index.xml
os.system(make_index_cmd)


#w2circ -i /localdata/data/radar/20180827/KMVX/netcdf/code_index.xml -o /localdata/data/radar/20180827/KMVX/netcdf
# -g "KMVX /localdata/data/terrain" -v Velocity -m -tot -az -div -S -D -sr -t -w -z ReflectivityQC --verbose=severe
w2circ_io = 'w2circ -i ' + ncIndexPath + ' -o '+ ncPath
w2circ_arg = ' -v Velocity -az -div -S -D -sr -t -w -z ReflectivityQC' + verb
w2circ_cmd = w2circ_io + w2circ_arg
os.system(w2circ_cmd)

#makeIndex.pl /localdata/data/radar/20180827/KMVX/netcdf code_index.xml
os.system(make_index_cmd)

#replaceIndex -i "/localdata/data/radar/20180827/KMVX/netcdf/code_index.xml /localdata/data/radar/20180827/NSE/code_index.xml" -o /localdata/data/radar/20180827/KMVX/netcdf/combinedIndex.xml --verbose=severe
os.system(rep_index_cmd)

#rm -r /localdata/data/radar/20180827/KMVX/netcdf/code_index.fam