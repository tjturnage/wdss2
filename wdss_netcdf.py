# -*- coding: utf-8 -*-
"""
Created on Wed May 29 16:37:23 2019

@author: tjtur
"""

import os

case_date = '20190519'
rda = 'KGRR'
topDir = '/data/radar'
verb = ' --verbose'

casePath = os.path.join(topDir,case_date,rda)
ncPath = os.path.join(casePath,'netcdf')
ncIndexPath  = os.path.join(ncPath,'code_index.xml')
comboPath = os.path.join(ncPath,'combinedIndex.xml')
rawPath = os.path.join(casePath,'raw')
nseIndexPath = os.path.join(topDir,case_date,'NSE','code_index.xml')
#terrainPath  = terrDir + "/" + radar + ".nc.gz"


#ldm2netcdf
#ldm2netcdf -i /localdata/data/radar/20180827/KMVX/raw -o /localdata/data/radar/20180827/KMVX/netcdf -a -1 -p KMVX -s KMVX -H -L --verbose=severe
ldm2nc_io = 'ldm2netcdf -i ' + rawPath + ' -o ' + ncPath # 
ldm2nc_arg = ' -a -1 -p ' + rda + ' -s ' + rda + ' -H -L ' + verb
ldm2nc_cmd = ldm2nc_io + ldm2nc_arg
os.system(ldm2nc_cmd)

#makeIndex.pl /localdata/data/radar/20180827/KMVX/netcdf code_index.xml
make_index_cmd = 'makeIndex.pl ' + ncPath + ' code_index.xml'
os.system(make_index_cmd)

#replaceIndex -i "/localdata/data/radar/20180827/KMVX/netcdf/code_index.xml /localdata/data/radar/20180827/NSE/code_index.xml" -o /localdata/data/radar/20180827/KMVX/netcdf/combinedIndex.xml --verbose=severe
rep_index_cmd = 'replaceIndex -i ' + ncIndexPath + ' ' + nseIndexPath + ' -o ' + comboPath + verb
os.system(rep_index_cmd)


#dealiasVel -i /localdata/data/radar/20180827/KMVX/netcdf/combinedIndex.xml -o /localdata/data/radar/20180827/KMVX/netcdf -R KMVX --verbose=severe
dealias_cmd = 'dealiasVel -i ' + comboPath + ' -o ' + ncPath + ' -R ' + rda + verb
print('---------------- ' + dealias_cmd)
os.system(dealias_cmd)

#makeIndex.pl /localdata/data/radar/20180827/KMVX/netcdf code_index.xml
os.system(make_index_cmd)


#w2qcnn -i /localdata/data/radar/20180827/KMVX/netcdf/combinedIndex.xml -o /localdata/data/radar/20180827/KMVX/netcdf -s KMVX -R 0.25x0.5x460 -u -T NearRadarEnvironmentTable --outputProducts="ReflectivityQC" --verbose=severe
#makeIndex.pl /localdata/data/radar/20180827/KMVX/netcdf code_index.xml

#w2circ -i /localdata/data/radar/20180827/KMVX/netcdf/code_index.xml -o /localdata/data/radar/20180827/KMVX/netcdf
# -g "KMVX /localdata/data/terrain" -v Velocity -m -tot -az -div -S -D -sr -t -w -z ReflectivityQC --verbose=severe
w2circ_io = 'w2circ -i ' + ncIndexPath + ' -o '+ ncPath
w2circ_arg = ' -v Velocity -az -div -S -D -sr -t -w -z Reflectivity' + verb
w2circ_cmd = w2circ_io + w2circ_arg
os.system(w2circ_cmd)

#makeIndex.pl /localdata/data/radar/20180827/KMVX/netcdf code_index.xml
os.system(make_index_cmd)

#replaceIndex -i "/localdata/data/radar/20180827/KMVX/netcdf/code_index.xml /localdata/data/radar/20180827/NSE/code_index.xml" -o /localdata/data/radar/20180827/KMVX/netcdf/combinedIndex.xml --verbose=severe
os.system(rep_index_cmd)

#rm -r /localdata/data/radar/20180827/KMVX/netcdf/code_index.fam

