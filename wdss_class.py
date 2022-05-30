import nexradaws
import pandas as pd
import sys
import os
from datetime import datetime, timezone

try:
    os.listdir('/usr')
    windows = False
    radar_dir = '/data/radar'
    sys.path.append('/data/scripts/resources')
    cmd1 = 'export PATH=/data/WDSS2/bin:$PATH'
    cmd2 = 'export LD_LIBRARY_PATH=/data/WDSS2/lib:$LD_LIBRARY_PATH'
    cmd3 = 'export W2_CONFIG_LOCATION=/data/w2config:/data/WDSS2/w2config'
    os.system(cmd1)
    os.system(cmd2)
    os.system(cmd3)

except:
    windows = True
    radar_dir = 'C:/data/events'
    sys.path.append('C:/data/scripts/resources')



from case_data import this_case

case_date = this_case['date']
rda = this_case['rda']
#radar_dir = '/data/radar'
verb = ' --verbose'

class WDSS():

    def __init__(self,case_date=this_case['date'],rda = this_case['rda']):
      
        self.case_date = case_date
        self.rda = rda
        self.casePath = os.path.join(radar_dir,self.case_date,self.rda)
        self.ncPath = os.path.join(self.casePath,'netcdf')
        self.ncIndexPath  = os.path.join(self.ncPath,'code_index.xml')
        self.comboPath = os.path.join(self.ncPath,'combinedIndex.xml')
        self.rawPath = os.path.join(self.casePath,'raw')
        self.nseIndexPath = os.path.join(radar_dir,self.case_date,'NSE','code_index.xml')

    def create_netcdfs(self):
        ldm2nc_io = 'ldm2netcdf -i ' + self.rawPath + ' -o ' + self.ncPath # 
        ldm2nc_arg = ' -a -1 -p ' + self.rda + ' -s ' + self.rda + ' -H -L ' + verb
        ldm2nc_cmd = ldm2nc_io + ldm2nc_arg
        os.system(ldm2nc_cmd)

        make_index_cmd = 'makeIndex.pl ' + self.ncPath + ' code_index.xml'
        os.system(make_index_cmd)

        rep_index_cmd = 'replaceIndex -i ' + self.ncIndexPath + ' ' + self.nseIndexPath + ' -o ' + self.comboPath + verb
        os.system(rep_index_cmd)

        dealias_cmd = 'dealiasVel -i ' + self.comboPath + ' -o ' + self.ncPath + ' -R ' + self.rda + verb
        os.system(dealias_cmd)
        os.system(make_index_cmd)

        qc_io = 'w2qcnn -i '  + self.comboPath + ' -o ' + self.ncPath
        qc_arg = ' -s ' + self.rda + ' -R 0.25x0.5x460 -u -T NearRadarEnvironmentTable --outputProducts="ReflectivityQC"' + verb
        qc_cmd = qc_io + qc_arg
        os.system(qc_cmd)
        os.system(make_index_cmd)

        w2circ_io = 'w2circ -i ' + self.ncIndexPath + ' -o '+ self.ncPath
        w2circ_arg = ' -v Velocity -az -div -S -D -sr -t -w -z ReflectivityQC' + verb
        w2circ_cmd = w2circ_io + w2circ_arg
        os.system(w2circ_cmd)

        os.system(make_index_cmd)
        os.system(rep_index_cmd)

