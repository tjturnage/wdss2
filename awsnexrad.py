import nexradaws
import pandas as pd
import sys
import os
from datetime import datetime, timezone

try:
    os.listdir('/usr')
    windows = False
    base_dir = '/data/radar'
    sys.path.append('/data/scripts/resources')
except:
    windows = True
    base_dir = 'C:/data/events'
    sys.path.append('C:/data/scripts/resources')


from case_data import this_case

dest_dir = os.path.join(base_dir,this_case['date'],this_case['rda'],'raw')
print(dest_dir)
os.makedirs(dest_dir, exist_ok=True)

start = pd.Timestamp(this_case['start_time']).tz_localize(tz='UTC') - pd.Timedelta(minutes=2) 
end = pd.Timestamp(this_case['end_time']).tz_localize(tz='UTC') + pd.Timedelta(minutes=2) 

#print(start,end)

conn = nexradaws.NexradAwsInterface()
scans = conn.get_avail_scans_in_range(start, end, this_case['rda'])
print("There are {} scans available between {} and {}\n".format(len(scans), start, end))
print(scans[0:4])

# download these files
#results = conn.download(scans[0:2], templocation)
results = conn.download(scans, dest_dir)